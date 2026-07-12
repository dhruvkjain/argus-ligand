//! The Cloudflare Worker HTTP handler.
//!
//! This is the only part of the crate that talks to the Worker runtime, and it
//! is compiled only for the `wasm32` target. It stays thin: it serves the UI
//! and forwards `POST /scan` bodies to [`crate::scan`]. All of the biology
//! lives in the engine modules.

use crate::ai::{decompose, Decomposition};
use crate::engine::run_request;
use crate::scan;
use crate::types::ScanRequest;
use serde::Deserialize;
use worker::*;

/// Body of a `POST /ask` request: a natural language prompt and a sequence.
#[derive(Deserialize)]
struct AskBody {
    /// The user's plain English request, for example "find TATA boxes".
    prompt: String,
    /// The DNA sequence to scan. Not sent to the model.
    #[serde(default)]
    sequence: String,
}

/// Build a JSON error response with the given HTTP status.
///
/// # Arguments
///
/// * `message` - The error text to place under an `error` field.
/// * `status` - The HTTP status code to set.
///
/// # Returns
///
/// A JSON [`Response`] carrying `{"error": message}` at `status`.
fn json_error(message: &str, status: u16) -> Result<Response> {
    Ok(Response::from_json(&serde_json::json!({ "error": message }))?.with_status(status))
}

/// Look up a cached decomposition for a prompt key.
///
/// Any problem (missing binding, KV error, stale value that no longer parses) is
/// treated as a cache miss, so caching never breaks the request.
///
/// # Arguments
///
/// * `cache` - The KV store, or `None` if the binding is absent.
/// * `key` - The normalized prompt key.
///
/// # Returns
///
/// `Some(decomposition)` on a hit, `None` on any miss.
async fn cache_get(cache: &Option<worker::kv::KvStore>, key: &str) -> Option<Decomposition> {
    let store = cache.as_ref()?;
    let text = store.get(key).text().await.ok()??;
    serde_json::from_str(&text).ok()
}

/// Store a decomposition under a prompt key. Best effort; errors are ignored.
///
/// # Arguments
///
/// * `cache` - The KV store, or `None` if the binding is absent.
/// * `key` - The normalized prompt key.
/// * `decomposition` - The value to cache.
async fn cache_put(cache: &Option<worker::kv::KvStore>, key: &str, decomposition: &Decomposition) {
    if let Some(store) = cache {
        if let Ok(text) = serde_json::to_string(decomposition) {
            if let Ok(builder) = store.put(key, text) {
                let _ = builder.execute().await;
            }
        }
    }
}

/// How many real AI calls one visitor may make before the limit resets.
const AI_CALL_LIMIT: u32 = 5;

/// Identify the caller for rate limiting. Uses the Cloudflare client IP header,
/// falling back to a shared key in local development where it is absent.
fn client_ip(req: &Request) -> String {
    req.headers()
        .get("cf-connecting-ip")
        .ok()
        .flatten()
        .unwrap_or_else(|| "local".to_string())
}

/// Read how many AI calls this visitor has used in the current window.
///
/// # Arguments
///
/// * `store` - The KV store holding the counter.
/// * `ip` - The visitor key from [`client_ip`].
///
/// # Returns
///
/// The used count, or 0 on any miss or parse failure.
async fn ai_calls_used(store: &worker::kv::KvStore, ip: &str) -> u32 {
    store
        .get(&format!("rl:{ip}"))
        .text()
        .await
        .ok()
        .flatten()
        .and_then(|t| t.parse().ok())
        .unwrap_or(0)
}

/// Record a used AI call, refreshing the 24 hour expiry. Best effort.
///
/// The counter carries a rolling one day TTL, so a visitor gets a fresh budget
/// a day after their last call. Cache hits never call this, so they are free.
///
/// # Arguments
///
/// * `store` - The KV store holding the counter.
/// * `ip` - The visitor key from [`client_ip`].
/// * `new_count` - The count to store.
async fn bump_ai_calls(store: &worker::kv::KvStore, ip: &str, new_count: u32) {
    if let Ok(builder) = store.put(&format!("rl:{ip}"), new_count.to_string()) {
        let _ = builder.expiration_ttl(86_400).execute().await;
    }
}

/// Build the `/ask` JSON response from a decomposition.
///
/// Runs the engine when the model returned scans, or returns the clarifying
/// question. Every response carries the visitor's remaining AI-call budget.
///
/// # Arguments
///
/// * `decomposition` - The model's interpretation of the prompt.
/// * `cached` - Whether it came from the cache (no AI call spent).
/// * `remaining` - AI calls left in the visitor's window.
/// * `sequence` - The DNA sequence to scan.
///
/// # Returns
///
/// A JSON [`Response`] with either `clarification` or `interpreted` + `result`.
fn ask_response(
    decomposition: Decomposition,
    cached: bool,
    remaining: u32,
    sequence: String,
) -> Result<Response> {
    match decomposition {
        Decomposition::NeedsClarification(question) => Response::from_json(&serde_json::json!({
            "clarification": question,
            "cached": cached,
            "remaining": remaining,
            "max": AI_CALL_LIMIT,
        })),
        Decomposition::Scans(scans) => {
            let interpreted = serde_json::to_value(&scans).unwrap_or_default();
            let request = ScanRequest { sequence, scans };
            match run_request(request) {
                Ok(result) => Response::from_json(&serde_json::json!({
                    "interpreted": interpreted,
                    "result": result,
                    "cached": cached,
                    "remaining": remaining,
                    "max": AI_CALL_LIMIT,
                })),
                Err(e) => json_error(&e.to_string(), 400),
            }
        }
    }
}

/// The UI page and its stylesheet and script, embedded at compile time.
///
/// Serving these from the Worker itself means there is no separate static asset
/// service to configure.
const INDEX_HTML: &str = include_str!("../public/index.html");
const INDEX_CSS: &str = include_str!("../public/index.css");
const INDEX_JS: &str = include_str!("../public/index.js");

/// Serve static text with an explicit content type.
///
/// # Arguments
///
/// * `body` - The file contents.
/// * `content_type` - The MIME type to send.
///
/// # Returns
///
/// A [`Response`] carrying `body` with the given content type.
fn text_asset(body: &str, content_type: &str) -> Result<Response> {
    let mut resp = Response::ok(body)?;
    resp.headers_mut().set("content-type", content_type)?;
    Ok(resp)
}

/// Handle every incoming HTTP request.
///
/// Routes `GET /` to the embedded UI and `POST /scan` to the engine. Any other
/// path falls through to a 404 from the router.
///
/// # Arguments
///
/// * `req` - The incoming request.
/// * `env` - Worker bindings and environment. Unused in milestone 1.
/// * `_ctx` - The Worker execution context. Unused here.
///
/// # Returns
///
/// The HTTP response to send back, or a [`worker::Error`] if the runtime fails
/// to build one.
#[event(fetch)]
async fn fetch(req: Request, env: Env, _ctx: Context) -> Result<Response> {
    Router::new()
        .get("/", |_req, _ctx| Response::from_html(INDEX_HTML))
        .get("/index.css", |_req, _ctx| {
            text_asset(INDEX_CSS, "text/css; charset=utf-8")
        })
        .get("/index.js", |_req, _ctx| {
            text_asset(INDEX_JS, "application/javascript; charset=utf-8")
        })
        .get("/samples", |_req, _ctx| {
            Response::from_json(&crate::samples::list())
        })
        .get("/samples/:id", |_req, ctx| {
            match ctx.param("id").and_then(|id| crate::samples::fasta(id)) {
                Some(fasta) => {
                    let mut resp = Response::ok(fasta)?;
                    resp.headers_mut()
                        .set("content-type", "text/plain; charset=utf-8")?;
                    Ok(resp)
                }
                None => Response::error("unknown sample", 404),
            }
        })
        .get_async("/quota", |req, ctx| async move {
            // Read-only: report the caller's remaining AI budget without spending it,
            // so the page can show the true count right after a refresh.
            let ip = client_ip(&req);
            let cache = ctx.env.kv("PROMPT_CACHE").ok();
            let remaining = match &cache {
                Some(s) => AI_CALL_LIMIT.saturating_sub(ai_calls_used(s, &ip).await),
                None => AI_CALL_LIMIT,
            };
            Response::from_json(&serde_json::json!({ "remaining": remaining, "max": AI_CALL_LIMIT }))
        })
        .post_async("/scan", |mut req, _ctx| async move {
            let body = req.text().await?;
            let out = scan(&body); // JSON string in, JSON string out
            let mut resp = Response::ok(out)?;
            resp.headers_mut()
                .set("content-type", "application/json; charset=utf-8")?;
            Ok(resp)
        })
        .post_async("/ask", |mut req, ctx| async move {
            let ip = client_ip(&req);
            let body: AskBody = match req.json().await {
                Ok(body) => body,
                Err(_) => {
                    return json_error("request body must be JSON with a prompt", 400);
                }
            };
            if body.prompt.trim().is_empty() {
                return json_error("prompt is empty", 400);
            }

            let key = crate::ai::cache_key(&body.prompt);
            let cache = ctx.env.kv("PROMPT_CACHE").ok();

            // A cache hit makes no AI call, so it does not spend the daily budget.
            if let Some(hit) = cache_get(&cache, &key).await {
                let remaining = match &cache {
                    Some(s) => AI_CALL_LIMIT.saturating_sub(ai_calls_used(s, &ip).await),
                    None => AI_CALL_LIMIT,
                };
                return ask_response(hit, true, remaining, body.sequence);
            }

            // Cache miss means a real AI call, so enforce the per-visitor limit.
            let used = match &cache {
                Some(s) => ai_calls_used(s, &ip).await,
                None => 0,
            };
            if used >= AI_CALL_LIMIT {
                return Response::from_json(&serde_json::json!({
                    "limit_exceeded": true,
                    "remaining": 0,
                    "max": AI_CALL_LIMIT,
                    "message": format!(
                        "Daily limit of {AI_CALL_LIMIT} AI requests reached. Repeated questions are free (cached), and manual mode has no limit."
                    ),
                }));
            }
            if let Some(s) = &cache {
                bump_ai_calls(s, &ip, used + 1).await;
            }
            let remaining = AI_CALL_LIMIT.saturating_sub(used + 1);

            let fresh = match decompose(&ctx.env, &body.prompt).await {
                Ok(fresh) => fresh,
                Err(e) => return json_error(&format!("AI decomposition failed: {e}"), 502),
            };
            // Cache only definitive interpretations, never clarifications: a
            // clarification should be re-asked, and must not get stuck if the tool
            // later gains the capability the prompt needs.
            if matches!(fresh, Decomposition::Scans(_)) {
                cache_put(&cache, &key, &fresh).await;
            }
            ask_response(fresh, false, remaining, body.sequence)
        })
        .run(req, env)
        .await
}
