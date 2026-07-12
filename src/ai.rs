//! Natural language decomposition using Workers AI (milestone 2).
//!
//! This module turns a plain English request into the same typed scans the
//! engine already runs. The model never sees the DNA sequence, only the user's
//! text, so calls stay small and cheap. The model is asked to return strict JSON
//! (via Workers AI JSON mode), and when the request is too vague it returns a
//! clarifying question instead of guessing.
//!
//! The split mirrors the rest of the crate: building the request, parsing the
//! reply, and validating it are pure functions that run under `cargo test`. Only
//! [`decompose`] talks to the Worker runtime, and it is compiled for `wasm32`.

// On a plain host build the only caller of these items is the wasm-only Worker,
// so they look dead even though tests exercise them. Allow that off-wasm only;
// on the real wasm target normal dead-code checks still apply.
#![cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]

use serde::{Deserialize, Serialize};

use crate::error::EngineError;
use crate::pattern::compile_iupac;
use crate::types::ScanSpec;

/// The Workers AI model used for decomposition.
///
/// Small and fast, supports JSON mode, and is more than enough for mapping a
/// short prompt to our small schema. Swap this constant to change models.
pub(crate) const MODEL: &str = "@cf/meta/llama-3.1-8b-instruct-fp8-fast";

/// Instructions given to the model on every decomposition call.
const SYSTEM_PROMPT: &str = "\
You translate a user's request into DNA scan instructions for a sequence scanner.

The scanner supports exactly these scan types:
1. \"motif\": search for a nucleotide pattern written in IUPAC codes.
   IUPAC codes: A C G T are bases; R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C,
   B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G, N=any base. Fields: \"pattern\" (string),
   \"both_strands\" (boolean, default true).
2. \"orf\": find open reading frames (genes) from a start codon to a stop codon.
   Fields: \"min_aa\" (integer minimum protein length, default 30),
   \"both_strands\" (boolean, default true).
3. \"gc\": find GC-rich regions or CpG islands. Fields: \"min_len\" (integer minimum
   region length), \"min_gc\" (number 0..1, default 0.5), \"min_cpg_ratio\" (number,
   use 0 for a plain GC-rich region, 0.6 for a CpG island).
4. \"gc_skew\": locate the replication origin from GC skew. Field: \"window\" (integer
   window size, default 100).
5. \"restriction\": map restriction enzyme sites. Fields: \"enzymes\" (array of enzyme
   names, empty array means all known enzymes), \"both_strands\" (boolean, default true).
   Known enzymes include EcoRI, BamHI, HindIII, NotI, XhoI, PstI, SmaI, KpnI, SalI,
   XbaI, NcoI, NdeI, EcoRV, HaeIII, AluI.
6. \"pwm\": score a probabilistic motif from example sites. Fields: \"sites\" (array of
   equal-length example sequences), \"threshold\" (number 0..1, default 0.8),
   \"both_strands\" (boolean, default true).

Useful knowledge:
- \"TATA box\" or \"promoter TATA\" means the motif pattern TATAWAWR.
- \"genes\", \"ORFs\", \"reading frames\", or \"coding regions\" mean an orf scan.
- \"GC rich regions\" mean a gc scan with min_cpg_ratio 0. \"CpG islands\" mean a gc
  scan with min_cpg_ratio about 0.6.
- \"GC skew\", \"replication origin\", or \"ori\" mean a gc_skew scan.
- \"restriction sites\", \"cut sites\", or a named enzyme like EcoRI mean a restriction
  scan (put named enzymes in \"enzymes\").
- \"transcription factor binding site\", \"PWM\", or \"scored motif\" mean a pwm scan. Only
  use pwm if you can give several concrete equal-length example sites; otherwise ask.
- \"both strands\" or no mention means both_strands true; \"forward only\" means false.
- Always read any number in the request into the matching field and never leave it
  at the default when a number is given. \"greater than 10 length\" sets min_len to 11;
  \"at least 25 bp\" sets min_len to 25; \"genes over 50 aa\" sets min_aa to 50.

Rules:
- Only produce scans you are confident about. If the request is vague, or you
  cannot map it to a concrete scan type or fill its required fields, do not guess.
  Return an empty \"scans\" array and put one short, specific question in \"clarification\".
- If the request is clear, return the scans and leave \"clarification\" empty.
- Never invent scan types or fields beyond those listed above.

Respond only with the required JSON object.";

/// The result of interpreting a user prompt.
///
/// Derives serde so it can be cached in KV as JSON, letting a repeated prompt
/// skip the AI call entirely.
#[derive(Serialize, Deserialize)]
pub(crate) enum Decomposition {
    /// The request was too vague. Carries a question to show the user.
    NeedsClarification(String),
    /// The request mapped to one or more scans to run.
    Scans(Vec<ScanSpec>),
}

/// Build a stable cache key for a prompt.
///
/// Normalizes so that prompts that differ only in case or spacing share one
/// cache entry: trims ends, lowercases, and collapses runs of whitespace to a
/// single space.
///
/// # Arguments
///
/// * `prompt` - The user's natural language request.
///
/// # Returns
///
/// A normalized key string suitable for use as a KV key.
pub(crate) fn cache_key(prompt: &str) -> String {
    prompt.split_whitespace().collect::<Vec<_>>().join(" ").to_lowercase()
}

// --- Request types (sent to Workers AI) ------------------------------------

/// The chat request body, matching the Workers AI JSON mode shape.
#[derive(Serialize)]
struct ChatRequest {
    messages: Vec<Message>,
    response_format: ResponseFormat,
}

/// One chat message.
#[derive(Serialize)]
struct Message {
    role: &'static str,
    content: String,
}

/// The JSON mode wrapper: `type` is `json_schema` and `json_schema` is the schema.
#[derive(Serialize)]
struct ResponseFormat {
    #[serde(rename = "type")]
    kind: &'static str,
    json_schema: serde_json::Value,
}

// --- Reply types (returned by the model) -----------------------------------

/// The model's decomposition, before validation.
#[derive(Deserialize)]
struct AiDecomposition {
    /// A question to ask when the request is unclear. Empty when scans are given.
    #[serde(default)]
    clarification: String,
    /// The scans the model chose. Empty when it needs clarification.
    #[serde(default)]
    scans: Vec<AiScan>,
}

/// One scan from the model, in a flat shape that is easy for the model to fill.
/// Converted and validated into a [`ScanSpec`] by [`interpret`].
#[derive(Deserialize)]
struct AiScan {
    #[serde(rename = "type")]
    kind: String,
    #[serde(default)]
    pattern: Option<String>,
    #[serde(default)]
    min_aa: Option<usize>,
    #[serde(default = "default_true")]
    both_strands: bool,
    #[serde(default)]
    min_len: Option<usize>,
    #[serde(default)]
    min_gc: Option<f64>,
    #[serde(default)]
    min_cpg_ratio: Option<f64>,
    #[serde(default)]
    window: Option<usize>,
    #[serde(default)]
    enzymes: Vec<String>,
    #[serde(default)]
    sites: Vec<String>,
    #[serde(default)]
    threshold: Option<f64>,
}

fn default_true() -> bool {
    true
}

/// Build the Workers AI request for a prompt.
///
/// Combines the fixed system instructions with the user's text and attaches the
/// JSON schema that constrains the reply.
///
/// # Arguments
///
/// * `prompt` - The user's natural language request.
///
/// # Returns
///
/// A serializable [`ChatRequest`] ready to pass to `Ai::run`.
fn build_request(prompt: &str) -> ChatRequest {
    ChatRequest {
        messages: vec![
            Message {
                role: "system",
                content: SYSTEM_PROMPT.to_string(),
            },
            Message {
                role: "user",
                content: prompt.to_string(),
            },
        ],
        response_format: ResponseFormat {
            kind: "json_schema",
            json_schema: output_schema(),
        },
    }
}

/// Return the JSON schema the model must fill.
///
/// # Returns
///
/// A JSON schema value describing a `clarification` string and a `scans` array.
fn output_schema() -> serde_json::Value {
    serde_json::json!({
        "type": "object",
        "properties": {
            "clarification": { "type": "string" },
            "scans": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "type": {
                            "type": "string",
                            "enum": ["motif", "orf", "gc", "gc_skew", "restriction", "pwm"]
                        },
                        "pattern": { "type": "string" },
                        "min_aa": { "type": "integer" },
                        "min_len": { "type": "integer" },
                        "min_gc": { "type": "number" },
                        "min_cpg_ratio": { "type": "number" },
                        "window": { "type": "integer" },
                        "enzymes": { "type": "array", "items": { "type": "string" } },
                        "sites": { "type": "array", "items": { "type": "string" } },
                        "threshold": { "type": "number" },
                        "both_strands": { "type": "boolean" }
                    },
                    "required": ["type"]
                }
            }
        },
        "required": ["clarification", "scans"]
    })
}

/// Turn the raw model reply into an [`AiDecomposition`].
///
/// Workers AI may return the structured result either as a JSON object or as a
/// JSON string. This accepts both.
///
/// # Arguments
///
/// * `response` - The model output, either an object or a string holding JSON.
///
/// # Returns
///
/// The parsed [`AiDecomposition`].
///
/// # Errors
///
/// Returns [`EngineError::InvalidRequest`] if the reply is not valid JSON or does
/// not match the expected shape.
fn parse_ai_response(response: serde_json::Value) -> Result<AiDecomposition, EngineError> {
    let value = match response {
        serde_json::Value::String(s) => serde_json::from_str(&s).map_err(|e| {
            EngineError::InvalidRequest(format!("model did not return valid JSON: {e}"))
        })?,
        other => other,
    };
    serde_json::from_value(value)
        .map_err(|e| EngineError::InvalidRequest(format!("model JSON did not match schema: {e}")))
}

/// Validate a decomposition and convert it into typed scans.
///
/// Empty scans mean the model wants clarification. Otherwise each scan is checked
/// and turned into a [`ScanSpec`]. Motif patterns are validated with
/// [`compile_iupac`] so an invalid pattern is caught here, not at scan time.
///
/// # Arguments
///
/// * `ai` - The parsed model reply.
///
/// # Returns
///
/// [`Decomposition::NeedsClarification`] with a question, or
/// [`Decomposition::Scans`] with validated scans.
///
/// # Errors
///
/// Returns [`EngineError::InvalidRequest`] if a scan has an unknown type or a
/// motif scan is missing its pattern. Returns [`EngineError::EmptyPattern`] or
/// [`EngineError::InvalidIupacCode`] if a motif pattern is not valid IUPAC.
fn interpret(ai: AiDecomposition) -> Result<Decomposition, EngineError> {
    if ai.scans.is_empty() {
        let question = ai.clarification.trim();
        let message = if question.is_empty() {
            "Could you say which motif or feature you want to search for?".to_string()
        } else {
            question.to_string()
        };
        return Ok(Decomposition::NeedsClarification(message));
    }

    let mut specs = Vec::with_capacity(ai.scans.len());
    for scan in ai.scans {
        match scan.kind.as_str() {
            "motif" => {
                let pattern = scan
                    .pattern
                    .filter(|p| !p.trim().is_empty())
                    .ok_or_else(|| {
                        EngineError::InvalidRequest("motif scan is missing a pattern".to_string())
                    })?;
                compile_iupac(&pattern)?; // validate; masks are recomputed at scan time
                specs.push(ScanSpec::Motif {
                    pattern,
                    both_strands: scan.both_strands,
                });
            }
            "orf" => specs.push(ScanSpec::Orf {
                min_aa: scan.min_aa.unwrap_or(30),
                both_strands: scan.both_strands,
            }),
            "gc" => specs.push(ScanSpec::Gc {
                min_len: scan.min_len.unwrap_or(200),
                min_gc: scan.min_gc.unwrap_or(0.5),
                min_cpg_ratio: scan.min_cpg_ratio.unwrap_or(0.0),
            }),
            "gc_skew" => specs.push(ScanSpec::GcSkew {
                window: scan.window.unwrap_or(100),
            }),
            "restriction" => specs.push(ScanSpec::Restriction {
                enzymes: scan.enzymes,
                both_strands: scan.both_strands,
            }),
            "pwm" => {
                if scan.sites.is_empty() {
                    return Err(EngineError::InvalidRequest(
                        "pwm scan needs example sites".to_string(),
                    ));
                }
                specs.push(ScanSpec::Pwm {
                    sites: scan.sites,
                    threshold: scan.threshold.unwrap_or(0.8),
                    both_strands: scan.both_strands,
                });
            }
            other => {
                return Err(EngineError::InvalidRequest(format!(
                    "unknown scan type '{other}'"
                )))
            }
        }
    }
    Ok(Decomposition::Scans(specs))
}

/// Decompose a natural language prompt into scans using Workers AI.
///
/// Calls the model, parses the reply, and validates it. The DNA sequence is not
/// sent to the model.
///
/// # Arguments
///
/// * `env` - The Worker environment, used to reach the `AI` binding.
/// * `prompt` - The user's natural language request.
///
/// # Returns
///
/// A [`Decomposition`]: either a clarifying question or validated scans.
///
/// # Errors
///
/// Returns a [`worker::Error`] if the AI call fails, the reply cannot be parsed,
/// or a scan fails validation.
#[cfg(target_arch = "wasm32")]
pub(crate) async fn decompose(
    env: &worker::Env,
    prompt: &str,
) -> std::result::Result<Decomposition, worker::Error> {
    let ai = env.ai("AI")?;
    let request = build_request(prompt);
    // The binding may return the model output directly or wrapped in `response`.
    let raw: serde_json::Value = ai.run(MODEL, request).await?;
    let payload = raw.get("response").cloned().unwrap_or(raw);
    let decomposition =
        parse_ai_response(payload).map_err(|e| worker::Error::RustError(e.to_string()))?;
    interpret(decomposition).map_err(|e| worker::Error::RustError(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(json: &str) -> AiDecomposition {
        parse_ai_response(serde_json::from_str(json).unwrap()).unwrap()
    }

    #[test]
    fn cache_key_normalizes_case_and_whitespace() {
        assert_eq!(cache_key("  Find   TATA  Boxes "), "find tata boxes");
        assert_eq!(cache_key("find TATA boxes"), cache_key("FIND tata   BOXES"));
    }

    #[test]
    fn decomposition_round_trips_through_json() {
        let d = Decomposition::Scans(vec![ScanSpec::Orf {
            min_aa: 30,
            both_strands: true,
        }]);
        let text = serde_json::to_string(&d).unwrap();
        let back: Decomposition = serde_json::from_str(&text).unwrap();
        assert!(matches!(back, Decomposition::Scans(_)));
    }

    #[test]
    fn builds_request_with_schema_and_prompt() {
        let json = serde_json::to_string(&build_request("find TATA boxes")).unwrap();
        assert!(json.contains("json_schema"));
        assert!(json.contains("find TATA boxes"));
        assert!(json.contains("\"scans\""));
        assert!(MODEL.starts_with("@cf/"));
    }

    #[test]
    fn parses_object_reply() {
        let d = parse(r#"{"clarification":"","scans":[{"type":"orf"}]}"#);
        assert_eq!(d.scans.len(), 1);
    }

    #[test]
    fn parses_string_reply() {
        // Workers AI sometimes returns the JSON as a string.
        let wrapped = serde_json::Value::String(
            r#"{"clarification":"","scans":[{"type":"orf"}]}"#.to_string(),
        );
        assert_eq!(parse_ai_response(wrapped).unwrap().scans.len(), 1);
    }

    #[test]
    fn interprets_gc_scan() {
        let d = parse(r#"{"clarification":"","scans":[
            {"type":"gc","min_len":11,"min_cpg_ratio":0.6}]}"#);
        match interpret(d).unwrap() {
            Decomposition::Scans(s) => {
                assert!(matches!(s[0], ScanSpec::Gc { min_len: 11, .. }));
            }
            _ => panic!("expected scans"),
        }
    }

    #[test]
    fn interprets_restriction_and_skew() {
        let d = parse(r#"{"clarification":"","scans":[
            {"type":"restriction","enzymes":["EcoRI"]},
            {"type":"gc_skew","window":50}]}"#);
        match interpret(d).unwrap() {
            Decomposition::Scans(s) => assert_eq!(s.len(), 2),
            _ => panic!("expected scans"),
        }
    }

    #[test]
    fn pwm_without_sites_is_rejected() {
        let d = parse(r#"{"clarification":"","scans":[{"type":"pwm"}]}"#);
        assert!(matches!(interpret(d), Err(EngineError::InvalidRequest(_))));
    }

    #[test]
    fn interprets_motif_scan() {
        let d = parse(r#"{"clarification":"","scans":[
            {"type":"motif","pattern":"TATAWAWR","both_strands":true}]}"#);
        match interpret(d).unwrap() {
            Decomposition::Scans(specs) => assert_eq!(specs.len(), 1),
            _ => panic!("expected scans"),
        }
    }

    #[test]
    fn defaults_orf_min_aa() {
        let d = parse(r#"{"clarification":"","scans":[{"type":"orf"}]}"#);
        match interpret(d).unwrap() {
            Decomposition::Scans(specs) => match &specs[0] {
                ScanSpec::Orf { min_aa, both_strands } => {
                    assert_eq!(*min_aa, 30);
                    assert!(*both_strands);
                }
                _ => panic!("expected orf"),
            },
            _ => panic!("expected scans"),
        }
    }

    #[test]
    fn empty_scans_becomes_clarification() {
        let d = parse(r#"{"clarification":"Which motif did you mean?","scans":[]}"#);
        match interpret(d).unwrap() {
            Decomposition::NeedsClarification(q) => assert_eq!(q, "Which motif did you mean?"),
            _ => panic!("expected clarification"),
        }
    }

    #[test]
    fn empty_scans_without_question_gets_default_question() {
        let d = parse(r#"{"clarification":"","scans":[]}"#);
        match interpret(d).unwrap() {
            Decomposition::NeedsClarification(q) => assert!(!q.is_empty()),
            _ => panic!("expected clarification"),
        }
    }

    #[test]
    fn invalid_pattern_is_rejected() {
        let d = parse(r#"{"clarification":"","scans":[{"type":"motif","pattern":"TAXQ"}]}"#);
        assert!(matches!(
            interpret(d),
            Err(EngineError::InvalidIupacCode('X'))
        ));
    }

    #[test]
    fn motif_without_pattern_is_rejected() {
        let d = parse(r#"{"clarification":"","scans":[{"type":"motif"}]}"#);
        assert!(matches!(interpret(d), Err(EngineError::InvalidRequest(_))));
    }

    #[test]
    fn unknown_scan_type_is_rejected() {
        let d = parse(r#"{"clarification":"","scans":[{"type":"gizmo"}]}"#);
        assert!(matches!(interpret(d), Err(EngineError::InvalidRequest(_))));
    }
}
