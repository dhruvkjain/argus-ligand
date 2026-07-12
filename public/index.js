// ---------------------------------------------------------------------------
// Demo data
// ---------------------------------------------------------------------------
const EXAMPLES = {
  tata: {
    seq: ">promoter_demo\nGGGCGCGCCTATATAAGGCTCGAGCCTATAAAACGTAGGGATATATATCGCG",
    scans: [{ type: "motif", pattern: "TATAWAWR", both_strands: true }],
  },
  orf: {
    seq: ">orf_demo\nCCGATGGCTAAAGATCTGCGTAAACATTAGCCCATGAAAGGGTTTTGATAA",
    scans: [{ type: "orf", min_aa: 3, both_strands: true }],
  },
};

const scansEl = document.getElementById("scans");

// ---------------------------------------------------------------------------
// Manual scan rows
// ---------------------------------------------------------------------------
function addScanRow(spec = { type: "motif", pattern: "GAATTC", both_strands: true }) {
  const div = document.createElement("div");
  div.className = "scan";
  div.innerHTML = `
    <select class="kind">
      <option value="motif">motif</option>
      <option value="orf">orf</option>
      <option value="gc">gc / cpg</option>
      <option value="gc_skew">gc skew</option>
      <option value="restriction">restriction</option>
      <option value="pwm">pwm</option>
    </select>
    <input class="param" />
    <label class="both"><input type="checkbox" checked /> both strands</label>
    <button class="rm" title="remove" type="button">×</button>`;
  const kind = div.querySelector(".kind");
  const param = div.querySelector(".param");
  const both = div.querySelector(".both input");
  const PLACEHOLDER = {
    motif: "IUPAC pattern e.g. TATAWAWR",
    orf: "min protein length (aa) e.g. 30",
    gc: "min region length e.g. 200",
    gc_skew: "window size e.g. 100",
    restriction: "enzyme names, blank = all e.g. EcoRI, BamHI",
    pwm: "example sites e.g. TATAAA TATATA TATAAT",
  };
  const DEFAULT = { motif: "TATAWAWR", orf: "30", gc: "200", gc_skew: "100", restriction: "", pwm: "TATAAA TATATA TATAAT" };
  const syncPlaceholder = () => { param.placeholder = PLACEHOLDER[kind.value]; };
  kind.value = spec.type;
  both.checked = spec.both_strands !== false;
  if (spec.type === "motif") param.value = spec.pattern || "";
  else if (spec.type === "orf") param.value = spec.min_aa ?? 30;
  else param.value = DEFAULT[spec.type] ?? "";
  syncPlaceholder();
  kind.onchange = () => { param.value = DEFAULT[kind.value] ?? ""; syncPlaceholder(); };
  div.querySelector(".rm").onclick = () => div.remove();
  scansEl.appendChild(div);
}

function readScans() {
  return [...scansEl.querySelectorAll(".scan")].map((div) => {
    const type = div.querySelector(".kind").value;
    const both_strands = div.querySelector(".both input").checked;
    const param = div.querySelector(".param").value.trim();
    const list = () => param.split(/[\s,]+/).filter(Boolean);
    switch (type) {
      case "motif": return { type, pattern: param, both_strands };
      case "orf": return { type, min_aa: parseInt(param || "30", 10), both_strands };
      case "gc": return { type, min_len: parseInt(param || "200", 10) };
      case "gc_skew": return { type, window: parseInt(param || "100", 10) };
      case "restriction": return { type, enzymes: list(), both_strands };
      case "pwm": return { type, sites: list(), both_strands };
      default: return { type };
    }
  });
}

// ---------------------------------------------------------------------------
// Result rendering
// ---------------------------------------------------------------------------
const strandClass = (s) => (s === "-" ? "strand-minus" : s === "±" ? "strand-pal" : "strand-plus");
const fmt = (n, d) => Number(n).toFixed(d);

function renderMotif(r) {
  const rows = r.matches
    .map((m) => `<tr><td>${m.start}-${m.end}</td><td class="${strandClass(m.strand)}">${m.strand}</td><td>${m.matched}</td></tr>`)
    .join("");
  const hasPal = r.matches.some((m) => m.strand === "±");
  const note = hasPal
    ? `<div class="hint"><span class="strand-pal">±</span> = palindromic site: reads the same on both strands, so it is counted once instead of twice.</div>`
    : "";
  const body = r.matches.length
    ? `<table><thead><tr><th>span</th><th>strand</th><th>match</th></tr></thead><tbody>${rows}</tbody></table>${note}`
    : `<div class="empty">no matches</div>`;
  return `<div class="result-card"><h3>motif <code>${r.pattern}</code><span class="badge">${r.count} hit(s)</span></h3>${body}</div>`;
}

function renderOrf(r) {
  const rows = r.matches
    .map((m) => `<tr><td>${m.start}-${m.end}</td><td class="${strandClass(m.strand)}">${m.strand}${m.frame}</td><td>${m.aa_length}</td><td class="prot">${m.protein}</td></tr>`)
    .join("");
  const body = r.matches.length
    ? `<table><thead><tr><th>span</th><th>strand/frame</th><th>aa</th><th>protein</th></tr></thead><tbody>${rows}</tbody></table>`
    : `<div class="empty">no ORFs found</div>`;
  return `<div class="result-card"><h3>ORFs <span class="badge">${r.count} found</span></h3>${body}</div>`;
}

function renderGc(r) {
  const rows = r.matches
    .map((m) => `<tr><td>${m.start}-${m.end}</td><td>${m.length}</td><td>${fmt(m.gc_percent, 1)}%</td><td>${fmt(m.cpg_ratio, 2)}</td></tr>`)
    .join("");
  const body = r.matches.length
    ? `<table><thead><tr><th>span</th><th>length</th><th>GC</th><th>CpG o/e</th></tr></thead><tbody>${rows}</tbody></table>`
    : `<div class="empty">no GC-rich regions found</div>`;
  return `<div class="result-card"><h3>GC-rich / CpG islands <span class="badge">${r.count} region(s)</span></h3>${body}</div>`;
}

function renderGcSkew(r) {
  return `<div class="result-card"><h3>GC skew <span class="badge">window ${r.window}</span></h3>
    <table><tbody>
      <tr><th>putative origin</th><td class="strand-minus">position ${r.origin}</td></tr>
      <tr><th>putative terminus</th><td class="strand-plus">position ${r.terminus}</td></tr>
      <tr><th>cumulative skew</th><td>${fmt(r.min_skew, 2)} to ${fmt(r.max_skew, 2)}</td></tr>
    </tbody></table>
    <div class="hint">Origin is where cumulative GC skew is lowest; terminus where it is highest.</div></div>`;
}

function renderRestriction(r) {
  const rows = r.matches
    .map((m) => `<tr><td>${m.enzyme}</td><td>${m.site}</td><td>${m.start}-${m.end}</td><td class="${strandClass(m.strand)}">${m.strand}</td></tr>`)
    .join("");
  const body = r.matches.length
    ? `<table><thead><tr><th>enzyme</th><th>site</th><th>span</th><th>strand</th></tr></thead><tbody>${rows}</tbody></table>`
    : `<div class="empty">no restriction sites found</div>`;
  return `<div class="result-card"><h3>Restriction sites <span class="badge">${r.count} found</span></h3>${body}</div>`;
}

function renderPwm(r) {
  const rows = r.matches
    .map((m) => `<tr><td>${m.start}-${m.end}</td><td class="${strandClass(m.strand)}">${m.strand}</td><td>${fmt(m.score, 2)}</td></tr>`)
    .join("");
  const body = r.matches.length
    ? `<table><thead><tr><th>span</th><th>strand</th><th>score</th></tr></thead><tbody>${rows}</tbody></table>`
    : `<div class="empty">no matches above threshold</div>`;
  return `<div class="result-card"><h3>PWM hits <span class="badge">${r.count} · max score ${fmt(r.max_score, 1)}</span></h3>${body}</div>`;
}

function renderOne(r) {
  switch (r.type) {
    case "motif": return renderMotif(r);
    case "orf": return renderOrf(r);
    case "gc": return renderGc(r);
    case "gc_skew": return renderGcSkew(r);
    case "restriction": return renderRestriction(r);
    case "pwm": return renderPwm(r);
    default: return `<div class="err">${r.message || "unknown result"}</div>`;
  }
}

function renderResults(data) {
  let html = `<div class="badge">sequence length: ${data.seq_length} nt</div>`;
  if (data.warnings?.length) html += `<div class="warn">⚠ ${data.warnings.join("; ")}</div>`;
  html += data.results.map(renderOne).join("");
  return html;
}

function renderInterpreted(scans) {
  const chips = scans
    .map((s) => {
      let label;
      switch (s.type) {
        case "motif": label = `motif ${s.pattern}`; break;
        case "orf": label = `ORF ≥ ${s.min_aa ?? 30} aa`; break;
        case "gc": label = `${s.min_cpg_ratio >= 0.5 ? "CpG islands" : "GC-rich"} ≥ ${s.min_len ?? 200} nt`; break;
        case "gc_skew": label = `GC skew · window ${s.window ?? 100}`; break;
        case "restriction": label = s.enzymes?.length ? `restriction: ${s.enzymes.join(", ")}` : "all restriction sites"; break;
        case "pwm": label = `PWM · ${(s.sites || []).length} sites`; break;
        default: label = s.type;
      }
      const strandable = ["motif", "orf", "restriction", "pwm"].includes(s.type);
      const strands = strandable ? (s.both_strands === false ? " · one strand" : " · both strands") : "";
      return `<span class="chip">${label}${strands}</span>`;
    })
    .join("");
  return `<div class="interp"><span class="badge">Interpreted your request as</span><div class="chips">${chips}</div></div>`;
}

// ---------------------------------------------------------------------------
// Actions
// ---------------------------------------------------------------------------
async function run() {
  const runBtn = document.getElementById("run");
  const out = document.getElementById("out");
  const sequence = document.getElementById("seq").value;
  const scans = readScans();
  if (!sequence.trim()) { out.innerHTML = `<div class="err">Paste a sequence first.</div>`; return; }
  if (!scans.length) { out.innerHTML = `<div class="err">Add at least one scan.</div>`; return; }
  runBtn.disabled = true;
  runBtn.textContent = "Scanning…";
  try {
    const res = await fetch("/scan", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ sequence, scans }),
    });
    const data = await res.json();
    if (data.error) { out.innerHTML = `<div class="err">Error: ${data.error}</div>`; return; }
    out.innerHTML = renderResults(data);
  } catch (e) {
    out.innerHTML = `<div class="err">Request failed: ${e.message}</div>`;
  } finally {
    runBtn.disabled = false;
    runBtn.textContent = "Scan";
  }
}

function updateQuota(remaining, max) {
  const q = document.getElementById("quota");
  if (!q || typeof remaining !== "number") return;
  q.textContent = `AI · ${remaining}/${max ?? 5} today`;
  const color = remaining === 0 ? "var(--danger)" : remaining <= 1 ? "var(--amber)" : "var(--muted)";
  q.style.color = color;
  q.style.borderColor = color;
}

async function ask() {
  const askBtn = document.getElementById("ask");
  const out = document.getElementById("out");
  const prompt = document.getElementById("prompt").value.trim();
  const sequence = document.getElementById("seq").value;
  if (!prompt) { out.innerHTML = `<div class="err">Type a request first.</div>`; return; }
  if (!sequence.trim()) { out.innerHTML = `<div class="err">Paste a sequence first.</div>`; return; }
  askBtn.disabled = true;
  askBtn.textContent = "Asking…";
  try {
    const res = await fetch("/ask", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ prompt, sequence }),
    });
    const data = await res.json();
    if (data.error) { out.innerHTML = `<div class="err">Error: ${data.error}</div>`; return; }
    if (data.limit_exceeded) {
      updateQuota(0, data.max);
      out.innerHTML = `<div class="clarify"><strong>Limit reached:</strong> ${data.message}</div>`;
      return;
    }
    updateQuota(data.remaining, data.max);
    const cachedBadge = data.cached ? `<div class="badge">⚡ answered from cache · 0 neurons</div>` : "";
    if (data.clarification) {
      out.innerHTML = cachedBadge + `<div class="clarify"><strong>I need a bit more detail:</strong> ${data.clarification}</div>`;
      return;
    }
    out.innerHTML = cachedBadge + renderInterpreted(data.interpreted) + renderResults(data.result);
  } catch (e) {
    out.innerHTML = `<div class="err">Request failed: ${e.message}</div>`;
  } finally {
    askBtn.disabled = false;
    askBtn.textContent = "Ask";
  }
}

// ---------------------------------------------------------------------------
// Wiring
// ---------------------------------------------------------------------------
document.getElementById("add").onclick = () => addScanRow();
document.getElementById("run").onclick = run;
document.getElementById("ask").onclick = ask;

document.querySelectorAll(".examples-ex [data-ex]").forEach((b) => {
  b.onclick = () => {
    const ex = b.dataset.ex;
    if (ex === "clear") { document.getElementById("seq").value = ""; return; }
    const e = EXAMPLES[ex];
    document.getElementById("seq").value = e.seq;
    scansEl.innerHTML = "";
    e.scans.forEach(addScanRow);
  };
});

document.querySelectorAll("[data-ask]").forEach((b) => {
  b.onclick = () => { document.getElementById("prompt").value = b.dataset.ask; };
});

const askPanel = document.getElementById("askPanel");
const manualPanel = document.getElementById("manualPanel");
const tabAsk = document.getElementById("tab-ask");
const tabManual = document.getElementById("tab-manual");
function setMode(mode) {
  const ai = mode === "ask";
  askPanel.hidden = !ai;
  manualPanel.hidden = ai;
  tabAsk.classList.toggle("active", ai);
  tabManual.classList.toggle("active", !ai);
}
tabAsk.onclick = () => setMode("ask");
tabManual.onclick = () => setMode("manual");

const sampleSel = document.getElementById("sample");
fetch("/samples")
  .then((r) => r.json())
  .then((list) => {
    list.forEach((s) => {
      const opt = document.createElement("option");
      opt.value = s.id;
      opt.textContent = s.label;
      sampleSel.appendChild(opt);
    });
  })
  .catch(() => {});
sampleSel.onchange = async () => {
  if (!sampleSel.value) return;
  sampleSel.disabled = true;
  try {
    const res = await fetch(`/samples/${sampleSel.value}`);
    document.getElementById("seq").value = await res.text();
  } catch (e) {
    document.getElementById("out").innerHTML = `<div class="err">Could not load sample: ${e.message}</div>`;
  } finally {
    sampleSel.disabled = false;
  }
};

// Show the real remaining AI budget on load. The limit is server-side (per IP),
// so a refresh does not reset it; this just makes the badge tell the truth.
fetch("/quota")
  .then((r) => r.json())
  .then((d) => updateQuota(d.remaining, d.max))
  .catch(() => {});

// Seed: TATA-box sequence, a matching manual scan, and a ready prompt.
document.getElementById("seq").value = EXAMPLES.tata.seq;
document.getElementById("prompt").value = "show me all CpG islands and any genes";
EXAMPLES.tata.scans.forEach(addScanRow);
