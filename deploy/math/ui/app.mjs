import {
  PAGE_SIZE,
  cycleText,
  passportText,
  status,
  filterEntries,
  readFilters,
  filtersURL,
  smoothness,
  source,
  SOURCE_LABELS,
  cuspGraphs,
  cuspDrawing,
} from "./catalog.mjs?v=20261109-v20";
const $ = (id) => document.getElementById(id);
// One page, two catalogs: /all-smooth (all-smooth only) and /all-one-piece (everything).
const MODE = document.documentElement.dataset.catalog === "one_piece" ? "one_piece" : "all_smooth";
const BASE = MODE === "one_piece" ? "/all-one-piece" : "/all-smooth";
const DEFAULT_ISO = MODE === "one_piece" ? "all" : "non-isotrivial";
const form = $("filters"),
  number = (n) => n == null ? "Not computed" : n.toLocaleString("en-US");
const escape = (text) =>
  String(text).replace(
    /[&<>"']/g,
    (c) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[
        c
      ],
  );
const cycleHTML = (a) =>
  escape(cycleText(a)).replace(/\^(\d+)/g, "<sup>$1</sup>");
const passport = (e) => "(" + (e.passport || (e.cycle_types || []).map(cycleText))
  .map((value) => escape(value).replace(/\^(\d+)/g, "<sup>$1</sup>"))
  .join(", ") + ")";
const orbit = (a) => a.map(number).join(" + ");
const kind = (e) =>
  smoothness(e) === "one_piece"
    ? `<span class="badge one-piece" title="One-piece at every state, not all-smooth">Not all-smooth</span>`
    : `<span class="badge all-smooth" title="Index additivity at every state">All-smooth</span>`;
const badge = (e) =>
  `<span class="badge ${status(e)}" title="${escape(status(e))}">${e.isotrivial === false ? "No" : e.isotrivial === true ? "Yes" : "Unknown"}</span>`;
let dataset,
  selected = [],
  page = 0,
  timer,
  activeId = null;
const detailShardCache = new Map(), detailById = new Map();
let modeEntries = null;
const entries = () =>
  (modeEntries ||= MODE === "one_piece"
    ? dataset.components
    : dataset.components.filter((e) => smoothness(e) === "all_smooth"));
function values() {
  const f = Object.fromEntries(new FormData(form));
  // The all-smooth page has no one-piece filters; keep its links clean.
  if (MODE === "all_smooth") {
    delete f.smoothness;
    delete f.source;
  }
  return f;
}
function setFilters(f) {
  HTMLFormElement.prototype.reset.call(form);
  for (const [key, v] of Object.entries(f)) {
    const field = form.elements.namedItem(key);
    if (field) field.value = v;
  }
}
function info(text) {
  $("notice").textContent = text;
  document
    .querySelectorAll(".detail-copy-note")
    .forEach((e) => (e.textContent = text));
}
async function copy(text) {
  try {
    await navigator.clipboard.writeText(text);
    info("Copied to clipboard.");
  } catch {
    info(
      "Clipboard unavailable. Select and copy the page URL or download the JSON instead.",
    );
  }
}
function blobDownload(text, name) {
  const url = URL.createObjectURL(
    new Blob(Array.isArray(text) ? text : [text], { type: "application/json" }),
  );
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function countStats() {
  const counts = { "non-isotrivial": 0, isotrivial: 0, unknown: 0 };
  let notSmooth = 0, primitiveNotSmooth = 0, states = 0;
  entries().forEach((e) => {
    counts[status(e)]++;
    states += e.states;
    if (smoothness(e) === "one_piece") {
      notSmooth++;
      if (e.primitive) primitiveNotSmooth++;
    }
  });
  const cards = [
    [MODE === "one_piece" ? "One-piece components" : "All-smooth components", entries().length, `${number(states)} quotient states`],
    ["Non-isotrivial", counts["non-isotrivial"], MODE === "one_piece" ? "Proposition 1" : "Shown by default"],
    ["Isotrivial", counts.isotrivial, MODE === "one_piece" ? "Proposition 1 or genus 1" : "Proposition 1"],
    ["Unknown", counts.unknown, "Criterion does not decide"],
  ];
  if (MODE === "one_piece")
    cards.push(["Not all-smooth", notSmooth,
      `<a href="${BASE}?smoothness=one_piece">Show them</a> · <a href="${BASE}?smoothness=one_piece&primitive=true">${number(primitiveNotSmooth)} primitive</a>`]);
  $("stats").innerHTML = cards
    .map(
      ([label, n, note]) =>
        `<div class="stat"><strong>${number(n)}</strong><span>${label}</span><small>${note}</small></div>`,
    )
    .join("");
  $("version").textContent =
    `Dataset v${dataset.dataset_version} · ${new Date(dataset.last_updated).toLocaleString("en-US", { dateStyle: "medium", timeStyle: "short", timeZone: "America/Chicago" })} (Chicago)`;
}
function render(updateURL = true) {
  if (!dataset) return;
  clearTimeout(timer);
  try {
    const f = values();
    selected = filterEntries(entries(), f);
    $("error").hidden = true;
    page = Math.min(
      page,
      Math.max(0, Math.ceil(selected.length / PAGE_SIZE) - 1),
    );
    $("result-count").textContent =
      `${number(selected.length)} ${selected.length === 1 ? "component" : "components"}`;
    $("result-note").textContent =
      `${number(selected.reduce((s, e) => s + e.states, 0))} quotient states · ${f.isotrivial === "all" ? "all classifications" : f.isotrivial}`;
    $("rows").innerHTML = selected
      .slice(page * PAGE_SIZE, (page + 1) * PAGE_SIZE)
      .map(
        (e) =>
          `<tr class="data-row" data-component="${e.id}" tabindex="0" aria-label="Degree ${e.degree}, genus ${e.genus}, ${e.states} states. Open example."><td data-label="Degree">${e.degree}</td><td data-label="Genus">${e.genus}</td><td class="passport" data-label="Passport">${passport(e)}</td><td data-label="States">${number(e.states)}</td><td data-label="Group order">${number(e.group_order)}</td><td data-label="Isotrivial?">${badge(e)}</td><td class="kind" data-label="Type">${kind(e)}</td><td class="row-action"><button class="open-example" data-id="${e.id}" aria-expanded="false">View tuple ↗</button></td></tr>`,
      )
      .join("");
    $("empty").hidden = selected.length !== 0;
    $("download").disabled = selected.length === 0;
    $("page-info").textContent = selected.length
      ? `${number(page * PAGE_SIZE + 1)}–${number(Math.min((page + 1) * PAGE_SIZE, selected.length))} of ${number(selected.length)}`
      : "0 results";
    $("previous").disabled = page === 0;
    $("next").disabled = (page + 1) * PAGE_SIZE >= selected.length;
    $("top-previous").disabled = $("previous").disabled;
    $("top-next").disabled = $("next").disabled;
    $("top-page-info").textContent =
      `Page ${selected.length ? page + 1 : 0} of ${Math.ceil(selected.length / PAGE_SIZE)} · ${PAGE_SIZE} per page`;
    $("mobile-sort").value = f.sort;
    $("mobile-direction").textContent = f.direction === "desc" ? "↓" : "↑";
    document.querySelectorAll("[data-sort]").forEach((b) => {
      b.innerHTML =
        escape(
          {
            degree: "Degree",
            genus: "Genus",
            states: "States",
            group_order: "Group order",
          }[b.dataset.sort],
        ) +
        (b.dataset.sort === f.sort
          ? ` <span>${f.direction === "desc" ? "↓" : "↑"}</span>`
          : "");
      b.closest("th").setAttribute(
        "aria-sort",
        b.dataset.sort === f.sort
          ? f.direction === "desc"
            ? "descending"
            : "ascending"
          : "none",
      );
    });
    if (updateURL) {
      const p = filtersURL({
        ...f,
        ...(activeId ? { example: activeId } : {}),
      });
      history.replaceState(null, "", `${location.pathname}?${p}`);
    }
    if (activeId) {
      const button = [...document.querySelectorAll("[data-id]")].find(
        (b) => b.dataset.id === activeId,
      );
      if (button) openDetail(button, false);
    }
  } catch (e) {
    $("error").textContent = e.message;
    $("error").hidden = false;
    $("rows").innerHTML = "";
    selected = [];
    $("download").disabled = true;
    $("next").disabled = true;
    $("previous").disabled = true;
    $("top-next").disabled = true;
    $("top-previous").disabled = true;
    $("empty").hidden = true;
    $("result-count").textContent = "Invalid query";
    $("result-note").textContent = "Correct the query to continue.";
    $("page-info").textContent = "";
    $("top-page-info").textContent = "";
  }
}
function quotientDetail(e, b, i) {
  const ws = (e.isotriviality.witnesses || []).filter(
    (w) => w.block_quotient_index === i,
  );
  return `<details class="quotient"><summary>${b.block_count} blocks × ${b.block_size} points · quotient genus h = ${b.genus}${ws.length ? " · classification witness" : ""}</summary>
    <p class="small">Quotient group order ${number(b.group_order)} · kernel order ${number(b.kernel_order)}${b.kernel_structure ? " (" + escape(b.kernel_structure) + ")" : ""} · quotient action ${b.primitive == null ? "not computed" : b.primitive ? "primitive" : "imprimitive"}</p>
    <p class="math">2g − 2 = m(2h − 2) + Σeᵢ: ${2 * e.genus - 2} = ${b.block_size}(${2 * b.genus - 2}) + ${b.ramification_excess.reduce((a, v) => a + v, 0)}</p>
    <p>Quotient indices: (${b.indices.join(", ")}). Ramification excess: (${b.ramification_excess.join(", ")}).</p>
    <div class="tuple math">${b.tuple.map((t, k) => `<span>ḡ<sub>${k + 1}</sub></span><p>= ${escape(t)}</p>`).join("")}</div>
    <p class="math">Blocks: ${b.blocks.map((a) => "{" + a.join(", ") + "}").join("; ")}</p>
    ${ws.map((w) => `<p class="small">Proposition 1(a)–(c): g<sub>${w.special_entry}</sub> fixes every block setwise and moves points only in block ${w.active_block}. Ordinary cycle-count equalities hold. ${w.claim === 6 ? "Claim 6 (h = 1)" : "Claim 7 (h ≥ 2)"} applies.</p>`).join("")}
  </details>`;
}
function criterionFailures(e) {
  const systems = e.isotriviality.failed_block_systems;
  if (e.isotrivial !== null || !systems) return "";
  if (!systems.length)
    return "<p>No proper invariant block system is recorded; the criterion requires one.</p>";
  return `<div class="criterion-failures"><h3>Why (a)–(c) do not all hold</h3>${systems.map((b) => {
    const title = `${b.block_count} blocks × ${b.block_size} points`;
    if (!b.a_passing_entries.length)
      return `<p><strong>${escape(title)}:</strong> (a) fails for every entry. The four quotient cycle counts are (${b.quotient_cycle_counts.join(", ")}); fixing every block requires ${b.block_count} cycles.</p>`;
    return `<p><strong>${escape(title)}:</strong></p><ul>${b.candidate_failures.map((c) => {
      const bText = c.failed_conditions.includes("b")
        ? `(b) fails: moves points in ${c.active_blocks.length} blocks (${c.active_blocks.join(", ") || "none"}); exactly one is required.`
        : `(b) passes: moves points only in block ${c.active_blocks[0]}.`;
      const cText = c.ordinary_cycle_count_failures.length
        ? "(c) fails: " + c.ordinary_cycle_count_failures.map((f) =>
          `g${f.entry} has ${f.actual_cycles} cycles, but ${b.block_size} × ${f.quotient_cycles} = ${f.required_cycles} are required (ramification excess ${f.ramification_excess})`).join("; ") + "."
        : "(c) passes for all ordinary entries.";
      return `<li>Special entry g<sub>${c.special_entry}</sub>: (a) passes. ${escape(bText)} ${escape(cText)}</li>`;
    }).join("")}</ul>`;
  }).join("")}</div>`;
}
function parsePerm(text, n) {
  const p = Array.from({ length: n }, (_, i) => i);
  for (const match of text.matchAll(/\(([^)]*)\)/g)) {
    const cycle = match[1].trim().split(/\s+/).map((x) => Number(x) - 1);
    if (cycle.length === 1 && cycle[0] === -1) continue;
    for (let i = 0; i < cycle.length; i++) p[cycle[i]] = cycle[(i + 1) % cycle.length];
  }
  return p;
}
function cycleData(p) {
  const seen = new Uint8Array(p.length), cycles = [];
  for (let root = 0; root < p.length; root++) {
    if (seen[root]) continue;
    const cycle = []; let x = root;
    while (!seen[x]) { seen[x] = 1; cycle.push(x + 1); x = p[x]; }
    cycles.push(cycle);
  }
  return cycles;
}
function product(a, b) { return b.map((x) => a[x]); }
function permutationText(p) {
  return cycleData(p).filter((c) => c.length > 1).map((c) => `(${c.join(" ")})`).join("") || "()";
}
function orbitSizes(generators, n) {
  const seen = new Uint8Array(n), sizes = [];
  for (let root = 0; root < n; root++) {
    if (seen[root]) continue;
    seen[root] = 1; const stack = [root]; let size = 0;
    while (stack.length) {
      const x = stack.pop(); size++;
      for (const g of generators) { const y = g[x]; if (!seen[y]) { seen[y] = 1; stack.push(y); } }
    }
    sizes.push(size);
  }
  return sizes.sort((a, b) => b - a);
}
function expandDetail(base, detail) {
  const e = { ...base, ...detail };
  if (detail.block_witness) {
    const b = detail.block_witness;
    const ramification_excess = b.ramification_excess || e.indices.map((i, j) => i - b.block_size * b.indices[j]);
    e.block_quotients = [{ ...b, ramification_excess, group_order: null,
      kernel_order: null, kernel_structure: null, primitive: null, orders_status: "not_computed" }];
    e.block_quotients_complete = false;
  }
  if (!e.representative_g1_g2 || !e.representative_pair_orbit_sizes) {
    const perms = e.tuple.map((s) => parsePerm(s, e.degree));
    const p = product(perms[0], perms[1]), cycles = cycleData(p);
    e.representative_g1_g2 ||= { permutation: permutationText(p), cycle_count: cycles.length, index: e.degree - cycles.length };
    e.representative_pair_orbit_sizes ||= {
      g1_g2: orbitSizes(perms.slice(0, 2), e.degree),
      g2_g3: orbitSizes(perms.slice(1, 3), e.degree),
      g3_g4: orbitSizes(perms.slice(2, 4), e.degree),
    };
  }
  e.defect_distribution ||= [{ defect: 0, states: e.states }];
  e.all_states_smooth = smoothness(e) === "all_smooth";
  e.detail_url ||= `/details/${e.id[0]}.jsonl.gz`;
  return e;
}
function gcdBig(a, b) {
  while (b) [a, b] = [b, a % b];
  return a < 0n ? -a : a;
}
function addFraction(a, b) {
  const n = a.n * b.d + b.n * a.d, d = a.d * b.d, g = gcdBig(n, d);
  return { n: n / g, d: d / g };
}
function rhoOfCycles(lengths) {
  return lengths.reduce((sum, n) => addFraction(sum, { n: 1n, d: BigInt(n) }), { n: 0n, d: 1n });
}
function fractionText(x) { return x.d === 1n ? String(x.n) : String(x.n) + "/" + String(x.d); }
function meanRho(distribution, states) {
  let sum = { n: 0n, d: 1n };
  for (const row of distribution) {
    const rho = rhoOfCycles(row.cycle_lengths);
    sum = addFraction(sum, { n: rho.n * BigInt(row.states), d: rho.d });
  }
  const denominator = sum.d * BigInt(states), common = gcdBig(sum.n, denominator);
  return fractionText({ n: sum.n / common, d: denominator / common });
}
function leanRecord(e) {
  const fields = [
    "id", "revision", "last_updated", "degree", "genus", "states", "tuple",
    "passport", "indices", "primitive", "contained_in", "group_order",
    "isotrivial", "representative_g1_g2", "representative_pair_orbit_sizes",
    "block_quotients_complete", "smoothness", "source", "smooth_states",
    "construction",
  ];
  const out = Object.fromEntries(fields.filter((key) => Object.hasOwn(e, key)).map((key) => [key, e[key]]));
  out.smoothness = smoothness(e);
  out.source = source(e);
  const iso = e.isotriviality || {};
  if (iso.method === "genus_one_j_invariant") {
    out.isotriviality = { method: "genus_one_j_invariant" };
  } else if (iso.witnesses?.length) {
    out.isotriviality = { method: "proposition_1", witnesses: iso.witnesses.map((w) => Object.fromEntries(["block_quotient_index", "special_entry", "active_block", "claim"].filter((key) => Object.hasOwn(w, key)).map((key) => [key, w[key]]))) };
  } else {
    out.isotriviality = { method: "undetermined", ...(iso.failed_block_systems ? { failed_block_systems: iso.failed_block_systems } : {}) };
  }
  const quotientFields = ["block_count", "block_size", "blocks", "tuple", "genus", "indices", "ramification_excess"];
  out.block_quotients = (e.block_quotients || []).map((b) =>
    Object.fromEntries(quotientFields.filter((key) => Object.hasOwn(b, key)).map((key) => [key, b[key]])),
  );
  return out;
}
async function downloadFilteredFull() {
  const rows = selected.slice();
  if (!rows.length) return;
  if (rows.length > 25000) {
    info(`This filter matches ${number(rows.length)} examples. Add filters to narrow it below 25,000, or download the complete JSON.`);
    return;
  }
  const buckets = new Map();
  for (const base of rows) {
    const prefix = base.id[0];
    if (!buckets.has(prefix)) buckets.set(prefix, []);
    buckets.get(prefix).push(base);
  }
  const metadata = {
    schema_version: 6,
    dataset_version: dataset.dataset_version,
    record_view: "lean",
    last_updated: dataset.last_updated,
    conventions: dataset.conventions,
    criteria: {
      proposition_1: "The global manifest links the proof and defines Claims 6 and 7.",
    },
    selection: values(),
    component_count: rows.length,
    total_state_count: rows.reduce((sum, e) => sum + e.states, 0),
  };
  const details = new Map();
  try {
    $("download").disabled = true;
    for (const [prefix, bases] of buckets) {
      info(`Loading details: ${number(details.size)} of ${number(rows.length)}…`);
      const shard = await fetchDetailShard(prefix);
      for (const base of bases) {
        const detail = shard.get(base.id);
        if (!detail) throw new Error(`Missing detail record ${base.id}.`);
        details.set(base.id, detail);
      }
      shard.clear();
    }
    const parts = [JSON.stringify(metadata).slice(0, -1) + ',"components":['];
    let written = 0;
    for (const base of rows) {
      const full = leanRecord(expandDetail(base, details.get(base.id)));
      parts.push((written ? "," : "") + JSON.stringify(full));
      written++;
    }
    parts.push("]}\n");
    blobDownload(parts, `${BASE.slice(1)}-v${dataset.dataset_version}-filtered.json`);
    info(`Downloaded ${number(written)} lean examples.`);
  } catch (error) {
    info(`Could not build filtered export: ${error.message}`);
  } finally {
    $("download").disabled = selected.length === 0;
  }
}
async function loadDetail(base) {
  if (detailById.has(base.id)) return detailById.get(base.id);
  const url = new URL(base.detail_url || `/details/${base.id[0]}.jsonl.gz`, location.origin).href;
  let pending = detailShardCache.get(url);
  if (!pending) { pending = fetchDetailShard(base.id[0]); detailShardCache.set(url, pending); }
  const map = await pending, detail = map.get(base.id);
  if (!detail) throw new Error("This example is missing from its detail shard.");
  detailById.set(base.id, detail);
  return detail;
}
async function fetchDetailShard(prefix) {
  const response = await fetch(`/details/${prefix}.jsonl.gz?v=${encodeURIComponent(dataset.dataset_version)}`);
  if (!response.ok) throw new Error(`Detail request failed (${response.status}).`);
  const bytes = await response.arrayBuffer();
  let text;
  if (response.headers.get("Content-Encoding")?.includes("gzip")) text = await response.text();
  else {
    if (typeof DecompressionStream === "undefined") throw new Error("This browser cannot expand the compressed example details.");
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
    text = await new Response(stream).text();
  }
  const map = new Map();
  for (const line of text.split("\n")) if (line) { const value = JSON.parse(line); map.set(value.id, value); }
  return map;
}
function detailHTML(e) {
  const p = e.representative_g1_g2,
    idx = e.indices;
  const c12 = e.representative_pair_orbit_sizes.g1_g2.length,
    c34 = e.representative_pair_orbit_sizes.g3_g4.length;
  const proof = e.isotriviality.source ? `<a href="${escape(e.isotriviality.source)}">Proof · Proposition 1${e.isotriviality.witnesses?.length ? " · Claim " + e.isotriviality.witnesses[0].claim : ""} ↗</a>` : "";
  const onePiece = smoothness(e) === "one_piece";
  const additive = idx[0] + idx[1] === p.index ? `ind(g₁g₂) = ${p.index} = ${idx[0]} + ${idx[1]}`
    : idx[2] + idx[3] === p.index ? `ind(g₁g₂) = ${p.index} = ${idx[2]} + ${idx[3]}`
    : `ind(g₁g₂) = ${p.index}; neither ${idx[0]} + ${idx[1]} nor ${idx[2]} + ${idx[3]} (not smooth at this state)`;
  const certificate = onePiece
    ? `One-piece at every state: defect 0 and at most one positive-genus piece among the orbits of ⟨g₁,g₂⟩ and ⟨g₃,g₄⟩, certified by full braid closure. Index additivity holds at ${number(e.smooth_states)} of ${number(e.states)} states, so the component is not all-smooth.`
    : "Defect 0 and smoothness at every state are certified by full braid closure.";
  const origin = `${escape(SOURCE_LABELS[source(e)] || source(e))}${e.construction ? ` · ${escape(e.construction)}` : ""}`;
  return `<div class="detail"><div class="detail-heading"><h3>Degree ${e.degree} · genus ${e.genus} · ${number(e.states)} states</h3><button class="close-detail" data-hide="${e.id}" aria-label="Close example" title="Close example">×</button></div><div class="copy-actions"><button data-copy-tuple="${e.id}">Copy tuple</button><button data-copy-latex="${e.id}">Copy tuple LaTeX</button><button data-copy-json="${e.id}">Copy example JSON</button><button data-download-json="${e.id}">Download example JSON</button></div><p class="detail-copy-note hint" aria-live="polite"></p>
    <p class="small">${onePiece ? "One-piece (not all-smooth)" : "All-smooth"} · ${origin}</p><p class="small">Product g₁g₂g₃g₄ = 1 · transitive generated group · ${e.primitive ? "primitive" : "imprimitive"} · contained in ${escape(e.contained_in)} · group order ${number(e.group_order)}</p>
    <div class="tuple math">${e.tuple.map((t, k) => `<span>g<sub>${k + 1}</sub></span><p>= ${escape(t)}</p>`).join("")}</div>
    <div class="detail-grid"><div><h3>Representative product & smoothness</h3><p class="math">g₁g₂ = ${escape(p.permutation)}</p><p>Cycle count: ${p.cycle_count}, including fixed points. Index: ${p.index}.</p><p class="math">${additive}</p><p class="math">δ = ${p.cycle_count} − ${c12} − ${c34} + 1 = 0.</p><p class="small">${certificate} Product and pair orbits here are for the displayed tuple.</p></div>
    <div><h3>Representative pair orbits</h3><table class="mini"><thead><tr><th>Subgroup</th><th>Point-orbit sizes</th></tr></thead><tbody>${["g1_g2", "g2_g3", "g3_g4"].map((k, i) => `<tr><td>⟨g<sub>${i + 1}</sub>, g<sub>${i + 2}</sub>⟩</td><td class="math">${orbit(e.representative_pair_orbit_sizes[k])}</td></tr>`).join("")}</tbody></table><p class="small">Indices (${idx.join(", ")}). Fixed points count as cycles; multiplication acts rightmost first.</p></div></div>
    ${cuspSection(e)}
    <section class="proof-note"><strong>Isotrivial: ${e.isotrivial === null ? "unknown" : e.isotrivial ? "true" : "false"}</strong><p>${escape(e.isotriviality.explanation)}</p>${criterionFailures(e)}${proof}</section>
    <h3>${e.block_quotients_complete === false ? "Recorded invariant block-system witness" : "All recorded invariant block systems"} (${e.block_quotients.length})</h3><p class="small">The quotient genus belongs to a particular partition; the source genus stays ${e.genus}.${e.block_quotients_complete === false ? " This witness list is not a complete block-lattice enumeration." : ""}</p>
    ${e.block_quotients.length ? e.block_quotients.map((b, i) => quotientDetail(e, b, i)).join("") : "<p>No proper block system.</p>"}
    <p class="metadata">Entry revision ${e.revision} · updated ${escape(e.last_updated)}<br>Stable component ID: ${escape(e.id)}</p></div>`;
}
// Pieces graph at each cusp, computed here from the displayed tuple.
function cuspSection(e) {
  let r;
  try { r = cuspGraphs(e.tuple, e.degree); }
  catch (error) { return `<section class="cusps"><h3>Pieces at the three cusps</h3><p class="error">Could not draw the graph: ${escape(error.message)}</p></section>`; }
  const first = Math.max(0, r.cusps.findIndex((c) => !c.smoothSide));
  const uid = "cg-" + e.id.slice(0, 10);
  const list = (vs) => vs.map((v) => v.name).join(", ");
  const panels = r.cusps.map((c, i) => {
    const { svg, nodes } = cuspDrawing(c);
    const pos = nodes.filter((v) => v.genus > 0), brk = nodes.filter((v) => v.breaks);
    const sideName = (s) => (s === "L" ? c.leftName : c.rightName);
    const why = brk.map((v) => v.genus > 0 ? `${v.name} has genus ${v.genus}` : `${v.name} has genus 0 but meets ${v.degree} cycles`).join("; ");
    const smooth = c.smoothSide
      ? `<li class="ok"><strong>Smooth here.</strong> Every ${sideName(c.smoothSide)} piece is a genus-0 leaf (meets one cycle), so the index equality holds on that side.</li>`
      : `<li class="warn"><strong>Not smooth here.</strong> Neither side is all genus-0 leaves. On the ${sideName(c.breakSide)} side: ${why}.</li>`;
    const onePiece = c.positive === 0 ? "No positive-genus piece" : `Positive-genus piece${c.positive > 1 ? "s" : ""}: ${pos.map((v) => `${v.name} (genus ${v.genus})`).join(", ")}`;
    return `<div class="cg-panel" id="${uid}-${i}" role="tabpanel"${i === first ? "" : " hidden"}>
      <figure class="cg-figure"><div class="cg-scroll">${svg}</div>
        <figcaption class="cg-legend"><span><i class="cg-key pos"></i>positive genus</span><span><i class="cg-key brk"></i>breaks smoothness</span>${c.smoothSide ? `<span><i class="cg-key leaf"></i>genus-0 leaf side</span>` : ""}<span><i class="cg-key len">4</i>cycle length</span></figcaption></figure>
      <ul class="cg-facts">
        <li>${c.left.length} piece${c.left.length === 1 ? "" : "s"} of ${c.leftName}, ${c.right.length} of ${c.rightName}, joined by ${c.cycleCount === 1 ? "the single cycle" : `the ${c.cycleCount} cycles`} of ${c.name}.</li>
        <li class="math">δ = ${c.cycleCount} − ${c.left.length} − ${c.right.length} + 1 = ${c.defect}${c.defect === 0 ? ": a tree" : ""}</li>
        <li>${onePiece}${c.onePiece ? " → one-piece here." : "."}</li>
        ${smooth}
        <li class="math">g = Σ piece genera + δ = ${c.pieceGenus} + ${c.defect} = ${r.genus}</li>
      </ul></div>`;
  }).join("");
  const tabs = r.cusps.map((c, i) => `<button type="button" role="tab" class="cg-tab" data-cusp-tab="${i}" aria-controls="${uid}-${i}" aria-selected="${i === first}"><span class="math">${c.name}</span><small class="${c.smoothSide ? "ok" : "warn"}">${c.smoothSide ? "smooth" : "not smooth"}</small></button>`).join("");
  return `<section class="cusps"><h3>Pieces at the three cusps</h3>
    <p class="small">At a cusp two branch points collide. Left: the pieces (point orbits) of the colliding pair; right: the pieces of the other pair; each line is a cycle of their product, labelled by its length. A cusp is smooth when one side is all genus-0 leaves. Computed in your browser from the tuple above; the g₁g₃ cusp uses (g₁, g₃, g₃⁻¹g₂g₃, g₄). Rows marked ×k merge k identical leaves.</p>
    <div class="cg-tabs" role="tablist" aria-label="Cusp">${tabs}</div>${panels}</section>`;
}
async function openDetail(button, toggle = true) {
  const id = button.dataset.id,
    old = button.closest("tr").nextElementSibling;
  if (toggle && old?.classList.contains("detail-row")) {
    old.remove();
    button.setAttribute("aria-expanded", "false");
    button.textContent = "View tuple ↗";
    activeId = null;
    render();
    return;
  }
  document.querySelectorAll(".detail-row").forEach((r) => r.remove());
  document.querySelectorAll("[data-id]").forEach((b) => {
    b.setAttribute("aria-expanded", "false");
    b.textContent = "View tuple ↗";
  });
  const e = entries().find((e) => e.id === id);
  if (!e) return;
  const row = document.createElement("tr");
  row.className = "detail-row";
  row.innerHTML = '<td colspan="8"><p class="small">Loading example details…</p></td>';
  button.closest("tr").after(row);
  button.setAttribute("aria-expanded", "true");
  button.textContent = "Loading…";
  button.disabled = true;
  activeId = id;
  if (toggle)
    history.replaceState(
      null,
      "",
      `${location.pathname}?${filtersURL({ ...values(), example: id })}`,
    );
  try {
    const detail = await loadDetail(e);
    if (!row.isConnected) return;
    row.firstElementChild.innerHTML = detailHTML(expandDetail(e, detail));
    button.textContent = "Close ×";
  } catch (error) {
    if (!row.isConnected) return;
    row.firstElementChild.innerHTML = `<p class="error">${escape(error.message)}</p>`;
    button.textContent = "Retry details";
  } finally { button.disabled = false; }
}
form.addEventListener("input", () => {
  clearTimeout(timer);
  timer = setTimeout(() => {
    page = 0;
    activeId = null;
    render();
  }, 120);
});
form.addEventListener("submit", (e) => {
  e.preventDefault();
  clearTimeout(timer);
  page = 0;
  activeId = null;
  render();
});
$("reset").addEventListener("click", () => {
  setFilters({ isotrivial: DEFAULT_ISO });
  page = 0;
  activeId = null;
  render();
  info("Filters reset.");
});
$("share").addEventListener("click", () => {
  render();
  if ($("error").hidden) copy(location.href);
});
$("download").addEventListener("click", downloadFilteredFull);
for (const id of ["previous", "top-previous"])
  $(id).addEventListener("click", () => {
    page--;
    activeId = null;
    render();
    $("results").scrollIntoView({ block: "start" });
  });
for (const id of ["next", "top-next"])
  $(id).addEventListener("click", () => {
    page++;
    activeId = null;
    render();
    $("results").scrollIntoView({ block: "start" });
  });
$("mobile-sort").addEventListener("change", () => {
  form.elements.sort.value = $("mobile-sort").value;
  page = 0;
  activeId = null;
  render();
});
$("mobile-direction").addEventListener("click", () => {
  form.elements.direction.value =
    form.elements.direction.value === "asc" ? "desc" : "asc";
  page = 0;
  activeId = null;
  render();
});
document.addEventListener("click", async (ev) => {
  const target = ev.target.closest("button");
  if (!target) {
    const row = ev.target.closest(".data-row");
    if (row && !ev.target.closest("a,input,select,textarea")) {
      const b = row.querySelector("[data-id]");
      openDetail(b);
    }
    return;
  }
  if (target.dataset.cuspTab != null) {
    const box = target.closest(".cusps");
    box.querySelectorAll("[data-cusp-tab]").forEach((b) => b.setAttribute("aria-selected", String(b === target)));
    box.querySelectorAll(".cg-panel").forEach((p, i) => { p.hidden = String(i) !== target.dataset.cuspTab; });
    return;
  }
  if (target.dataset.sort) {
    const k = target.dataset.sort;
    form.elements.direction.value =
      form.elements.sort.value === k && form.elements.direction.value === "asc"
        ? "desc"
        : "asc";
    form.elements.sort.value = k;
    page = 0;
    activeId = null;
    render();
  }
  if (target.dataset.id) openDetail(target);
  if (target.dataset.copyJson) {
    const base=entries().find((e)=>e.id===target.dataset.copyJson);
    try { copy(JSON.stringify(leanRecord(expandDetail(base,await loadDetail(base))),null,2)); }
    catch(error) { info(error.message); }
  }
  if (target.dataset.downloadJson) {
    const base = entries().find((e) => e.id === target.dataset.downloadJson);
    try {
      const full = leanRecord(expandDetail(base, await loadDetail(base)));
      blobDownload(JSON.stringify({
        schema_version: 6,
        dataset_version: dataset.dataset_version,
        last_updated: dataset.last_updated,
        record_view: "lean",
        component: full,
      }, null, 2) + "\n", `smooth-${full.degree}-g${full.genus}-${full.id.slice(0, 12)}.json`);
      info("Downloaded the complete example JSON.");
    } catch (error) { info(error.message); }
  }
  if (target.dataset.hide) {
    const button = document.querySelector(`[data-id="${target.dataset.hide}"]`);
    if (button) {
      openDetail(button);
      document.querySelector(`[data-component="${target.dataset.hide}"]`)?.focus();
    }
  }
  if (target.dataset.copyTuple || target.dataset.copyLatex) {
    const base = entries().find(
      (e) => e.id === (target.dataset.copyTuple || target.dataset.copyLatex),
    );
    let e;
    try { e=expandDetail(base,await loadDetail(base)); }
    catch(error) { info(error.message); return; }
    const text = target.dataset.copyTuple
      ? e.tuple.map((p, i) => `g${i + 1} = ${p}`).join("\n")
      : "\\[\n\\begin{aligned}\n" +
        e.tuple
          .map(
            (p, i) =>
              `g_${i + 1} &= ${p.replace(/ /g, "\\,")}${i < 3 ? " \\\\" : ""}`,
          )
          .join("\n") +
        "\n\\end{aligned}\n\\]";
    copy(text);
  }
});
document.addEventListener("keydown", (ev) => {
  if (
    ev.target.classList.contains("data-row") &&
    ["Enter", " "].includes(ev.key)
  ) {
    ev.preventDefault();
    const b = ev.target.querySelector("[data-id]");
    openDetail(b);
  }
});
window.addEventListener("popstate", () => {
  const f = readFilters(location.search, DEFAULT_ISO);
  setFilters(f);
  activeId = f.example || null;
  page = 0;
  render(false);
});
async function fetchJSON(url, cache = "default") {
  const controller = new AbortController(),
    timeout = setTimeout(() => controller.abort(), 30000);
  try {
    const r = await fetch(url, { signal: controller.signal, cache });
    if (!r.ok) throw new Error(`Dataset request failed (${r.status}).`);
    if (url.endsWith('.gz') && !r.headers.get('Content-Encoding')) {
      if (typeof DecompressionStream === "undefined") throw new Error("This browser cannot expand the compressed dataset.");
      const stream = r.body.pipeThrough(new DecompressionStream("gzip"));
      return await new Response(stream).json();
    }
    return await r.json();
  } finally {
    clearTimeout(timeout);
  }
}
async function load() {
  try {
    const manifest = await fetchJSON("/manifest.json", "no-store"),
      last = manifest.versions.at(-1);
    dataset = await fetchJSON(last.snapshot.url);
    if (
      dataset.dataset_version !== manifest.latest_version ||
      dataset.components.length !== dataset.component_count
    )
      throw new Error("Dataset version/count mismatch. Please reload.");
    if (!dataset.components.every((e) => Object.hasOwn(e, "isotrivial") && e.classification_method))
      throw new Error(
        "This snapshot predates isotriviality classification. Please retry after the new release.",
      );
    countStats();
    const f = readFilters(location.search, DEFAULT_ISO);
    setFilters(f);
    activeId = f.example || null;
    // Shared example links reveal the actual matching page, not a hidden row.
    selected = filterEntries(entries(), values());
    const ix = selected.findIndex((e) => e.id === activeId);
    if (ix >= 0) page = Math.floor(ix / PAGE_SIZE);
    render();
  } catch (e) {
    $("stats").textContent = "Catalog could not be loaded.";
    $("result-count").textContent = "Unable to load";
    $("error").textContent = e.message + " The raw JSON remains linked above.";
    $("error").hidden = false;
  }
}
load();
