// Pure query functions: browser and Node regression tests use the same code.
export const PAGE_SIZE = 500;
export const fields = ["degree", "genus", "states", "group_order"];
export const status = (e) =>
  e.isotrivial === false
    ? "non-isotrivial"
    : e.isotrivial === true
      ? "isotrivial"
      : "unknown";
// Records published before v18 carry no smoothness/source field: they are all-smooth.
export const smoothness = (e) => e.smoothness || "all_smooth";
export const source = (e) => e.source || "all_smooth_search";
export const SOURCE_LABELS = {
  all_smooth_search: "All-smooth search",
  one_piece_census: "One-piece census (d ≤ 10)",
  one_piece_lift: "Lift construction",
  group_search: "Group search",
};
export const cycleText = (lengths) => {
  const counts = new Map();
  for (const n of lengths) counts.set(n, (counts.get(n) || 0) + 1);
  return [...counts].map(([n, m]) => `${n}${m > 1 ? "^" + m : ""}`).join(" ");
};
export const passportText = (e) =>
  Array.isArray(e.passport)
    ? e.passport.join("; ")
    : (e.cycle_types || []).map(cycleText).join("; ");
const searchCache = new WeakMap();
function searchable(e) {
  if (!searchCache.has(e))
    searchCache.set(
      e,
      `${(e.tuple || []).join(" ")} ${passportText(e)} ${e.contained_in}`
        .toLowerCase()
        .replace(/\s+/g, " "),
    );
  return searchCache.get(e);
}

export function numeric(value, name) {
  if (value === "" || value == null) return null;
  const n = Number(value);
  if (!Number.isSafeInteger(n) || n < 0)
    throw new Error(`${name} must be a nonnegative whole number.`);
  return n;
}
export function advancedQuery(text) {
  if (!text.trim()) return () => true;
  const clauses = text
    .trim()
    .split(/\s+and\s+/i)
    .map((part) => {
      const match = part.match(
        /^(degree|genus|states|group_order)\s*(>=|<=|!=|=|>|<)\s*(\d+)$/,
      );
      if (!match)
        throw new Error(
          "Use comparisons joined by “and”: degree >= 15 and genus <= 6. Fields: degree, genus, states, group_order.",
        );
      const [, field, op, raw] = match,
        n = numeric(raw, field);
      return (e) => typeof e[field] === "number" &&
        ({
          ">=": e[field] >= n,
          "<=": e[field] <= n,
          "!=": e[field] !== n,
          "=": e[field] === n,
          ">": e[field] > n,
          "<": e[field] < n,
        })[op];
    });
  return (e) => clauses.every((f) => f(e));
}
export function filterEntries(entries, f) {
  const range = {};
  for (const name of fields) {
    const low = numeric(f[name + "_min"], name),
      high = numeric(f[name + "_max"], name);
    if (low != null && high != null && low > high)
      throw new Error(`${name}: minimum exceeds maximum.`);
    range[name] = [low, high];
  }
  const qg = numeric(f.quotient_genus, "Quotient genus"),
    bs = numeric(f.block_size, "Block size");
  const excess = f.excess?.trim()
    ? f.excess
        .trim()
        .split(/[\s,]+/)
        .map((x) => numeric(x, "Excess"))
    : null;
  if (excess && excess.length !== 4)
    throw new Error("Ramification excess needs four integers, e.g. 2,0,0,0.");
  const test = advancedQuery(f.query || "");
  const text = (f.search || "").trim().toLowerCase().replace(/\s+/g, " ");
  const filtered = entries.filter((e) => {
    if (f.isotrivial && f.isotrivial !== "all" && status(e) !== f.isotrivial)
      return false;
    if (
      f.primitive &&
      f.primitive !== "all" &&
      String(e.primitive) !== f.primitive
    )
      return false;
    if (f.smoothness && f.smoothness !== "all" && smoothness(e) !== f.smoothness)
      return false;
    if (f.source && f.source !== "all" && source(e) !== f.source) return false;
    if (f.parity && f.parity !== "all" && e.contained_in[0] !== f.parity)
      return false;
    if (
      !fields.every(
        (k) =>
          (range[k][0] == null && range[k][1] == null) ||
          (typeof e[k] === "number" &&
          (range[k][0] == null || e[k] >= range[k][0]) &&
          (range[k][1] == null || e[k] <= range[k][1])),
      )
    )
      return false;
    if (text && !searchable(e).includes(text)) return false;
    // An existential conjunction on ONE quotient, never independent matches.
    if (
      (qg != null || bs != null || excess) &&
      !(e.block_summary || e.block_quotients || []).some((raw) => {
        const b = Array.isArray(raw)
          ? { block_size: raw[0], genus: raw[1], ramification_excess: raw.slice(2, 6) }
          : raw;
        return (qg == null || b.genus === qg) &&
          (bs == null || b.block_size === bs) &&
          (!excess || excess.every((x, i) => b.ramification_excess[i] === x));
      })
    )
      return false;
    return test(e);
  });
  const sort = fields.includes(f.sort) ? f.sort : "degree",
    dir = f.direction === "desc" ? -1 : 1;
  return filtered.sort(
    (a, b) =>
      ((a[sort] == null) - (b[sort] == null)) ||
      dir * ((a[sort] ?? 0) - (b[sort] ?? 0)) ||
      a.degree - b.degree ||
      a.genus - b.genus ||
      a.states - b.states ||
      a.id.localeCompare(b.id),
  );
}
export function readFilters(search, defaultIsotrivial = "non-isotrivial") {
  const p = new URLSearchParams(search),
    f = Object.fromEntries(p);
  if (!p.has("isotrivial")) f.isotrivial = defaultIsotrivial;
  return f;
}
export function filtersURL(f) {
  return new URLSearchParams(
    Object.entries(f).filter(([, v]) => v !== "" && v != null),
  ).toString();
}

// ---- Pieces graph at the three cusps (computed in the browser from the tuple) ----
// Vertices: point orbits of the colliding pair (L) and of the other pair (R).
// Edges: cycles of the colliding pair's product. Multiplication is rightmost first.
const SUB = "₀₁₂₃₄₅₆₇₈₉";
const sub = (n) => String(n).replace(/\d/g, (d) => SUB[d]);
export function parsePermutation(text, n) {
  const p = Array.from({ length: n }, (_, i) => i);
  for (const m of String(text).matchAll(/\(([^)]*)\)/g)) {
    const c = m[1].trim().split(/[\s,]+/).filter(Boolean).map((x) => Number(x) - 1);
    for (let i = 0; i < c.length; i++) p[c[i]] = c[(i + 1) % c.length];
  }
  return p;
}
const compose = (a, b) => b.map((x) => a[x]); // a∘b: apply b first
function inverse(a) { const r = new Array(a.length); a.forEach((y, x) => { r[y] = x; }); return r; }
function cyclesOf(p) {
  const seen = new Uint8Array(p.length), out = [];
  for (let s = 0; s < p.length; s++) {
    if (seen[s]) continue;
    const c = []; let x = s;
    while (!seen[x]) { seen[x] = 1; c.push(x); x = p[x]; }
    out.push(c);
  }
  return out;
}
function orbitLabels(gens, n) {
  const lab = new Int32Array(n).fill(-1), sizes = [];
  for (let s = 0; s < n; s++) {
    if (lab[s] >= 0) continue;
    const k = sizes.length, st = [s]; lab[s] = k; let size = 0;
    while (st.length) { const x = st.pop(); size++; for (const g of gens) { const y = g[x]; if (lab[y] < 0) { lab[y] = k; st.push(y); } } }
    sizes.push(size);
  }
  return { lab, sizes };
}
// Genus of each orbit of <a,b>, a three-point cover branched over a, b, (ab)^-1.
function pieceGenera(a, b, ab, lab, sizes) {
  const cyc = sizes.map(() => 0);
  for (const p of [a, b, ab]) for (const c of cyclesOf(p)) cyc[lab[c[0]]]++;
  // ind_a + ind_b + ind_ab = 3n - cycles; 2g - 2 = -2n + that.
  return sizes.map((n, k) => (n - cyc[k]) / 2 + 1);
}
const CUSPS = [
  { key: "12", name: "g₁g₂", left: "⟨g₁, g₂⟩", right: "⟨g₃, g₄⟩", pick: (g) => [g[0], g[1], g[2], g[3]] },
  { key: "23", name: "g₂g₃", left: "⟨g₂, g₃⟩", right: "⟨g₄, g₁⟩", pick: (g) => [g[1], g[2], g[3], g[0]] },
  { key: "13", name: "g₁g₃", left: "⟨g₁, g₃⟩", right: "⟨g₃⁻¹g₂g₃, g₄⟩",
    pick: (g) => [g[0], g[2], compose(compose(inverse(g[2]), g[1]), g[2]), g[3]] },
];
export function tupleGenus(perms) {
  const n = perms[0].length, ind = perms.reduce((s, p) => s + n - cyclesOf(p).length, 0);
  return ind / 2 - n + 1;
}
export function cuspGraph(perms, cusp) {
  const n = perms[0].length, [a, b, c, d] = cusp.pick(perms);
  const ab = compose(a, b), cd = compose(c, d);
  if (!compose(ab, cd).every((y, x) => y === x)) throw new Error("The tuple's product is not 1.");
  const L = orbitLabels([a, b], n), R = orbitLabels([c, d], n);
  const gL = pieceGenera(a, b, ab, L.lab, L.sizes), gR = pieceGenera(c, d, cd, R.lab, R.sizes);
  const cycles = cyclesOf(ab), edges = new Map();
  for (const cy of cycles) {
    const l = L.lab[cy[0]], r = R.lab[cy[0]];
    if (cy.some((x) => R.lab[x] !== r)) throw new Error("A product cycle leaves its piece.");
    const k = l + ":" + r;
    if (!edges.has(k)) edges.set(k, { l, r, lengths: [] });
    edges.get(k).lengths.push(cy.length);
  }
  const deg = (side, k) => cycles.reduce((s, cy) => s + ((side === "L" ? L : R).lab[cy[0]] === k ? 1 : 0), 0);
  const mk = (side, sizes, genera) => sizes.map((size, k) => ({ side, k, size, genus: genera[k], degree: deg(side, k) }));
  const left = mk("L", L.sizes, gL), right = mk("R", R.sizes, gR);
  const defect = cycles.length - left.length - right.length + 1;
  const isLeafSide = (vs) => vs.every((v) => v.genus === 0 && v.degree === 1);
  const offenders = (vs) => vs.filter((v) => v.genus > 0 || v.degree > 1);
  const smoothSide = isLeafSide(right) ? "R" : isLeafSide(left) ? "L" : null;
  let breakSide = null;
  if (!smoothSide) {
    const oL = offenders(left), oR = offenders(right), posL = oL.some((v) => v.genus > 0), posR = oR.some((v) => v.genus > 0);
    breakSide = oL.length < oR.length || (oL.length === oR.length && !posL && posR) ? "L" : "R";
    for (const v of breakSide === "L" ? oL : oR) v.breaks = true;
  }
  const positive = [...left, ...right].filter((v) => v.genus > 0);
  return {
    key: cusp.key, name: cusp.name, leftName: cusp.left, rightName: cusp.right,
    left, right, edges: [...edges.values()], cycleCount: cycles.length, defect,
    pieceGenus: positive.reduce((s, v) => s + v.genus, 0), positive: positive.length,
    onePiece: defect === 0 && positive.length <= 1, smoothSide, breakSide,
  };
}
export function cuspGraphs(tuple, degree) {
  const perms = tuple.map((t) => parsePermutation(t, degree));
  return { genus: tupleGenus(perms), cusps: CUSPS.map((c) => cuspGraph(perms, c)) };
}

// Layout: leaves that look the same (side, genus, size, neighbour, cycle length) are
// merged when a side is crowded; rows are ordered by the barycentre heuristic.
function layoutGraph(G) {
  const nodes = [...G.left, ...G.right].map((v) => ({ ...v, count: 1, id: v.side + v.k }));
  const byId = new Map(nodes.map((v) => [v.id, v]));
  let links = G.edges.map((e) => ({ a: "L" + e.l, b: "R" + e.r, lengths: e.lengths }));
  for (const side of ["L", "R"]) {
    const vs = nodes.filter((v) => v.side === side);
    if (vs.length <= 6) continue;
    const groups = new Map();
    for (const v of vs) {
      if (v.degree !== 1 || v.breaks) continue;
      const ln = links.find((x) => x.a === v.id || x.b === v.id);
      const other = ln.a === v.id ? ln.b : ln.a;
      const key = [v.genus, v.size, other, ln.lengths.join(",")].join("|");
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key).push(v);
    }
    for (const g of groups.values()) {
      if (g.length < 2) continue;
      g[0].count = g.length;
      for (const v of g.slice(1)) { byId.delete(v.id); links = links.filter((x) => x.a !== v.id && x.b !== v.id); }
    }
  }
  const col = (side) => [...byId.values()].filter((v) => v.side === side);
  let left = col("L").sort((x, y) => y.genus - x.genus || y.degree - x.degree || y.size - x.size);
  let right = col("R");
  const nbr = new Map([...byId.keys()].map((id) => [id, []]));
  for (const x of links) { nbr.get(x.a).push(x.b); nbr.get(x.b).push(x.a); }
  const pos = new Map();
  const place = (vs) => vs.forEach((v, i) => pos.set(v.id, (i + 0.5) / vs.length));
  const bary = (v) => { const ns = nbr.get(v.id); return ns.length ? ns.reduce((s, u) => s + pos.get(u), 0) / ns.length : 0.5; };
  place(left);
  for (let it = 0; it < 6; it++) {
    right = right.map((v, i) => [bary(v), i, v]).sort((x, y) => x[0] - y[0] || x[1] - y[1]).map((x) => x[2]); place(right);
    left = left.map((v, i) => [bary(v), i, v]).sort((x, y) => x[0] - y[0] || x[1] - y[1]).map((x) => x[2]); place(left);
  }
  const label = (vs, letter) => { let i = 1; for (const v of vs) { v.name = v.count > 1 ? `${letter}${sub(i)}–${letter}${sub(i + v.count - 1)}` : letter + sub(i); i += v.count; } };
  label(left, "L"); label(right, "R");
  for (const x of links) x.count = byId.get(x.a).count * byId.get(x.b).count;
  return { left, right, links, byId };
}
function svgText(s) { return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]); }
export const cuspSVG = (G) => cuspDrawing(G).svg;
export function cuspDrawing(G) {
  const { left, right, links, byId } = layoutGraph(G);
  const W = 360, NW = 112, NH = 40, ROW = 52, TOP = 30;
  const rows = Math.max(left.length, right.length, 1), H = TOP + rows * ROW + 6;
  const xL = 6 + NW / 2, xR = W - 6 - NW / 2;
  const yOf = (vs) => new Map(vs.map((v, i) => [v.id, TOP + ((i + 0.5) * rows * ROW) / vs.length]));
  const y = new Map([...yOf(left), ...yOf(right)]);
  const bez = (t, p0, p1, p2, p3) => (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3;
  const placed = [], edgeSVG = [], labelSVG = [];
  const x0 = xL + NW / 2, x3 = xR - NW / 2, xm = (x0 + x3) / 2;
  for (const ln of links) {
    const y0 = y.get(ln.a), y3 = y.get(ln.b);
    const multi = ln.lengths.length > 1, thick = ln.count > 1 || multi;
    edgeSVG.push(`<path class="cg-edge${thick ? " cg-multi" : ""}" d="M${x0} ${y0} C${xm} ${y0} ${xm} ${y3} ${x3} ${y3}"/>`);
    const text = (multi ? cycleText(ln.lengths).replace(/\^(\d+)/g, (_, m) => "×" + m) : String(ln.lengths[0])) + (ln.count > 1 ? " ×" + ln.count : "");
    const w = 8 + 6.4 * text.length;
    let best = null;
    // Labels go where the lines fan out: near the side with more rows.
    const ts = left.length > right.length ? [0.3, 0.42, 0.2, 0.52, 0.62] : left.length < right.length ? [0.7, 0.58, 0.8, 0.48, 0.38] : [0.5, 0.36, 0.64, 0.26, 0.74];
    for (const t of ts) {
      const px = bez(t, x0, xm, xm, x3), py = bez(t, y0, y0, y3, y3);
      const hit = placed.some((q) => Math.abs(q.x - px) < (q.w + w) / 2 + 2 && Math.abs(q.y - py) < 15);
      if (!hit) { best = { x: px, y: py, w }; break; }
      best ||= { x: px, y: py, w };
    }
    placed.push(best);
    labelSVG.push(`<g class="cg-len"><rect x="${(best.x - best.w / 2).toFixed(1)}" y="${(best.y - 8).toFixed(1)}" width="${best.w.toFixed(1)}" height="16" rx="8"/><text x="${best.x.toFixed(1)}" y="${(best.y + 4).toFixed(1)}">${svgText(text)}</text></g>`);
  }
  const node = (v, cx) => {
    const cy = y.get(v.id), cls = ["cg-node", v.genus > 0 ? "cg-pos" : "", v.breaks ? "cg-break" : "", G.smoothSide === v.side ? "cg-leaf" : ""].join(" ").trim();
    const title = `${v.name}: ${v.count > 1 ? v.count + " pieces, each " : ""}${v.size} points, genus ${v.genus}, meets ${v.degree} cycle${v.degree === 1 ? "" : "s"} of the product${v.breaks ? " — breaks smoothness here" : ""}`;
    return `<g class="${cls}"><title>${svgText(title)}</title><rect x="${cx - NW / 2}" y="${cy - NH / 2}" width="${NW}" height="${NH}" rx="9"/>` +
      `<text class="cg-name" x="${cx}" y="${cy - 3}">${svgText(v.name)}${v.count > 1 ? `<tspan class="cg-count"> ×${v.count}</tspan>` : ""}</text>` +
      `<text class="cg-sub" x="${cx}" y="${cy + 12}">${v.size} pt${v.size === 1 ? "" : "s"} · genus ${v.genus}</text></g>`;
  };
  const svg = `<svg class="cg-svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="Pieces graph at cusp ${svgText(G.name)}: ${G.left.length} pieces on the left, ${G.right.length} on the right, ${G.cycleCount} product cycles.">` +
    `<text class="cg-head" x="${xL}" y="14">${svgText(G.leftName)}</text><text class="cg-head" x="${xR}" y="14">${svgText(G.rightName)}</text>` +
    edgeSVG.join("") + labelSVG.join("") + left.map((v) => node(v, xL)).join("") + right.map((v) => node(v, xR)).join("") + `</svg>`;
  return { svg, nodes: [...left, ...right] };
}
