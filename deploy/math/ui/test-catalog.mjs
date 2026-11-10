import assert from "node:assert/strict";
import fs from "node:fs";
import {
  PAGE_SIZE,
  filterEntries,
  readFilters,
  filtersURL,
  advancedQuery,
  cycleText,
  cuspGraphs,
  cuspDrawing,
} from "./catalog.mjs";
assert.equal(PAGE_SIZE, 500);
const sample = (id, b) => ({
  id,
  degree: 24,
  genus: 4,
  states: 576,
  group_order: 98304,
  primitive: false,
  contained_in: "A24",
  isotrivial: false,
  tuple: ["(1 2)", "()", "()", "()"],
  cycle_types: [[2, 1, 1], [1], [1], [1]],
  block_quotients: b,
});
const a = sample("a", [
  { genus: 2, block_size: 2, ramification_excess: [2, 0, 0, 0] },
  { genus: 1, block_size: 4, ramification_excess: [2, 0, 2, 2] },
]);
assert.equal(
  filterEntries([a], { quotient_genus: "2", block_size: "4" }).length,
  0,
);
assert.equal(
  filterEntries([a], {
    quotient_genus: "2",
    block_size: "2",
    excess: "2,0,0,0",
  }).length,
  1,
);
assert.equal(
  filterEntries([a], { query: "degree >= 20 and states = 576" }).length,
  1,
);
assert.throws(() => advancedQuery("alert(1)"));
assert.throws(() => filterEntries([a], { degree_min: "25", degree_max: "20" }));
assert.throws(() => filterEntries([a], { degree_min: "1.5" }));
assert.throws(() => filterEntries([a], { excess: "2,0" }));
assert.equal(readFilters("").isotrivial, "non-isotrivial");
assert.equal(readFilters("?isotrivial=all").isotrivial, "all");
assert.equal(readFilters("", "all").isotrivial, "all");
assert.equal(readFilters("?isotrivial=unknown", "all").isotrivial, "unknown");
assert.equal(
  readFilters(
    "?" +
      filtersURL({
        query: "states >= 100 and genus < 9",
        isotrivial: "unknown",
      }),
  ).query,
  "states >= 100 and genus < 9",
);
assert.equal(cycleText([3, 3, 3, 1, 1]), "3^3 1^2");
const missing = {...a, id: "missing", group_order: null};
assert.equal(filterEntries([missing], {group_order_max: "1"}).length, 0);
assert.equal(filterEntries([missing], {query: "group_order != 1"}).length, 0);
assert.equal(filterEntries([missing], {}).length, 1);
const compact = {
  id: "compact", degree: 24, genus: 4, states: 576, group_order: null,
  primitive: false, contained_in: "A24", isotrivial: false,
  passport: ["2 1^2", "1", "1", "1"],
  block_summary: [[2, 2, 2, 0, 0, 0], [4, 1, 2, 0, 2, 2]],
};
assert.equal(filterEntries([compact], {
  search: "2 1^2", quotient_genus: "2", block_size: "2", excess: "2,0,0,0",
}).length, 1);
assert.equal(filterEntries([compact], {
  quotient_genus: "2", block_size: "4", excess: "2,0,0,0",
}).length, 0);
// v18: smoothness/source filters; old records without the fields are all-smooth.
const onePiece = {...compact, id: "op", isotrivial: true, primitive: true,
  smoothness: "one_piece", source: "group_search"};
assert.equal(filterEntries([compact, onePiece], {isotrivial: "all", smoothness: "all_smooth"})[0].id, "compact");
assert.equal(filterEntries([compact, onePiece], {isotrivial: "all", smoothness: "one_piece"})[0].id, "op");
assert.equal(filterEntries([compact, onePiece], {isotrivial: "all", smoothness: "all"}).length, 2);
assert.equal(filterEntries([compact, onePiece], {isotrivial: "all", source: "all_smooth_search"})[0].id, "compact");
assert.equal(filterEntries([compact, onePiece], {isotrivial: "all", source: "group_search", primitive: "true"})[0].id, "op");
assert.equal(filterEntries([compact, onePiece], {isotrivial: "all", source: "one_piece_lift"}).length, 0);
for (const direction of ["asc", "desc"]) {
  assert.equal(filterEntries([missing, a], {sort: "group_order", direction}).at(-1).id, "missing");
}
if (process.argv[2]) {
  const d = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
  const count = (value) => d.components.filter((e) => e.isotrivial === value).length;
  assert.equal(filterEntries(d.components, readFilters("")).length, count(false));
  assert.equal(
    filterEntries(d.components, { isotrivial: "isotrivial" }).length,
    count(true),
  );
  assert.equal(
    filterEntries(d.components, { isotrivial: "unknown" }).length,
    count(null),
  );
  assert.equal(filterEntries(d.components, { isotrivial: "all" }).length, d.component_count);
  assert.equal(filterEntries(d.components, { primitive: "true" }).length,
    d.components.filter((e) => e.primitive).length);
  const kinds = filterEntries(d.components, { isotrivial: "all", smoothness: "one_piece" }).length
    + filterEntries(d.components, { isotrivial: "all", smoothness: "all_smooth" }).length;
  assert.equal(kinds, d.component_count);
  const result = filterEntries(d.components, {
    isotrivial: "all",
    sort: "states",
    direction: "desc",
  });
  assert.ok(result.every((e, i) => !i || result[i - 1].states >= e.states));
}
console.log(
  "Catalog query tests passed (same-quotient conjunction, tri-state, safe parser, sorting, URL round trip).",
);

// Pieces graph: PSL(3,4) on 21 points, genus 1, one-piece, not smooth at the g1g2 cusp.
{
  const t = ["(1 11)(2 20)(3 9)(5 14)(7 21)(8 16)(10 13)(15 18)", "(1 20)(3 8)(4 12)(5 16)(6 17)(9 10)(13 14)(18 21)",
    "(2 21)(3 18)(4 20)(5 19)(10 16)(11 14)(12 15)(13 17)", "(1 14 10 11 4 15 2)(3 16 18 12 20 21 7)(5 17 6 13 9 8 19)"];
  const r = cuspGraphs(t, 21);
  assert.equal(r.genus, 1);
  const [c12, c23, c13] = r.cusps;
  assert.deepEqual([c12.left.length, c12.right.length, c12.cycleCount, c12.defect], [6, 2, 7, 0]);
  assert.deepEqual(c12.right.map((v) => [v.size, v.genus]).sort(), [[14, 1], [7, 0]].sort());
  assert.equal(c12.smoothSide, null);
  assert.deepEqual(c12.left.filter((v) => v.breaks).map((v) => [v.size, v.degree]), [[8, 2]]);
  assert.equal(c23.smoothSide, "L");
  assert.ok(c12.onePiece && c23.onePiece && c13.onePiece);
  const d = cuspDrawing(c12);
  assert.ok(d.svg.startsWith("<svg") && d.nodes.length === 8);
  assert.throws(() => cuspGraphs(["(1 2)", "(1 2)", "(1 2)", "()"], 2), /product is not 1/);
  console.log("Cusp graph tests passed.");
}
