// Paired warm-query microbenchmark, not browser rendering or search throughput.
import fs from "node:fs";
import { performance } from "node:perf_hooks";
import assert from "node:assert/strict";
import { filterEntries } from "./catalog.mjs";
const code = fs.readFileSync(new URL("./catalog.mjs", import.meta.url), "utf8");
const original = code.replace(
  "!searchable(e).includes(text)",
  "!`${e.tuple.join(' ')} ${passportText(e)} ${e.contained_in}`.toLowerCase().replace(/\\s+/g, ' ').includes(text)",
);
assert.notEqual(original, code);
const before = await import(
  "data:text/javascript;base64," + Buffer.from(original).toString("base64")
);
const entries = JSON.parse(fs.readFileSync(process.argv[2], "utf8")).components;
const queries = [
  { isotrivial: "all", search: "3^3" },
  { isotrivial: "all", search: "(1 2)" },
  { isotrivial: "all", search: "1^" },
  { isotrivial: "all", search: "A25" },
];
for (const q of queries)
  assert.deepEqual(before.filterEntries(entries, q), filterEntries(entries, q));
const runs = [];
for (let i = 0; i < 31; i++) {
  const pair = {};
  for (const label of i % 2 ? ["after", "before"] : ["before", "after"]) {
    const start = performance.now();
    for (const q of queries)
      (label === "before" ? before.filterEntries : filterEntries)(entries, q);
    pair[label] = performance.now() - start;
  }
  runs.push(pair);
}
const median = (x) => [...x].sort((a, b) => a - b)[Math.floor(x.length / 2)];
console.log(
  JSON.stringify(
    {
      scope:
        "31 paired warm-query repeats, four textual filters across frozen 8095 records; Node only, no DOM/network",
      runs,
      median_before_ms: median(runs.map((x) => x.before)),
      median_after_ms: median(runs.map((x) => x.after)),
      median_paired_ratio: median(runs.map((x) => x.before / x.after)),
    },
    null,
    2,
  ),
);
