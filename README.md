# braid-orbits

Code for studying **braid orbits of four-tuples of permutations**: product-one generating tuples
`(g1, g2, g3, g4)` in `S_d` or `A_d`. Each orbit is a one-parameter family of curves. The code finds and
checks orbits whose curves stay smooth everywhere in the family (*all-smooth* and *one-piece* orbits). It also
decides whether the family is isotrivial, and computes braid-graph spectra.

All examples found so far are on the website, where you can search and filter them:
<https://math.xupingz.com/all-smooth> and <https://math.xupingz.com/all-one-piece>.

## Start here

| I want to... | go to |
|---|---|
| look at examples | the website above. Nothing to install. |
| download every example | [The database](#the-database) |
| understand the words used here | [The mathematics in one page](#the-mathematics-in-one-page) |
| run the code on my own tuple | [Setup](#setup), then [Try it on your own tuple](#try-it-on-your-own-tuple) |
| understand one particular search | [What is in each folder](#what-is-in-each-folder). Each folder has its own README. |

## Contents

1. [The mathematics in one page](#the-mathematics-in-one-page)
2. [The database](#the-database)
3. [Setup](#setup)
4. [Try it on your own tuple](#try-it-on-your-own-tuple)
5. [What is in each folder](#what-is-in-each-folder)
6. [What the searches found](#what-the-searches-found)
7. [How to read the results](#how-to-read-the-results)
8. [Running larger searches](#running-larger-searches)
9. [Notes for programmers](#notes-for-programmers)

## The mathematics in one page

**Tuples.** A tuple is four permutations `g1, g2, g3, g4` of `{1, ..., d}` with `g1 g2 g3 g4 = 1` that together
move every point (the group they generate is transitive). It describes a degree-`d` cover of the sphere branched
over four points. The *index* of a permutation is `ind g = d - (number of cycles, counting fixed points)`. The
genus of the cover comes from Riemann-Hurwitz: `2g - 2 = -2d + ind g1 + ind g2 + ind g3 + ind g4`. The *passport*
is the list of the four cycle types, such as `3^2 2 1^4`, which means two 3-cycles, one 2-cycle and four fixed
points.

**Braid moves and orbits.** Moving one branch point around another changes the tuple by a braid move. The move
`sigma_j` replaces the neighbours `(a, b)` in positions `j, j+1` by `(a b a^-1, a)`. The *braid orbit* is every tuple
you can reach by braid moves. Each tuple is taken up to relabelling the points `1..d` (simultaneous conjugation by
`S_d`), and one tuple in an orbit is called a *state*. An orbit is one connected component of the Hurwitz space: a
one-parameter family of curves.

**Pieces and defect.** When branch points 1 and 2 collide, the curve breaks into pieces:
- the left pieces are the orbits of `<g1, g2>`;
- the right pieces are the orbits of `<g3, g4>`;
- the pieces are joined by the cycles of `g1 g2`.

These form a graph, which the website draws for every example. Its *defect* is its number of independent loops,
`cycles(g1 g2) - orbits<g1,g2> - orbits<g3,g4> + 1`.

**All-smooth.** At every state, one side of this graph is all genus-0 leaves. Equivalently,
`ind(g1 g2) = ind g1 + ind g2` or `ind(g3 g4) = ind g3 + ind g4`.

**One-piece.** At every state, the defect is 0 and at most one piece has positive genus. This means the limit curve
is still smooth (`projects/geosmooth/README.md`, Lemma G1). All-smooth implies one-piece.

**Isotrivial.** All curves in the family are the same curve. The central question is whether a *non*-isotrivial
all-smooth family exists in genus 3.

**Primitive.** No proper block system: there is no way to split `{1..d}` into equal blocks that the four
permutations only permute among themselves. Imprimitive orbits come from a smaller cover, so primitive ones are
the most interesting.

**Conventions used everywhere.**
- Permutations are written in cycle notation, starting from 1.
- Products are composed **rightmost first**: `(x y)(i) = x(y(i))`.
- "Orbit size" or "states" counts tuples up to relabelling by `S_d`, unless a tool says otherwise.

## The database

Every certified example is in one of two files. Each file is the complete database, not a sample.

| what | download link | download size | after unzipping | examples |
|---|---|---|---|---|
| All-smooth | <https://math.xupingz.com/all-smooth.json.gz> | 67.4 MB | 690.2 MB | 256,560 |
| All one-piece (includes all the all-smooth ones) | <https://math.xupingz.com/all-one-piece.json.gz> | 94.4 MB | 957.5 MB | 328,033 |

Both files are compressed (`.gz`), so they download about ten times faster. On a Mac, double-click one to unzip it.
On Linux, run `gunzip all-smooth.json.gz`. Or read it directly with Python as shown below, without unzipping. You do
not need this repository to use the database.

Each file is one JSON object. Its `conventions` field explains every field, and `components` is the list of
examples, one per braid orbit. Each example has fields such as:
- `degree`, `genus`;
- `states`, the orbit size;
- `passport`;
- `primitive`;
- `isotrivial`;
- `tuple`, one state of the orbit.

To load it in Python. This takes under a minute and about 5 GB of memory:

```python
import gzip, json
with gzip.open("all-one-piece.json.gz", "rt") as f:   # or "all-smooth.json.gz"
    db = json.load(f)
examples = db["components"]
print(len(examples), "examples")
first = examples[0]
print(first["degree"], first["genus"], first["passport"])
print(first["tuple"])
# for example, all examples of degree at most 8, and how many are isotrivial:
small = [e for e in examples if e["degree"] <= 8]
print(len(small), sum(e["isotrivial"] is True for e in small))
```

`isotrivial` is `true`, `false`, or `null`, which means the test used could not decide. In Python these become
`True`, `False` and `None`.

## Setup

You need Rust, Python 3.11 or newer, `uv` (fetches SymPy and NumPy for some tests) and Node 18 or newer (website
tests only). Setup takes about 20 minutes the first time, mostly downloads.

### macOS

```sh
xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"
brew install python node uv
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
```

`xcode-select --install` opens a window; let it finish before the next line. On an Intel Mac, Homebrew lives in
`/usr/local` instead of `/opt/homebrew`; change the two `brew shellenv` lines accordingly.

### Windows

The code uses Unix features (`make`, `bash`, `fcntl`, `/dev/stdout`), so it runs inside WSL (Linux on Windows).
In PowerShell as administrator:

```powershell
wsl --install -d Ubuntu-24.04
```

Restart, open the "Ubuntu" app, create a Linux user when asked, then follow the Linux steps below inside it. Keep
the repository inside the Linux home folder (for example `~/braid-orbits`), not under `/mnt/c`, which is much
slower.

### Linux (Ubuntu or Debian, and WSL)

```sh
sudo apt update
sudo apt install -y build-essential git curl python3 nodejs
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
curl -LsSf https://astral.sh/uv/install.sh | sh
source "$HOME/.local/bin/env"
```

Ubuntu 24.04 ships Python 3.12 and Node 18. On older releases, run `uv python install 3.12` and use `python3.12`
in place of `python3`; Node is only needed for the website tests.

### Check, build, test

```sh
cargo --version && python3 --version && uv --version && node --version
cd braid-orbits
make build
make test
```

`make build` takes 1 to 5 minutes the first time (50 seconds on a recent Mac). `make test` is optional and takes
5 to 15 minutes; it ends with `All Python tests passed.` and the website test results. A few tests are skipped on
purpose; see [Notes for programmers](#notes-for-programmers).

**If something fails:**
- `command not found` right after installing: open a new terminal, or run the `source` lines above again.
- `linker 'cc' not found`: the C toolchain is missing (`xcode-select --install` on macOS,
  `sudo apt install build-essential` on Linux).
- `make test` fails at a SymPy or NumPy test: `uv` could not download the package. Check the connection and rerun.

## Try it on your own tuple

Every command here runs from the `braid-orbits` folder and prints its answer in the terminal. All the tuples
below are real entries in the database. Everything here except the degree-9 search takes under a second.

### 1. Close one orbit

`onepiece --orbit` takes a tuple, computes its whole braid orbit, and checks every state. Put the tuple between
the quotes: four permutations in cycle notation, separated by `;`, with product 1.

**The smallest all-smooth family:** degree 9, genus 2.

```sh
projects/onepiece/target/release/onepiece --orbit "(5 6 7);(1 2 5)(3 7 9)(4 6 8);(1 3 7)(2 6 9)(4 5 8);(1 4 5)(2 7 8)(3 6 9)"
```

It prints one long line. The fields to look at first:

| field | here | meaning |
|---|---|---|
| `d`, `genus` | 9, 2 | degree and genus of the cover |
| `states` | 12 | orbit size, up to relabelling |
| `all_smooth` | true | all-smooth at every state |
| `all_one_piece` | true | one-piece at every state |
| `all_zero` | true | defect 0 at every state |
| `comb_smooth_states` | 12 | states where the index equality holds |
| `bad_clusters` | 0 | places in the family where the curve becomes singular (0 means none) |
| `primitive`, `block_sizes` | false, [3] | imprimitive, with blocks of size 3 |
| `types`, `tuple` | | the passport and the tuple in a standard form |

**The smallest-degree non-isotrivial all-smooth family:** degree 15, genus 5. It prints `"states":540` and
`"all_smooth":true`; the database marks it non-isotrivial.

```sh
projects/onepiece/target/release/onepiece --orbit "(1 2 3);(1 4 7 10 13)(2 5 8 11 14)(3 6 9 12 15);(1 5 9 12 15)(2 6 8 11 14)(3 4 7 10 13);(1 12 5 14 8)(2 10 4 15 9)(3 11 6 13 7)"
```

**One-piece but not all-smooth:** degree 6, genus 1. It prints `"states":24`, `"all_one_piece":true` and
`"all_smooth":false`.

```sh
projects/onepiece/target/release/onepiece --orbit "(1 2 3 6 5 4);(1 3 5);(1 3 5);(1 4)(2 3)(5 6)"
```

Add `--clusters` to any of these to see the orbit cluster by cluster. `projects/onepiece/README.md` explains the
other fields. `onepiece` handles tuples on at most 16 points (24 if built with
`cargo build --release --features wide` inside `projects/onepiece`).

### 2. Search a whole degree

`onepiece --census` finds **every** orbit of a given degree and genus that is one-piece at every state. It prints
one line per orbit found, then a summary line (`"summary":true`, with `"orbits_reported"` and `"complete":true`).

**Degree 6, genus 1:** there are exactly 2, found instantly. One of them is the degree-6 example above.

```sh
projects/onepiece/target/release/onepiece --census --d 6 --genus 1
```

**Degree 9, genus 2:** there is exactly one, and it is the smallest all-smooth family from section 1. This takes
about 5 seconds on a recent Mac.

```sh
projects/onepiece/target/release/onepiece --census --d 9 --genus 2
```

**Every orbit, smooth or not:** add `--all-orbits`. In degree 5 and genus 2 this lists all 38 orbits; none is
one-piece everywhere, and the fewest bad clusters any of them has is 1.

```sh
projects/onepiece/target/release/onepiece --census --d 5 --genus 2 --all-orbits
```

Larger degrees take much longer; `projects/onepiece/README.md` has the timings. For a long run, save the results
to a file with `--out FILE`. The program then also shows its progress.

## What is in each folder

```
braid-orbits/
├── README.md          this file
├── Makefile           make build, make test
├── src/               main Rust engine: hurwitz_breaker, plus 24 programs in src/bin/
├── scripts/           Python: launchers, independent checks, catalogue builders (6 topic folders)
├── projects/          separate, self-contained studies, each with its own README
├── deploy/math/ui/    the website (plain HTML and JavaScript)
├── evidence/          80 saved input files that the code and tests read
└── data/              created when you run a search; holds its output
```

**Projects.** Each one answers a specific question, and its README states the question, the exact set of tuples
searched, and what was found. Start with these.

| folder | question | language |
|---|---|---|
| `projects/onepiece/` | Which orbits are one-piece everywhere? Complete censuses for small degrees, lifts, and searches over primitive groups such as PSL(3,q). `allzero/` holds the primitive one-piece examples and their independent checks. | Rust, Python |
| `projects/geosmooth/` | How do geometric smoothness (one-piece) and combinatorial smoothness (all-smooth) differ? A graph criterion for one-piece. | Rust, Python |
| `projects/frames/` | One branch point moving over a fixed elliptic curve: which collisions occur. | Rust |
| `projects/quotient_exhaust/` | Every genus-3 tuple over a fixed genus-0 or genus-1 quotient, enumerated exhaustively. | Rust |
| `projects/genus4_structure/` | Helpers for block systems and quotients, used by the geosmooth checks. | Python |
| `projects/common/` | `braid.py`: shared Python helpers (parse cycle notation, compose, braid moves, close an orbit). | Python |

**Scripts** (`scripts/README.md` gives an overview, and each topic folder has a table of every file):

| folder | contents |
|---|---|
| `scripts/smooth/` | all-smooth searches, small-degree censuses, and audits that re-check every result |
| `scripts/genus3/` | the search for a non-isotrivial genus-3 all-smooth orbit, and its finite exhaustive slices |
| `scripts/weak/` | the earlier project on weak pair-generation failures, and the defect questions |
| `scripts/spectra/` | braid-graph spectra and defect statistics |
| `scripts/catalog/` | builds the website's database: isotriviality test, stable IDs, versioned releases |
| `scripts/gpu/` | an Apple GPU prototype. It was slower than the CPU code, so nothing uses it. |

Script names follow one pattern:
- `run_*` starts a search;
- `verify_*`, `audit_*`, `check_*` re-check results with separately written code;
- `test_X.py` tests `X.py`.

## What the searches found

The details and exact scopes are in each project's README and on the website.

- **Catalogue.** 256,560 all-smooth orbits and 71,473 more that are one-piece but not all-smooth (328,033 in
  total), with tuples, passports, block systems and isotriviality where it is decided.
- **Genus 2.** Every all-smooth family is isotrivial, in every degree, because `M_2` contains no complete curves.
- **Genus 4.** Non-isotrivial all-smooth families exist; the catalogue contains them, along with examples in
  higher genus.
- **Genus 3.** Still open: no non-isotrivial all-smooth orbit is known. Many families and degrees have been
  excluded exhaustively (`scripts/genus3/README.md`).
- **Primitive orbits.** Primitive one-piece orbits exist. PSL(3,4) acting on the 21 points of PG(2,4) gives
  genus-1 ones, and PSL(3,5) on the 31 points of PG(2,5) gives genus-3 ones. There are 29 in all, in genus 1 and
  3 to 7, from projective planes PG(2,q) and a 5 x 5 product action; all are isotrivial
  (`projects/onepiece/allzero/README.md`). No primitive one-piece orbit of genus 2 was found. A primitive
  all-smooth orbit is yet to be found: all 256,560 all-smooth orbits in the catalogue are imprimitive.

## How to read the results

- **Every example is checked twice.** Each positive example was re-verified by separately written code.
- **A cap is not a result.** A search that stopped at a time or size limit proves nothing about that case.
- **Samples are not censuses.** Random samples are not distinct orbits. An exhausted slice (one passport, one
  quotient, one degree) is not a whole-degree census. Each README states its exact population.

## Running larger searches

- Run every command from the `braid-orbits` folder. Paths such as `data/...`, `src/...` and `target/...` are
  relative to it.
- Searches write to `data/` here. Git ignores it. To keep large output on another disk, make `data` a link to a
  folder there:

  ```sh
  ln -s /path/to/big/disk/braid-orbits-data data
  ```

- Most launchers print their options with `--help`, for example `python3 scripts/smooth/run_forest_fast.py --help`.
  Some can resume an earlier run and name one by default (for example `data/smooth_frontier_2026...`). Those runs are
  not included here, so pass another folder or start a new run.
- Long searches can use every CPU core for hours or days. Check a launcher's README and `--help` for its stop
  command before starting it.

## Notes for programmers

**Build details.**
- `make build` builds the main crate into `target/smart/` with `-C target-cpu=native`, and the four Rust projects
  into their own `target/` folders.
- There are no Rust dependencies.
- Canonical state keys support degree at most 25 by default. Build with `RUSTFLAGS="--cfg smooth_wide_keys"` for
  degrees up to 63.

**Single tools.**

```sh
./target/smart/release/hurwitz_breaker selftest
projects/onepiece/target/release/onepiece --selftest
python3 scripts/weak/run_campaign.py --help
node deploy/math/ui/test-catalog.mjs
```

**Python imports.** Scripts import each other by module name across the six topic folders through
`_scripts_path.py`. The details are in `scripts/README.md`.

**Tests skipped by `make test`.** It runs everything that works from this repository alone, and skips:
- One Rust test, `triangle_draws_have_the_declared_geodesies_and_budget`. It also fails in the original project.
  The TRIANGLE sampler in `src/smooth_structural_novelty.rs` was rewritten on Oct 31 to apply index bounds while
  listing passports, and this older test still expects the earlier draws.
- Eight Python tests that need files not included here:
  - These need historical runs or website exports: `test_defect_coin_model`, `test_genus_three_double_marks`,
    `test_genus_three_frame_retirement`, `test_genus_three_saturated_launch`, `test_mathematical_patterns`,
    `test_negative_audit`.
  - These need the proof PDF that every website release pins by hash: `test_build_public_smooth_catalog`,
    `test_package_current_math_release`. For the same reason, building a new website release needs
    `deploy/math/sources/isotrivial-proof.pdf`, which you can download from the live site at
    `/sources/isotrivial-proof.pdf`.

**Do not edit `evidence/`.** The code checks several of those files by SHA-256 hash.

**What is not here.** These stayed in the original project:
- the written notes (PDF and TeX);
- the raw run output, hundreds of GB;
- the per-run frozen copies of the source;
- the per-result check scripts.

Ask if you need any of them.
