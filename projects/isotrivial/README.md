# isotrivial: is the family of a braid orbit isotrivial?

Input: a product-one tuple `(g1, g2, g3, g4)` in `S_d` (cycle notation, separated by `;`). Products are composed
right to left. If the product is 1 only in the other order, the tuple is reversed (the same family, mirrored).

```sh
cargo build --release
./target/release/isotrivial "(1 2 4)(3 5 7)(6 8 9);(1 3 6);(1 2 7)(3 8 4)(5 9 6);(1 8 4)(2 9 3)(5 7 6)"
./target/release/isotrivial --selftest
```

## What it computes

1. The genus `g`, from `2g - 2 = -2d + ind g1 + ind g2 + ind g3 + ind g4`.
2. The braid orbit: breadth-first search with the six moves, each state stored in a canonical form up to
   simultaneous conjugation.
3. One-piece at every state: at the collision of `g1, g2`, the defect
   `cycles(g1 g2) - orbits<g1,g2> - orbits<g3,g4> + 1` is 0 and at most one piece has positive genus.
4. The Hodge degree per sheet, exactly:

   `12L = d - (rho(g1) + rho(g2) + rho(g3) + rho(g4)) + 3 * (average over the orbit of rho(g1 g2))`,

   where `rho(p)` is the sum over the cycles of `p` of `1/length`.
5. The verdict:
   - not one-piece: a limit curve is singular, so **not isotrivial**;
   - one-piece and `g >= 2`: **isotrivial** exactly when `L = 0` (otherwise `L > 0`);
   - one-piece and `g = 1`, or `g = 0`: **isotrivial**.

## Why it is correct

- One-piece means every limit curve of the family is smooth: at a collision the admissible-cover limit has the
  pieces as components and the cycles of `g1 g2` as nodes (Harris and Mumford), and its stable model (Deligne and
  Mumford) is smooth exactly when the defect is 0 and at most one piece has positive genus.
- For a family of smooth curves of genus at least 2 over a compact curve, the Hodge degree is at least 0, with
  equality exactly when the family is isotrivial: it is the degree of an ample class (Baily and Borel) pulled back
  along the period map, which is constant exactly when the curves are isomorphic (Torelli).
- The formula for `12L` is Mumford's relation `12 lambda = kappa + delta` on the admissible-cover family: `kappa`
  gives `d - S` and the nodes give the average. This follows the method of van der Geer and Kouvidakis,
  <https://arxiv.org/html/1010.3335v3>.

As a check, the program also averages `rho(g2 g3)` and `rho(g1 g2 g3 g2^-1)`. The three averages must agree; if
they do not, the report says so.

## Checks

- `--selftest`: three orbits (degree 9, 15 and 18) with known state counts and `12L = 0`, `16/15`, `2/3`.
- On 240 entries of the database (60 labelled isotrivial, 60 labelled non-isotrivial, 60 unknown, 60 one-piece but
  not all-smooth), the checker reproduced every state count and agreed with all 120 labelled entries.
- The largest all-smooth orbit in the database (degree 25, 447,120 states) takes about 15 seconds and 110 MB.

## Options

- `--json`: one JSON line per tuple.
- `--file FILE`: one tuple per line; lines starting with `#` are skipped.
- `--max-states N`: stop above `N` states (default 20,000,000).
