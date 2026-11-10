//! Complete (4 2 1^9),(3^5)^3 target, modulo simultaneous S15 conjugation.
//! No geometric conjecture or Hodge filter is used to reject a state.
//! Completion propagation enforces c^3=d^3=1 and cd=p^-1 simultaneously.
#[path = "genus_three_factors.rs"]
mod factors;
use std::{
    collections::{HashMap, VecDeque},
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    time::Instant,
};
const N: usize = 15;
const NIL: u8 = 255;
type P = [u8; N];
type T = [P; 4];
const WORDS: [&[u8]; 6] = [&[0, 0], &[1], &[2], &[3, 3], &[4], &[5]];
fn id() -> P {
    std::array::from_fn(|i| i as u8)
}
fn mul(a: P, b: P) -> P {
    std::array::from_fn(|i| a[b[i] as usize])
}
fn inv(p: P) -> P {
    let mut q = id();
    for i in 0..N {
        q[p[i] as usize] = i as u8;
    }
    q
}
fn cycle(p: &mut P, c: &[usize]) {
    for i in 0..c.len() {
        p[c[i]] = c[(i + 1) % c.len()] as u8;
    }
}
fn pack(p: P) -> u64 {
    p.iter()
        .enumerate()
        .fold(0, |s, (i, &x)| s | (x as u64) << (4 * i))
}
fn braid(t: &mut T, m: u8) {
    let i = m as usize % 3;
    let (a, b) = (t[i], t[i + 1]);
    if m < 3 {
        t[i] = mul(mul(a, b), inv(a));
        t[i + 1] = a;
    } else {
        t[i] = b;
        t[i + 1] = mul(mul(inv(b), a), b);
    }
}
fn apply(t: &mut T, word: &[u8]) {
    for &m in word {
        braid(t, m);
    }
}
fn full(p: P) -> bool {
    let mut x = 0;
    for _ in 0..N - 1 {
        x = p[x] as usize;
        if x == 0 {
            return false;
        }
    }
    p[x] == 0
}
fn order_three(p: P, n: usize) -> bool {
    (0..n).all(|i| p[i] as usize != i && p[p[p[i] as usize] as usize] as usize == i)
}

#[derive(Clone, Copy)]
struct Partial {
    image: [[u8; N]; 2],
    inverse: [[u8; N]; 2],
}
impl Partial {
    fn empty() -> Self {
        Self {
            image: [[NIL; N]; 2],
            inverse: [[NIL; N]; 2],
        }
    }
    // Every successful assignment consumes an unset edge; queue depth <= 2n.
    // A length-two path forces its closing edge in an order-three permutation.
    fn force(&mut self, side: usize, u: usize, v: usize, n: usize) -> bool {
        if self.image[side][u] != NIL {
            return self.image[side][u] as usize == v;
        }
        if u == v || self.inverse[side][v] != NIL {
            return false;
        }
        self.image[side][u] = v as u8;
        self.inverse[side][v] = u as u8;
        // c*d=R, R(i)=i-1. c(u)=v iff d(v+1)=u.
        let paired = if side == 0 {
            self.force(1, (v + 1) % n, u, n)
        } else {
            self.force(0, v, (u + n - 1) % n, n)
        };
        if !paired {
            return false;
        }
        let next = self.image[side][v];
        if next != NIL && !self.force(side, next as usize, u, n) {
            return false;
        }
        let prev = self.inverse[side][u];
        if prev != NIL && !self.force(side, v, prev as usize, n) {
            return false;
        }
        true
    }
}
#[derive(Default)]
struct Counts {
    branches: u64,
    leaves: u64,
}
fn propagate(s: Partial, n: usize, out: &mut Vec<P>, counts: &mut Counts) {
    let Some(u) = (0..n).find(|&i| s.image[0][i] == NIL) else {
        counts.leaves += 1;
        let mut p = id();
        p[..n].copy_from_slice(&s.image[0][..n]);
        out.push(p);
        return;
    };
    for v in 0..n {
        if v == u || s.inverse[0][v] != NIL {
            continue;
        }
        counts.branches += 1;
        let mut next = s;
        if next.force(0, u, v, n) {
            propagate(next, n, out, counts);
        }
    }
}
fn completions(n: usize) -> (Vec<P>, Counts) {
    let mut out = vec![];
    let mut counts = Counts::default();
    propagate(Partial::empty(), n, &mut out, &mut counts);
    out.sort_unstable();
    assert!(out.windows(2).all(|w| w[0] != w[1]));
    (out, counts)
}
// Deliberately independent, unpruned conjugacy-class enumeration for timing
// and completeness controls. The least unused point anchors its next 3-cycle.
fn baseline_rec(c: &mut P, used: u16, n: usize, out: &mut Vec<P>, leaves: &mut u64) {
    let Some(u) = (0..n).find(|&i| used & (1 << i) == 0) else {
        *leaves += 1;
        let ci = inv(*c);
        let mut d = id();
        for i in 0..n {
            d[i] = ci[(i + n - 1) % n];
        }
        if order_three(d, n) {
            out.push(*c);
        }
        return;
    };
    let used = used | 1 << u;
    for v in u + 1..n {
        if used & (1 << v) != 0 {
            continue;
        }
        for w in u + 1..n {
            if w == v || used & (1 << w) != 0 {
                continue;
            }
            c[u] = v as u8;
            c[v] = w as u8;
            c[w] = u as u8;
            baseline_rec(c, used | 1 << v | 1 << w, n, out, leaves);
        }
    }
}
fn baseline(n: usize) -> (Vec<P>, u64) {
    let mut out = vec![];
    let mut leaves = 0;
    baseline_rec(&mut id(), 0, n, &mut out, &mut leaves);
    out.sort_unstable();
    (out, leaves)
}
fn profile(p: P) -> Vec<usize> {
    let mut used = 0u16;
    let mut out = vec![];
    for root in 0..N {
        if used & (1 << root) != 0 {
            continue;
        }
        let mut x = root;
        let mut len = 0;
        while used & (1 << x) == 0 {
            used |= 1 << x;
            len += 1;
            x = p[x] as usize;
        }
        out.push(len);
    }
    out.sort_unstable();
    out
}
// The smaller uniform class suffices for BOTH mixed completion orientations.
// This baseline costs 44,844,800 leaves once, rather than enumerating the much
// larger mixed class. A c,d -> d,(R^-1 c R) bijection swaps the two types.
fn target_counts(lengths: &[usize]) -> [u8; N + 1] {
    let mut counts = [0; N + 1];
    for &r in lengths {
        counts[r] += 1;
    }
    counts
}
fn mixed_rec(c: &mut P, used: u16, out: &mut Vec<P>, target: &[u8; N + 1]) {
    let Some(u) = (0..N).find(|&i| used & (1 << i) == 0) else {
        let r = std::array::from_fn(|i| ((i + N - 1) % N) as u8);
        let d = mul(inv(*c), r);
        let mut seen = 0u16;
        let mut counts = [0u8; N + 1];
        for root in 0..N {
            if seen & (1 << root) != 0 {
                continue;
            }
            let mut x = root;
            let mut len = 0;
            while seen & (1 << x) == 0 {
                seen |= 1 << x;
                len += 1;
                x = d[x] as usize;
            }
            counts[len] += 1;
            if counts[len] > target[len] {
                return;
            }
        }
        if &counts == target {
            out.push(*c);
        }
        return;
    };
    let used = used | 1 << u;
    for v in u + 1..N {
        if used & (1 << v) != 0 {
            continue;
        }
        for w in u + 1..N {
            if w == v || used & (1 << w) != 0 {
                continue;
            }
            c[u] = v as u8;
            c[v] = w as u8;
            c[w] = u as u8;
            mixed_rec(c, used | 1 << v | 1 << w, out, target);
        }
    }
}
fn relabel(p: P, labels: P) -> P {
    let mut q = id();
    for i in 0..N {
        q[labels[i] as usize] = labels[p[i] as usize];
    }
    q
}
fn product_labels(p: P) -> P {
    assert!(full(p));
    let mut labels = id();
    let mut x = 0;
    for k in 0..N {
        labels[x] = k as u8;
        x = p[x] as usize;
    }
    labels
}
struct Frames {
    pairs: Vec<[P; 2]>,
    lookup: HashMap<(u64, u64), (usize, usize)>,
}
impl Frames {
    fn add(&mut self, a: P, b: P) {
        let labels = product_labels(mul(a, b));
        let pair = [relabel(a, labels), relabel(b, labels)];
        if self.lookup.contains_key(&(pack(pair[0]), pack(pair[1]))) {
            return;
        }
        let f = self.pairs.len();
        for shift in 0..N {
            let rot = std::array::from_fn(|i| ((i + shift) % N) as u8);
            let key = (pack(relabel(pair[0], rot)), pack(relabel(pair[1], rot)));
            assert!(
                self.lookup.insert(key, (f, shift)).is_none(),
                "nontrivial pair stabilizer/duplicate frame"
            );
        }
        self.pairs.push(pair);
    }
    fn new(mixed: bool) -> Self {
        let profiles = if mixed {
            vec![vec![3, 3, 3, 3, 3], vec![4, 3, 3, 3, 2]]
        } else {
            vec![vec![3, 3, 3, 3, 3]]
        };
        let out = Self::profiles(profiles);
        assert_eq!(out.pairs.len(), if mixed { 42 } else { 2 });
        out
    }
    fn profiles(profiles: Vec<Vec<usize>>) -> Self {
        let mut out = Self {
            pairs: vec![],
            lookup: HashMap::new(),
        };
        for lengths in profiles {
            let mut b = id();
            let mut cycles = vec![];
            let mut start = 0;
            for len in lengths {
                let c = (start..start + len).collect::<Vec<_>>();
                cycle(&mut b, &c);
                cycles.push(c);
                start += len;
            }
            for s in 0..5 {
                for t in 0..5 {
                    if t == s {
                        continue;
                    }
                    let rest = (0..5).filter(|&i| i != s && i != t).collect::<Vec<_>>();
                    for &u in &rest {
                        for &v in &rest {
                            if u == v {
                                continue;
                            }
                            let w = *rest.iter().find(|&&i| i != u && i != v).unwrap();
                            for offset in 1..cycles[s].len() {
                                let mut a = id();
                                cycle(
                                    &mut a,
                                    &[cycles[s][0], cycles[u][0], cycles[v][0], cycles[w][0]],
                                );
                                cycle(&mut a, &[cycles[s][offset], cycles[t][0]]);
                                out.add(a, b);
                            }
                        }
                    }
                }
            }
        }
        out
    }
    fn key(&self, t: T) -> Option<(usize, u64)> {
        let p = mul(t[0], t[1]);
        if !full(p) {
            return None;
        }
        let l = product_labels(p);
        let pair = (pack(relabel(t[0], l)), pack(relabel(t[1], l)));
        let &(f, shift) = self.lookup.get(&pair).expect("complete tree-frame lookup");
        let c = relabel(t[2], l);
        let c = std::array::from_fn(|i| ((c[(i + shift) % N] as usize + N - shift) % N) as u8);
        Some((f, pack(c)))
    }
    fn tuple(&self, f: usize, c: P) -> T {
        let [a, b] = self.pairs[f];
        [a, b, c, inv(mul(mul(a, b), c))]
    }
}
fn preserves_cyclic_blocks(c: P, count: usize) -> bool {
    (0..count).all(|r| {
        (r..N)
            .step_by(count)
            .all(|i| c[i] as usize % count == c[r] as usize % count)
    })
}
fn images(p: P) -> String {
    p.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}
fn word_string(w: &[u8]) -> String {
    w.iter().map(u8::to_string).collect::<Vec<_>>().join(",")
}

fn short_negative(t: T, depth: usize, last: Option<u8>, word: &mut Vec<u8>) -> bool {
    if depth == 0 {
        return false;
    }
    for m in [1, 4, 0, 3, 2, 5] {
        if last.is_some_and(|p| m == (p + 3) % 6) {
            continue;
        }
        let mut u = t;
        braid(&mut u, m);
        word.push(m);
        // Passport indices are ALWAYS a permutation of (4,10,10,10).
        // Both complementary products have equal index; smooth iff full15.
        if !full(mul(u[0], u[1])) || short_negative(u, depth - 1, Some(m), word) {
            return true;
        }
        word.pop();
    }
    false
}
fn selected_profile(code: &str) -> Result<(Vec<usize>, usize), String> {
    Ok(match code {
        "A" => (vec![3, 3, 3, 3, 3], 2),
        "B" => (vec![4, 3, 3, 3, 2], 40),
        "C" => (vec![4, 4, 3, 2, 2], 60),
        "D" => (vec![5, 3, 3, 2, 2], 60),
        "E" => (vec![5, 4, 2, 2, 2], 40),
        "F" => (vec![6, 3, 2, 2, 2], 40),
        _ => return Err("slice must be A, B, C, D, E or F".into()),
    })
}
fn read_library(
    path: &Path,
    first: &[usize],
    last: &[usize],
    expected: usize,
) -> Result<Vec<P>, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut lines = text.lines();
    if lines.next() != Some("completion\tc") {
        return Err("library header".into());
    }
    let mut out = vec![];
    let r = std::array::from_fn(|i| ((i + N - 1) % N) as u8);
    let mut last = last.to_vec();
    last.sort_unstable();
    let mut first = first.to_vec();
    first.sort_unstable();
    for line in lines {
        let (rank, raw) = line.split_once('\t').ok_or("library columns")?;
        if rank.parse::<usize>().map_err(|e| e.to_string())? != out.len() {
            return Err("library rank".into());
        }
        let values = raw
            .split(',')
            .map(|s| s.parse::<u8>().map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        let c: P = values.try_into().map_err(|_| "library length")?;
        let mut sorted = c;
        sorted.sort_unstable();
        if sorted != id() {
            return Err("library images".into());
        }
        if profile(c) != first || profile(mul(inv(c), r)) != last {
            return Err("library passport/product".into());
        }
        out.push(c);
    }
    out.sort_unstable();
    if out.len() != expected || !out.windows(2).all(|w| w[0] != w[1]) {
        return Err("library coverage/duplicates".into());
    }
    Ok(out)
}
fn proof_words() -> Vec<Vec<u8>> {
    let mut out = vec![vec![1], vec![4], vec![0, 1], vec![0, 4], vec![3, 1]];
    fn visit(prefix: &mut Vec<u8>, depth: usize, out: &mut Vec<Vec<u8>>) {
        if depth == 0 {
            if [1, 4].contains(prefix.last().unwrap()) && !out.contains(prefix) {
                out.push(prefix.clone());
            }
            return;
        }
        for m in [1, 4, 0, 3, 2, 5] {
            if prefix.last().is_some_and(|&p| m == (p + 3) % 6) {
                continue;
            }
            prefix.push(m);
            visit(prefix, depth - 1, out);
            prefix.pop();
        }
    }
    for depth in 1..=3 {
        visit(&mut vec![], depth, &mut out);
    }
    out
}
fn full_composed(a: &P, c: &P) -> bool {
    let mut x = 0;
    for _ in 0..N - 1 {
        x = a[c[x] as usize] as usize;
        if x == 0 {
            return false;
        }
    }
    a[c[x] as usize] == 0
}
fn prefix_factors(a: P, b: P) -> [P; 5] {
    let p = mul(a, b);
    [
        mul(mul(inv(b), a), b),
        a,
        b,
        mul(mul(a, b), inv(a)),
        mul(mul(inv(p), b), p),
    ]
}
fn search_slice(frames: &Frames, cs: &[P], dictionary: &[Vec<u8>], threads: usize) -> Vec<u16> {
    let mut result = vec![0u16; frames.pairs.len() * cs.len()];
    let tokens = dictionary
        .iter()
        .enumerate()
        .map(|(i, w)| (w.clone(), i as u16 + 1))
        .collect::<HashMap<_, _>>();
    let chunk_frames = frames.pairs.len().div_ceil(threads).max(1);
    std::thread::scope(|scope| {
        for (part, chunk) in result.chunks_mut(chunk_frames * cs.len()).enumerate() {
            let tokens = &tokens;
            scope.spawn(move || {
                for (local, out) in chunk.chunks_mut(cs.len()).enumerate() {
                    let frame = part * chunk_frames + local;
                    let [a, b] = frames.pairs[frame];
                    let factors = prefix_factors(a, b);
                    for (&c, token) in cs.iter().zip(out) {
                        if let Some(j) = factors.iter().position(|f| !full_composed(f, &c)) {
                            *token = j as u16 + 1;
                            continue;
                        }
                        let t = frames.tuple(frame, c);
                        let mut word = vec![];
                        if (1..=3).any(|depth| short_negative(t, depth, None, &mut word)) {
                            *token = *tokens
                                .get(&word)
                                .expect("complete reduced short-word dictionary");
                        }
                    }
                }
            });
        }
    });
    result
}
fn slice(
    path: &Path,
    code: &str,
    first_code: &str,
    remaining: &str,
    library: Option<&Path>,
    threads: usize,
) -> Result<(), String> {
    let (lengths, expected_frames) = selected_profile(code)?;
    let (first_lengths, _) = selected_profile(first_code)?;
    let (remaining_lengths, _) = selected_profile(remaining)?;
    let expected_completions = match (first_code, remaining) {
        ("A", "A") => 2800,
        ("A", "B") => 63000,
        ("A", "C") => 106050,
        ("A", "D") => 111720,
        ("B", "B") => 1428000,
        ("B", "C") => 2425500,
        _ => return Err("supported factors: AA, AB, AC, AD, BB, BC".into()),
    };
    fs::create_dir(path).map_err(|e| e.to_string())?;
    let now = Instant::now();
    let (cs, counts) = if let Some(file) = library {
        (
            read_library(
                file,
                &first_lengths,
                &remaining_lengths,
                expected_completions,
            )?,
            Counts::default(),
        )
    } else if first_code == "A" && remaining == "A" {
        completions(15)
    } else {
        let (cs, stats) = factors::enumerate_parallel(&first_lengths, &remaining_lengths, threads);
        (
            cs,
            Counts {
                branches: stats.branches,
                leaves: stats.solutions,
            },
        )
    };
    assert_eq!(cs.len(), expected_completions);
    let frames = Frames::profiles(vec![lengths.clone()]);
    assert_eq!(frames.pairs.len(), expected_frames);
    let mut f = BufWriter::new(File::create(path.join("FRAMES.tsv")).map_err(|e| e.to_string())?);
    writeln!(f, "frame\ta\tb").map_err(|e| e.to_string())?;
    for (i, pair) in frames.pairs.iter().enumerate() {
        writeln!(f, "{i}\t{}\t{}", images(pair[0]), images(pair[1])).map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    let mut f =
        BufWriter::new(File::create(path.join("COMPLETIONS.tsv")).map_err(|e| e.to_string())?);
    writeln!(f, "completion\tc").map_err(|e| e.to_string())?;
    for (i, &c) in cs.iter().enumerate() {
        writeln!(f, "{i}\t{}", images(c)).map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    let construction_seconds = now.elapsed().as_secs_f64();
    let now = Instant::now();
    let dictionary = proof_words();
    let mut unresolved = 0;
    let mut max_word = 0;
    let result = search_slice(&frames, &cs, &dictionary, threads);
    let mut f =
        BufWriter::new(File::create(path.join("WITNESSES.u16le")).map_err(|e| e.to_string())?);
    let mut bytes = vec![0u8; 131072];
    for block in result.chunks(bytes.len() / 2) {
        for (i, &token) in block.iter().enumerate() {
            if token == 0 {
                unresolved += 1;
            } else {
                max_word = max_word.max(dictionary[token as usize - 1].len());
            }
            let pair = token.to_le_bytes();
            bytes[2 * i] = pair[0];
            bytes[2 * i + 1] = pair[1];
        }
        f.write_all(&bytes[..2 * block.len()])
            .map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    let search_seconds = now.elapsed().as_secs_f64();
    let mut f = BufWriter::new(File::create(path.join("WORDS.tsv")).map_err(|e| e.to_string())?);
    writeln!(f, "token\tword").map_err(|e| e.to_string())?;
    for (i, w) in dictionary.iter().enumerate() {
        writeln!(f, "{}\t{}", i + 1, word_string(w)).map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    let library_reused = library.is_some();
    let summary=format!("{{\"schema\":2,\"population\":\"degree15 special4x2, first ordinary {code}, remaining {first_code},{remaining}\",\"first_ordinary_profile\":{lengths:?},\"remaining_first_profile\":{first_lengths:?},\"remaining_last_profile\":{remaining_lengths:?},\"quotient\":\"simultaneous S15 conjugation\",\"degree\":15,\"genus\":3,\"frames\":{},\"completions\":{},\"candidate_states\":{},\"nonsmooth\":{},\"unresolved\":{unresolved},\"maximum_negative_word_length\":{max_word},\"word_dictionary_size\":{},\"constructor_branches_or_class_leaves\":{},\"validated_library_reused\":{library_reused},\"workers\":{threads},\"construction_seconds\":{construction_seconds},\"search_seconds\":{search_seconds}}}\n",frames.pairs.len(),cs.len(),frames.pairs.len()*cs.len(),frames.pairs.len()*cs.len()-unresolved,dictionary.len(),counts.branches);
    fs::write(path.join("SUMMARY.json"), &summary).map_err(|e| e.to_string())?;
    print!("{summary}");
    Ok(())
}

pub fn main(args: &[String]) -> Result<(), String> {
    if args.len() < 2 || args[1] == "--help" {
        println!("genus_three_finite OUT_DIR [--mixed | --slice A..F [--first-factor A|B] [--remaining A..D] [--library COMPLETIONS.tsv] [--threads 1..10]] [--baseline] [--repeat N]\nExact degree15 passport; fresh output directory required. No live-search mutation.");
        return Ok(());
    }
    let path = Path::new(&args[1]);
    let mut run_baseline = false;
    let mut repeats = 3;
    let mut mixed = false;
    let mut slice_code = None;
    let mut remaining = "A";
    let mut first_factor = "A";
    let mut threads = 1;
    let mut library = None;
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--baseline" => run_baseline = true,
            "--mixed" => mixed = true,
            "--slice" => {
                i += 1;
                slice_code = Some(args.get(i).ok_or("missing slice code")?.as_str());
            }
            "--remaining" => {
                i += 1;
                remaining = args.get(i).ok_or("missing remaining code")?.as_str();
            }
            "--first-factor" => {
                i += 1;
                first_factor = args.get(i).ok_or("missing first factor")?;
            }
            "--threads" => {
                i += 1;
                threads = args
                    .get(i)
                    .ok_or("missing thread count")?
                    .parse::<usize>()
                    .map_err(|e| e.to_string())?;
            }
            "--library" => {
                i += 1;
                library = Some(Path::new(args.get(i).ok_or("missing library path")?));
            }
            "--repeat" => {
                i += 1;
                repeats = args
                    .get(i)
                    .ok_or("missing repeat count")?
                    .parse::<usize>()
                    .map_err(|e| e.to_string())?;
            }
            _ => return Err(format!("unknown option {}", args[i])),
        }
        i += 1;
    }
    if repeats == 0 {
        return Err("repeat must be positive".into());
    }
    if let Some(code) = slice_code {
        if mixed || run_baseline {
            return Err("slice is separate from --mixed/--baseline".into());
        }
        if !["A", "B", "C", "D"].contains(&remaining) {
            return Err("remaining must be A, B, C or D".into());
        }
        if !(1..=10).contains(&threads) {
            return Err("threads1..10".into());
        }
        return slice(path, code, first_factor, remaining, library, threads);
    }
    if remaining != "A" {
        return Err("--remaining needs --slice".into());
    }
    if library.is_some() {
        return Err("--library needs --slice".into());
    }
    if first_factor != "A" || threads != 1 {
        return Err("first-factor/threads need --slice".into());
    }
    fs::create_dir(path).map_err(|e| e.to_string())?;
    let mut times = vec![];
    let mut base_times = vec![];
    let mut cset = vec![];
    let mut branches = 0;
    let mut base_leaves = 0;
    // Alternate order to reduce temperature/order bias under the live CPU load.
    for r in 0..repeats {
        if run_baseline && r % 2 == 0 {
            let now = Instant::now();
            let (b, n) = baseline(N);
            base_times.push(now.elapsed().as_secs_f64());
            cset = b;
            base_leaves = n;
        }
        let now = Instant::now();
        let (c, k) = completions(N);
        times.push(now.elapsed().as_secs_f64());
        branches = k.branches;
        if !cset.is_empty() {
            assert_eq!(c, cset, "different completion sets");
        }
        cset = c;
        if run_baseline && r % 2 == 1 {
            let now = Instant::now();
            let (b, n) = baseline(N);
            base_times.push(now.elapsed().as_secs_f64());
            assert_eq!(b, cset);
            base_leaves = n;
        }
    }
    assert_eq!(cset.len(), 2800, "Frobenius completeness count");
    let completion_seconds_start = Instant::now();
    let mut mixed_ab = vec![];
    let mut mixed_ba = vec![];
    if mixed {
        mixed_rec(
            &mut id(),
            0,
            &mut mixed_ab,
            &target_counts(&[4, 3, 3, 3, 2]),
        );
        mixed_ab.sort_unstable();
        let r = std::array::from_fn(|i| ((i + N - 1) % N) as u8);
        mixed_ba = mixed_ab.iter().map(|&c| mul(inv(c), r)).collect();
        mixed_ba.sort_unstable();
        assert_eq!(mixed_ab.len(), 63000);
        assert!(mixed_ab.windows(2).all(|w| w[0] != w[1]));
        assert!(mixed_ba.windows(2).all(|w| w[0] != w[1]));
    }
    let mixed_completion_seconds = completion_seconds_start.elapsed().as_secs_f64();
    let now = Instant::now();
    let frames = Frames::new(mixed);
    let mut keys = vec![];
    let mut tuples = vec![];
    let mut index = HashMap::new();
    let mut imprimitive = 0;
    for f in 0..frames.pairs.len() {
        let completions: Vec<&Vec<P>> = if !mixed || f >= 2 {
            vec![&cset]
        } else {
            vec![&mixed_ab, &mixed_ba]
        };
        for cs in completions {
            for &c in cs {
                let t = frames.tuple(f, c);
                let key = (f, pack(c));
                assert_eq!(frames.key(t), Some(key));
                assert!(index.insert(key, keys.len()).is_none());
                keys.push(key);
                tuples.push(t);
                let [a, b] = frames.pairs[f];
                if [3, 5]
                    .iter()
                    .any(|&k| [a, b, c].iter().all(|&p| preserves_cyclic_blocks(p, k)))
                {
                    imprimitive += 1;
                }
            }
        }
    }
    assert_eq!(keys.len(), if mixed { 364000 } else { 5600 });
    if !mixed {
        assert_eq!(imprimitive, 400);
    }
    let mut edges = vec![[None; 6]; keys.len()];
    let mut reverse = vec![vec![]; keys.len()];
    let mut failure: Vec<Option<(usize, usize)>> = vec![None; keys.len()];
    let mut queue = VecDeque::new();
    let mut immediate = 0;
    for (u, &t) in tuples.iter().enumerate() {
        for (m, word) in WORDS.iter().enumerate() {
            let mut t = t;
            apply(&mut t, word);
            if let Some(key) = frames.key(t) {
                let v = *index
                    .get(&key)
                    .expect("candidate set must be closed except nonsmooth exits");
                edges[u][m] = Some(v);
                reverse[v].push((u, m));
            } else if failure[u].is_none() {
                failure[u] = Some((u, m));
                queue.push_back(u);
                immediate += 1;
            }
        }
    }
    // Every internal mixed move has its literal inverse. Check before propagation.
    for u in 0..keys.len() {
        for m in 0..6 {
            if let Some(v) = edges[u][m] {
                assert_eq!(edges[v][(m + 3) % 6], Some(u));
            }
        }
    }
    while let Some(v) = queue.pop_front() {
        for &(u, m) in &reverse[v] {
            if failure[u].is_none() {
                failure[u] = Some((v, m));
                queue.push_back(u);
            }
        }
    }
    let graph_seconds = now.elapsed().as_secs_f64();
    let mut f = BufWriter::new(File::create(path.join("STATES.tsv")).map_err(|e| e.to_string())?);
    writeln!(f, "state\tframe\ta\tb\tc\td\tstatus\tword\tedges").map_err(|e| e.to_string())?;
    let mut max_word = 0;
    let mut safe = 0;
    for u in 0..keys.len() {
        let mut word = vec![];
        let mut at = u;
        while let Some((next, m)) = failure[at] {
            word.extend_from_slice(WORDS[m]);
            if at == next {
                break;
            }
            at = next;
            assert!(word.len() < 2 * keys.len() + 3, "failure predecessor cycle");
        }
        let status = if failure[u].is_some() {
            let mut t = tuples[u];
            apply(&mut t, &word);
            assert!(!full(mul(t[0], t[1])));
            "nonsmooth"
        } else {
            safe += 1;
            "closed_mixed_requires_full_audit"
        };
        max_word = max_word.max(word.len());
        let t = tuples[u];
        let es = edges[u]
            .iter()
            .map(|e| e.map_or("-".into(), |v| v.to_string()))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            f,
            "{u}\t{}\t{}\t{}\t{}\t{}\t{status}\t{}\t{es}",
            keys[u].0,
            images(t[0]),
            images(t[1]),
            images(t[2]),
            images(t[3]),
            word_string(&word)
        )
        .map_err(|e| e.to_string())?;
    }
    f.flush().map_err(|e| e.to_string())?;
    let population = if mixed {
        "degree15 mixed (4 2 1^9),(4 2 3^3),(3^5)^2"
    } else {
        "degree15 pure (4 2 1^9),(3^5)^3"
    };
    let result=format!("{{\"schema\":1,\"population\":\"{population}\",\"quotient\":\"simultaneous S15 conjugation\",\"degree\":15,\"genus\":3,\"frames\":{},\"uniform_completions\":{},\"mixed_completions_per_orientation\":{},\"candidate_states\":{},\"imprimitive\":{imprimitive},\"primitive_states\":{},\"immediate_escape_states\":{immediate},\"nonsmooth\":{},\"safe_mixed_states\":{safe},\"unresolved\":0,\"maximum_negative_word_length\":{max_word},\"propagation_branches\":{branches},\"baseline_class_leaves\":{base_leaves},\"propagation_seconds\":{times:?},\"baseline_seconds\":{base_times:?},\"mixed_completion_seconds\":{mixed_completion_seconds},\"graph_seconds\":{graph_seconds}}}\n",frames.pairs.len(),cset.len(),mixed_ab.len(),keys.len(),keys.len()-imprimitive,keys.len()-safe);
    fs::write(path.join("SUMMARY.json"), &result).map_err(|e| e.to_string())?;
    print!("{result}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn propagation_matches_unpruned_classes() {
        for n in [3, 6, 9, 12] {
            let (a, _) = completions(n);
            let (b, _) = baseline(n);
            assert_eq!(a, b, "degree {n}");
        }
    }
    #[test]
    fn exact_fifteen_count() {
        let (c, k) = completions(15);
        assert_eq!(c.len(), 2800);
        assert_eq!(k.leaves, 2800);
        let r = std::array::from_fn(|i| ((i + 14) % 15) as u8);
        for c in c {
            let d = mul(inv(c), r);
            assert!(order_three(c, 15) && order_three(d, 15));
            assert_eq!(mul(c, d), r);
        }
    }
    #[test]
    fn specialized_key_survives_every_product_rotation_and_relabeling() {
        let frames = Frames::new(false);
        let (cs, _) = completions(15);
        for f in 0..2 {
            for &c in cs.iter().step_by(137) {
                let t = frames.tuple(f, c);
                let key = frames.key(t);
                for root in 0..15 {
                    let l = std::array::from_fn(|i| ((7 * i + root) % 15) as u8);
                    assert_eq!(frames.key(t.map(|p| relabel(p, l))), key);
                }
                for m in 0..6 {
                    let mut q = t;
                    braid(&mut q, m);
                    braid(&mut q, (m + 3) % 6);
                    assert_eq!(q, t);
                }
            }
        }
    }
    #[test]
    fn general_tree_frames_and_mixed_completion_bijection() {
        let f = Frames::new(true);
        assert_eq!(f.lookup.len(), 42 * 15);
        for pair in &f.pairs {
            assert!(full(mul(pair[0], pair[1])));
            assert_eq!(profile(pair[0]), [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4]);
        }
    }
    #[test]
    fn path_solver_matches_complete_unpruned_mixed_class() {
        let mut literal = vec![];
        mixed_rec(&mut id(), 0, &mut literal, &target_counts(&[4, 3, 3, 3, 2]));
        literal.sort_unstable();
        let (propagated, _) = factors::enumerate(&[3; 5], &[4, 3, 3, 3, 2]);
        assert_eq!(literal.len(), 63000);
        assert_eq!(propagated, literal);
    }
    #[test]
    #[ignore = "bounded paired performance experiment; not a routine test"]
    fn benchmark_path_solver_against_literal_class() {
        let mut reference_times = vec![];
        let mut path_times = vec![];
        for repetition in 0..5 {
            let literal = || {
                let mut out = vec![];
                mixed_rec(&mut id(), 0, &mut out, &target_counts(&[4, 3, 3, 3, 2]));
                out.sort_unstable();
                out
            };
            let propagated = || factors::enumerate(&[3; 5], &[4, 3, 3, 3, 2]).0;
            let (a, b) = if repetition % 2 == 0 {
                let now = Instant::now();
                let a = literal();
                reference_times.push(now.elapsed().as_secs_f64());
                let now = Instant::now();
                let b = propagated();
                path_times.push(now.elapsed().as_secs_f64());
                (a, b)
            } else {
                let now = Instant::now();
                let b = propagated();
                path_times.push(now.elapsed().as_secs_f64());
                let now = Instant::now();
                let a = literal();
                reference_times.push(now.elapsed().as_secs_f64());
                (a, b)
            };
            assert_eq!(a.len(), 63000);
            assert_eq!(a, b);
        }
        println!("FACTOR_BENCHMARK {{\"same_complete_set\":63000,\"literal_class_leaves\":44844800,\"literal_seconds\":{reference_times:?},\"path_seconds\":{path_times:?}}}");
    }
    #[test]
    fn short_word_oracle_keeps_smooth_control_and_exposes_negative() {
        // Published positive genus3 degree15 single-5 elliptic control. It has
        // the same indices, not the type4x2 passport being ruled out here.
        let mut smooth = [id(); 4];
        cycle(&mut smooth[0], &[0, 2, 3, 1, 4]);
        for (j, cycles) in [
            [[0, 5, 10], [1, 6, 11], [2, 7, 12], [3, 8, 13], [4, 9, 14]],
            [[0, 8, 11], [1, 7, 14], [2, 9, 13], [3, 6, 10], [4, 5, 12]],
            [[0, 7, 13], [1, 9, 12], [2, 6, 14], [3, 5, 11], [4, 8, 10]],
        ]
        .iter()
        .enumerate()
        {
            for c in cycles {
                cycle(&mut smooth[j + 1], c);
            }
        }
        assert_eq!(
            mul(mul(smooth[0], smooth[1]), mul(smooth[2], smooth[3])),
            id()
        );
        for depth in 1..=3 {
            assert!(!short_negative(smooth, depth, None, &mut vec![]));
        }
        let f = Frames::new(false);
        let (cs, _) = completions(15);
        let mut word = vec![];
        assert!(short_negative(f.tuple(0, cs[0]), 3, None, &mut word));
        let mut t = f.tuple(0, cs[0]);
        apply(&mut t, &word);
        assert!(!full(mul(t[0], t[1])));
    }
}
