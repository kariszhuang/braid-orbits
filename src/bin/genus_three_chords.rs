//! Exact full-product slices with a transposition or double transposition.
//! The third/fourth entries each have three cycles. Standalone:
//! no discovery imports, geometric filters, caps treated as negatives, or GPU.
use std::{
    collections::VecDeque,
    fs,
    io::{BufWriter, Write},
    path::Path,
    time::Instant,
};
const MAX: usize = 11;
type P = [u8; MAX];
type T = [P; 4];
fn id(n: usize) -> P {
    let mut p = [0; MAX];
    for i in 0..n {
        p[i] = i as u8;
    }
    p
}
fn mul(a: &P, b: &P, n: usize) -> P {
    let mut p = [0; MAX];
    for i in 0..n {
        p[i] = a[b[i] as usize];
    }
    p
}
fn inv(a: &P, n: usize) -> P {
    let mut p = [0; MAX];
    for i in 0..n {
        p[a[i] as usize] = i as u8;
    }
    p
}
fn labels(p: &P, n: usize) -> (usize, P) {
    let mut seen = 0u16;
    let mut lab = [0; MAX];
    let mut count = 0;
    for r in 0..n {
        if seen & (1 << r) != 0 {
            continue;
        }
        let mut x = r;
        while seen & (1 << x) == 0 {
            seen |= 1 << x;
            lab[x] = count;
            x = p[x] as usize;
        }
        count += 1;
    }
    (count as usize, lab)
}
fn smooth(t: &T, n: usize) -> bool {
    let k = n - labels(&mul(&t[0], &t[1], n), n).0;
    let ix = t.map(|p| n - labels(&p, n).0);
    k == ix[0] + ix[1] || k == ix[2] + ix[3]
}
fn braid(t: &T, m: u8, n: usize) -> T {
    let mut u = *t;
    let i = m as usize % 3;
    let (a, b) = (t[i], t[i + 1]);
    if m < 3 {
        u[i] = mul(&mul(&a, &b, n), &inv(&a, n), n);
        u[i + 1] = a;
    } else {
        u[i] = b;
        u[i + 1] = mul(&mul(&inv(&b, n), &a, n), &b, n);
    }
    u
}
fn replay(t: &T, w: &[u8], n: usize) -> T {
    w.iter().fold(*t, |u, &m| braid(&u, m, n))
}
fn notation(p: &P, n: usize) -> String {
    let mut seen = 0u16;
    let mut s = String::new();
    for r in 0..n {
        if seen & (1 << r) != 0 || p[r] as usize == r {
            continue;
        }
        let mut x = r;
        let mut cycle = vec![];
        while seen & (1 << x) == 0 {
            seen |= 1 << x;
            cycle.push((x + 1).to_string());
            x = p[x] as usize;
        }
        s += &format!("({})", cycle.join(" "));
    }
    if s.is_empty() {
        s = "()".into();
    }
    s
}
fn seed(t: &T, n: usize) -> String {
    format!(
        "[{}]",
        t.iter()
            .map(|p| notation(p, n))
            .collect::<Vec<_>>()
            .join(", ")
    )
}
fn twist(k: usize, n: usize) -> Vec<u8> {
    if k <= n - k {
        vec![0; 2 * k]
    } else {
        vec![3; 2 * (n - k)]
    }
}
fn first_failure(t: &T, n: usize, bank: &[Vec<u8>]) -> Option<Vec<u8>> {
    for w in bank {
        if !smooth(&replay(t, w, n), n) {
            return Some(w.clone());
        }
    }
    let mut queue = VecDeque::from([(*t, Vec::<u8>::new())]);
    while let Some((u, w)) = queue.pop_front() {
        if w.len() == 4 {
            continue;
        }
        for m in 0..6u8 {
            if w.last().is_some_and(|&last| (last + 3) % 6 == m) {
                continue;
            }
            let v = braid(&u, m, n);
            let mut path = w.clone();
            path.push(m);
            if !smooth(&v, n) {
                return Some(path);
            }
            queue.push_back((v, path));
        }
    }
    None // A bounded local test never certifies an all-smooth component.
}
struct Run {
    n: usize,
    frame: usize,
    edges: Vec<(usize, usize)>,
    t: P,
    e: P,
    p: P,
    partitions: u64,
    constructed: u64,
    qualified: u64,
    negative: u64,
    pending: u64,
    bank: Vec<Vec<u8>>,
    leads: BufWriter<fs::File>,
    proofs: BufWriter<fs::File>,
    survivors: BufWriter<fs::File>,
}
impl Run {
    fn compatible(&self, lab: &P, assigned: usize, k: usize) -> bool {
        let mut mask = 0u8;
        let mut ready = true;
        for &(a, b) in &self.edges {
            let (a, b) = ((a + k) % self.n, (b + k) % self.n);
            if a < assigned && b < assigned && lab[a] == lab[b] {
                return false;
            }
            ready &= a < assigned && b < assigned;
            if a < assigned {
                mask |= 1 << lab[a];
            }
            if b < assigned {
                mask |= 1 << lab[b];
            }
        }
        self.edges.len() == 1 || !ready || mask.count_ones() == 3
    }
    fn leaf(&mut self, b: &P) {
        self.constructed += 1;
        let c = inv(&mul(&self.p, b, self.n), self.n);
        let (cycles, lab) = labels(&c, self.n);
        if cycles != 3 {
            return;
        }
        self.qualified += 1;
        let t = [self.t, self.e, *b, c];
        assert!(smooth(&t, self.n));
        let bad = (0..self.n).find(|&k| !self.compatible(&lab, self.n, k));
        let word = if let Some(k) = bad {
            let mut w = twist(k, self.n);
            w.extend([5, 4]);
            Some(w)
        } else {
            first_failure(&t, self.n, &self.bank)
        };
        let text = seed(&t, self.n);
        if let Some(w) = word {
            assert!(!smooth(&replay(&t, &w, self.n), self.n));
            self.negative += 1;
            if self.bank.len() < 64 && !self.bank.contains(&w) {
                self.bank.push(w.clone());
            }
            writeln!(self.leads, "{}\t{}", self.negative, text).unwrap();
            writeln!(
                self.proofs,
                "{}\tnonsmooth\t{}\t\t{}",
                self.negative,
                w.iter()
                    .map(|m| m.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                text
            )
            .unwrap();
        } else {
            self.pending += 1;
            writeln!(
                self.survivors,
                "{}\t{}\t{}",
                self.frame, self.qualified, text
            )
            .unwrap();
        }
    }
    fn cycle_orders(&mut self, blocks: &[Vec<usize>], j: usize, b: &mut P) {
        if j == 3 {
            self.leaf(b);
            return;
        }
        fn permute(
            run: &mut Run,
            blocks: &[Vec<usize>],
            j: usize,
            b: &mut P,
            order: &mut Vec<usize>,
            k: usize,
        ) {
            if k == order.len() {
                for i in 0..order.len() {
                    b[order[i]] = order[(i + 1) % order.len()] as u8;
                }
                run.cycle_orders(blocks, j + 1, b);
                return;
            }
            for l in k..order.len() {
                order.swap(k, l);
                permute(run, blocks, j, b, order, k + 1);
                order.swap(k, l);
            }
        }
        let mut order = blocks[j].clone();
        permute(self, blocks, j, b, &mut order, 1);
    }
    fn coloring(&mut self, x: usize, used: u8, color: &mut P) {
        if x == self.n {
            if used != 3 {
                return;
            }
            self.partitions += 1;
            let blocks = (0..3)
                .map(|k| (0..self.n).filter(|&i| color[i] == k).collect())
                .collect::<Vec<Vec<usize>>>();
            let mut b = id(self.n);
            self.cycle_orders(&blocks, 0, &mut b);
            return;
        }
        for k in 0..=used.min(2) {
            // Restricted-growth colors represent UNLABELLED partitions once.
            color[x] = k;
            if (0..self.n).all(|shift| self.compatible(color, x + 1, shift)) {
                self.coloring(x + 1, used.max(k + 1), color);
            }
        }
    }
}
fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    assert!([3, 4].contains(&args.len()));
    let mode = args.get(3).map(String::as_str).unwrap_or("transposition");
    assert!(["transposition", "double"].contains(&mode));
    let n = args[1].parse::<usize>().unwrap();
    assert!([7, 9, 11].contains(&n));
    let out = Path::new(&args[2]);
    fs::create_dir(out).unwrap();
    let start = Instant::now();
    let p = {
        let mut p = id(n);
        for i in 0..n {
            p[i] = ((i + 1) % n) as u8;
        }
        p
    };
    let mut run = Run {
        n,
        frame: 0,
        edges: vec![],
        t: id(n),
        e: id(n),
        p,
        partitions: 0,
        constructed: 0,
        qualified: 0,
        negative: 0,
        pending: 0,
        bank: vec![],
        leads: BufWriter::new(fs::File::create(out.join("LEADS.tsv")).unwrap()),
        proofs: BufWriter::new(fs::File::create(out.join("RESOLUTIONS.tsv")).unwrap()),
        survivors: BufWriter::new(fs::File::create(out.join("SURVIVORS.tsv")).unwrap()),
    };
    writeln!(run.proofs, "candidate\tstatus\tword\tfolder\tseed").unwrap();
    writeln!(run.survivors, "frame\tqualifying\tseed").unwrap();
    let mut frames = vec![];
    if mode == "transposition" {
        for r in 1..=n / 2 {
            let mut t = id(n);
            t.swap(0, r);
            frames.push(t);
        }
    } else {
        for i in 0..n {
            for j in i + 1..n {
                for k in j + 1..n {
                    for l in k + 1..n {
                        for edges in [[(i, j), (k, l)], [(i, l), (j, k)]] {
                            let mut t = id(n);
                            for (a, b) in edges {
                                t.swap(a, b);
                            }
                            assert_eq!(labels(&mul(&t, &p, n), n).0, 3);
                            let minimum = (0..n)
                                .map(|shift| {
                                    let mut u = id(n);
                                    for x in 0..n {
                                        u[(x + shift) % n] = ((t[x] as usize + shift) % n) as u8;
                                    }
                                    u
                                })
                                .min()
                                .unwrap();
                            if t == minimum {
                                frames.push(t);
                            }
                        }
                    }
                }
            }
        }
        frames.sort();
        assert_eq!(frames.len(), 2 * (n * (n - 1) * (n - 2) * (n - 3) / 24) / n);
    }
    let mut rows = vec![];
    for (frame, t) in frames.into_iter().enumerate() {
        run.frame = frame;
        run.t = t;
        run.edges = (0..n)
            .filter_map(|i| {
                let j = t[i] as usize;
                (i < j).then_some((i, j))
            })
            .collect();
        run.e = mul(&run.t, &run.p, n);
        let old = (
            run.partitions,
            run.constructed,
            run.qualified,
            run.negative,
            run.pending,
        );
        run.coloring(0, 0, &mut [0; MAX]);
        rows.push(format!("{{\"frame\":{frame},\"special\":[{}],\"partitions\":{},\"constructed\":{},\"qualified\":{},\"negative\":{},\"pending\":{}}}",t[..n].iter().map(u8::to_string).collect::<Vec<_>>().join(","),run.partitions-old.0,run.constructed-old.1,run.qualified-old.2,run.negative-old.3,run.pending-old.4));
    }
    run.leads.flush().unwrap();
    run.proofs.flush().unwrap();
    run.survivors.flush().unwrap();
    // Guided journal counts ONLY persisted literal negatives. Pending local
    // survivors have their own durable table and remain undecided in COVERAGE.
    fs::write(out.join("COUNTS.tsv"),format!("degree\tlane\tseed\tconstructed\tscope\ttransitive\tpassport\tshort_reject\tduplicate\tclosed\tunresolved\tseconds\n{n}\tsolved\t0\t{}\t{}\t{}\t{}\t0\t0\t0\t0\t{}\n",run.negative,run.negative,run.negative,run.negative,start.elapsed().as_secs_f64())).unwrap();
    let rank = if mode == "transposition" { 1 } else { 2 };
    fs::write(out.join("COVERAGE.json"),format!("{{\"degree\":{n},\"genus\":{},\"mode\":\"{mode}\",\"ordered_indices\":[{rank},{},{},{}],\"frames\":[{}],\"constructed\":{},\"qualified\":{},\"negative\":{},\"pending\":{},\"seconds\":{}}}\n",(n-5)/2,n-rank-1,n-3,n-3,rows.join(","),run.constructed,run.qualified,run.negative,run.pending,start.elapsed().as_secs_f64())).unwrap();
}
