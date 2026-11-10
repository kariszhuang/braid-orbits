//! Exact fixed-full-cycle connection coefficients by two-sided path propagation.
//! Each partial permutation is a union of closed cycles and directed paths.
//! Joining paths either consumes a requested cycle or forces its last edge.
//! The second permutation is not guessed: c(u)=v iff d(v+1)=u.
const MAX: usize = 15;
const UNSET: u8 = 255;
type P = [u8; MAX];
#[derive(Clone, Copy)]
#[cfg(test)]
struct State {
    image: [[u8; MAX]; 2],
    inverse: [[u8; MAX]; 2],
    remaining: [[u8; MAX + 1]; 2],
}
#[derive(Default, Debug)]
pub struct Stats {
    pub branches: u64,
    pub solutions: u64,
}
#[cfg(test)]
impl State {
    fn assign(&mut self, side: usize, u: usize, v: usize, n: usize) -> bool {
        if self.image[side][u] != UNSET {
            return self.image[side][u] as usize == v;
        }
        if self.inverse[side][v] != UNSET {
            return false;
        }
        // Before insertion a chain v -> ... -> u closes exactly one cycle.
        let mut x = v;
        let mut len = 1;
        while x != u && self.image[side][x] != UNSET {
            x = self.image[side][x] as usize;
            len += 1;
        }
        if x == u {
            if self.remaining[side][len] == 0 {
                return false;
            }
            self.remaining[side][len] -= 1;
        }
        self.image[side][u] = v as u8;
        self.inverse[side][v] = u as u8;
        let paired = if side == 0 {
            self.assign(1, (v + 1) % n, u, n)
        } else {
            self.assign(0, v, (u + n - 1) % n, n)
        };
        if !paired {
            return false;
        }
        // Re-evaluate path capacities after cycle consumption as well as after
        // joins. A path of maximal permitted size has a forced closing edge.
        for k in 0..2 {
            let max = (1..=n)
                .rev()
                .find(|&r| self.remaining[k][r] != 0)
                .unwrap_or(0);
            let mut paths = [0u8; MAX + 1];
            let mut forced = None;
            for start in 0..n {
                if self.inverse[k][start] != UNSET {
                    continue;
                }
                let mut end = start;
                let mut size = 1;
                while self.image[k][end] != UNSET {
                    end = self.image[k][end] as usize;
                    size += 1;
                }
                if size > max {
                    return false;
                }
                paths[size] += 1;
                if size == max {
                    forced = Some((start, end));
                    break;
                }
            }
            if let Some((start, end)) = forced {
                if !self.assign(k, end, start, n) {
                    return false;
                }
            } else if max <= 4 {
                // Exact necessary capacity inequalities for remaining 2/3/4
                // cycles. Prefer placing length3 paths in 3-cycles, retaining
                // 4-cycles for two length2 paths. Singletons fill unused slots.
                let r = self.remaining[k];
                if paths[4] > r[4] {
                    return false;
                }
                let fours = r[4] - paths[4];
                if paths[3] > r[3] + fours {
                    return false;
                }
                let threes = paths[3].min(r[3]);
                let long_in_four = paths[3] - threes;
                if paths[2] > r[2] + r[3] - threes + 2 * (fours - long_in_four) {
                    return false;
                }
            }
        }
        true
    }
}
#[cfg(test)]
fn visit(s: State, n: usize, out: &mut Vec<P>, stats: &mut Stats) {
    let Some(u) = (0..n).find(|&i| s.image[0][i] == UNSET) else {
        assert!(s.remaining.iter().all(|r| r.iter().all(|&c| c == 0)));
        let mut p = std::array::from_fn(|i| i as u8);
        p[..n].copy_from_slice(&s.image[0][..n]);
        out.push(p);
        stats.solutions += 1;
        return;
    };
    for v in 0..n {
        if s.inverse[0][v] != UNSET {
            continue;
        }
        stats.branches += 1;
        let mut next = s;
        if next.assign(0, u, v, n) {
            visit(next, n, out, stats);
        }
    }
}
#[cfg(test)]
pub fn enumerate_reference(mu: &[usize], nu: &[usize]) -> (Vec<P>, Stats) {
    let n = mu.iter().sum::<usize>();
    assert!(n > 0 && n <= MAX && n == nu.iter().sum::<usize>());
    let mut s = State {
        image: [[UNSET; MAX]; 2],
        inverse: [[UNSET; MAX]; 2],
        remaining: [[0; MAX + 1]; 2],
    };
    for (k, partition) in [mu, nu].iter().enumerate() {
        for &r in partition.iter() {
            assert!(r > 0 && r <= n);
            s.remaining[k][r] += 1;
        }
    }
    let mut out = vec![];
    let mut stats = Stats::default();
    visit(s, n, &mut out, &mut stats);
    out.sort_unstable();
    assert!(out.windows(2).all(|w| w[0] != w[1]));
    (out, stats)
}

// Only path heads and tails need endpoint information: a missing outgoing edge
// is at a tail, and a missing incoming edge is at a head. Internal entries may
// remain stale. This replaces repeated walks over all15 vertices per edge.
#[derive(Clone, Copy)]
struct Cached {
    image: [[u8; MAX]; 2],
    inverse: [[u8; MAX]; 2],
    head_at_tail: [[u8; MAX]; 2],
    tail_at_head: [[u8; MAX]; 2],
    size_at_head: [[u8; MAX]; 2],
    heads: [[u16; MAX + 1]; 2],
    path_lengths: [u16; 2],
    open_heads: [u16; 2],
    available_lengths: [u16; 2],
    remaining: [[u8; MAX + 1]; 2],
}
impl Cached {
    fn remove(&mut self, k: usize, len: usize, head: usize) {
        debug_assert!(self.heads[k][len] & (1 << head) != 0);
        self.heads[k][len] &= !(1 << head);
        if self.heads[k][len] == 0 {
            self.path_lengths[k] &= !(1 << len);
        }
    }
    fn assign(&mut self, k: usize, u: usize, v: usize, n: usize) -> bool {
        if self.image[k][u] != UNSET {
            return self.image[k][u] as usize == v;
        }
        if self.inverse[k][v] != UNSET {
            return false;
        }
        let head = self.head_at_tail[k][u] as usize;
        let left = self.size_at_head[k][head] as usize;
        if head == v {
            if self.remaining[k][left] == 0 {
                return false;
            }
            self.remaining[k][left] -= 1;
            if self.remaining[k][left] == 0 {
                self.available_lengths[k] &= !(1 << left);
            }
            self.remove(k, left, head);
            self.open_heads[k] &= !(1 << head);
        } else {
            let right = self.size_at_head[k][v] as usize;
            let size = left + right;
            if size > MAX || (self.available_lengths[k] as u32 >> size) == 0 {
                return false;
            }
            let tail = self.tail_at_head[k][v] as usize;
            self.remove(k, left, head);
            self.remove(k, right, v);
            self.heads[k][size] |= 1 << head;
            self.path_lengths[k] |= 1 << size;
            self.open_heads[k] &= !(1 << v);
            self.head_at_tail[k][tail] = head as u8;
            self.tail_at_head[k][head] = tail as u8;
            self.size_at_head[k][head] = size as u8;
        }
        self.image[k][u] = v as u8;
        self.inverse[k][v] = u as u8;
        let paired = if k == 0 {
            self.assign(1, (v + 1) % n, u, n)
        } else {
            self.assign(0, v, (u + n - 1) % n, n)
        };
        if !paired {
            return false;
        }
        for side in 0..2 {
            let max = if self.available_lengths[side] == 0 {
                0
            } else {
                15 - self.available_lengths[side].leading_zeros() as usize
            };
            if (self.path_lengths[side] as u32 >> (max + 1)) != 0 {
                return false;
            }
            let forced = self.heads[side][max];
            if forced != 0 {
                let start = forced.trailing_zeros() as usize;
                let end = self.tail_at_head[side][start] as usize;
                if !self.assign(side, end, start, n) {
                    return false;
                }
            } else if max <= 4 {
                let p4 = self.heads[side][4].count_ones() as u8;
                let p3 = self.heads[side][3].count_ones() as u8;
                let p2 = self.heads[side][2].count_ones() as u8;
                let r = self.remaining[side];
                if p4 > r[4] {
                    return false;
                }
                let fours = r[4] - p4;
                if p3 > r[3] + fours {
                    return false;
                }
                let threes = p3.min(r[3]);
                if p2 > r[2] + r[3] - threes + 2 * (fours - (p3 - threes)) {
                    return false;
                }
            }
        }
        true
    }
    fn tail(&self) -> Option<usize> {
        let mask = self.path_lengths[0];
        if mask == 0 {
            return None;
        }
        let size = 15 - mask.leading_zeros() as usize;
        let head = self.heads[0][size].trailing_zeros() as usize;
        Some(self.tail_at_head[0][head] as usize)
    }
}
fn cached_visit(s: Cached, n: usize, out: &mut Vec<P>, stats: &mut Stats) {
    let Some(u) = s.tail() else {
        assert_eq!(s.available_lengths, [0, 0]);
        let mut p = std::array::from_fn(|i| i as u8);
        p[..n].copy_from_slice(&s.image[0][..n]);
        out.push(p);
        stats.solutions += 1;
        return;
    };
    let mut heads = s.open_heads[0];
    while heads != 0 {
        let v = heads.trailing_zeros() as usize;
        heads &= heads - 1;
        stats.branches += 1;
        let mut next = s;
        if next.assign(0, u, v, n) {
            cached_visit(next, n, out, stats);
        }
    }
}
#[allow(dead_code)] // single-core diagnostic surface used by independent set tests
pub fn enumerate(mu: &[usize], nu: &[usize]) -> (Vec<P>, Stats) {
    enumerate_parallel(mu, nu, 1)
}
pub fn enumerate_parallel(mu: &[usize], nu: &[usize], threads: usize) -> (Vec<P>, Stats) {
    let n = mu.iter().sum::<usize>();
    assert!(n > 0 && n <= MAX && n == nu.iter().sum::<usize>() && threads > 0);
    let id: P = std::array::from_fn(|i| i as u8);
    let mut s = Cached {
        image: [[UNSET; MAX]; 2],
        inverse: [[UNSET; MAX]; 2],
        head_at_tail: [id; 2],
        tail_at_head: [id; 2],
        size_at_head: [[1; MAX]; 2],
        heads: [[0; MAX + 1]; 2],
        path_lengths: [2; 2],
        open_heads: [(1 << n) - 1; 2],
        available_lengths: [0; 2],
        remaining: [[0; MAX + 1]; 2],
    };
    for k in 0..2 {
        s.heads[k][1] = (1 << n) - 1;
    }
    for (k, partition) in [mu, nu].iter().enumerate() {
        for &r in partition.iter() {
            assert!(r > 0 && r <= n);
            s.remaining[k][r] += 1;
            s.available_lengths[k] |= 1 << r;
        }
    }
    let u = s.tail().unwrap();
    let mut tasks = vec![];
    for v in 0..n {
        let mut next = s;
        if next.assign(0, u, v, n) {
            tasks.push(next);
        }
    }
    let cursor = std::sync::atomic::AtomicUsize::new(0);
    let mut out = vec![];
    let mut stats = Stats {
        branches: n as u64,
        solutions: 0,
    };
    std::thread::scope(|scope| {
        let mut handles = vec![];
        for _ in 0..threads.min(tasks.len().max(1)) {
            let tasks = &tasks;
            let cursor = &cursor;
            handles.push(scope.spawn(move || {
                let mut local = vec![];
                let mut counts = Stats::default();
                loop {
                    let at = cursor.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if at >= tasks.len() {
                        break;
                    }
                    cached_visit(tasks[at], n, &mut local, &mut counts);
                }
                (local, counts)
            }));
        }
        for handle in handles {
            let (mut local, counts) = handle.join().expect("factor worker panic");
            out.append(&mut local);
            stats.branches += counts.branches;
            stats.solutions += counts.solutions;
        }
    });
    out.sort_unstable();
    assert!(out.windows(2).all(|w| w[0] != w[1]));
    (out, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn profile(p: P, n: usize) -> Vec<usize> {
        let mut seen = [false; MAX];
        let mut ls = vec![];
        for root in 0..n {
            if seen[root] {
                continue;
            }
            let mut x = root;
            let mut len = 0;
            while !seen[x] {
                seen[x] = true;
                len += 1;
                x = p[x] as usize;
            }
            ls.push(len);
        }
        ls.sort_unstable();
        ls
    }
    fn literal(mu: &[usize], nu: &[usize]) -> Vec<P> {
        let n = mu.iter().sum::<usize>();
        let mut p: P = std::array::from_fn(|i| i as u8);
        let mut out = vec![];
        loop {
            if profile(p, n) == mu {
                let mut inverse = p;
                for i in 0..n {
                    inverse[p[i] as usize] = i as u8;
                }
                let mut d = p;
                for i in 0..n {
                    d[i] = inverse[(i + n - 1) % n];
                }
                if profile(d, n) == nu {
                    out.push(p);
                }
            }
            let Some(i) = (0..n - 1).rev().find(|&i| p[i] < p[i + 1]) else {
                break;
            };
            let j = (i + 1..n).rev().find(|&j| p[j] > p[i]).unwrap();
            p.swap(i, j);
            p[i + 1..n].reverse();
        }
        out
    }
    #[test]
    fn generic_matches_literal_small_connections() {
        for (a, b) in [
            (vec![1, 2, 3], vec![2, 2, 2]),
            (vec![2, 2, 3], vec![2, 2, 3]),
            (vec![2, 2, 4], vec![2, 3, 3]),
            (vec![2, 3, 4], vec![3, 3, 3]),
        ] {
            assert_eq!(enumerate(&a, &b).0, literal(&a, &b));
        }
    }
    #[test]
    fn generic_matches_uniform_character_count() {
        assert_eq!(enumerate(&[3; 5], &[3; 5]).0.len(), 2800);
    }
    #[test]
    fn cached_parallel_and_scanning_paths_are_identical() {
        for (a, b) in [
            (vec![1, 2, 3], vec![2, 2, 2]),
            (vec![2, 2, 3], vec![2, 2, 3]),
            (vec![3; 5], vec![3; 5]),
            (vec![3; 5], vec![2, 3, 3, 3, 4]),
        ] {
            let reference = enumerate_reference(&a, &b).0;
            assert_eq!(enumerate(&a, &b).0, reference);
            assert_eq!(enumerate_parallel(&a, &b, 4).0, reference);
        }
    }
    #[test]
    #[ignore = "bounded alternating-order exact-set performance experiment"]
    fn benchmark_cached_endpoint_solver() {
        let mut reference = vec![];
        let mut one = vec![];
        let mut four = vec![];
        for repetition in 0..5 {
            let mut answers = vec![];
            for variant in if repetition % 2 == 0 {
                [0, 1, 2]
            } else {
                [2, 1, 0]
            } {
                let now = std::time::Instant::now();
                let (value, counts) = match variant {
                    0 => enumerate_reference(&[3; 5], &[2, 3, 3, 3, 4]),
                    1 => enumerate_parallel(&[3; 5], &[2, 3, 3, 3, 4], 1),
                    _ => enumerate_parallel(&[3; 5], &[2, 3, 3, 3, 4], 4),
                };
                let seconds = now.elapsed().as_secs_f64();
                match variant {
                    0 => reference.push(seconds),
                    1 => one.push(seconds),
                    _ => four.push(seconds),
                }
                assert_eq!(counts.solutions, 63000);
                answers.push(value);
            }
            assert!(answers.windows(2).all(|w| w[0] == w[1]));
        }
        println!("CACHED_BENCHMARK {{\"reference_one_core\":{reference:?},\"cached_one_core\":{one:?},\"cached_four_cores\":{four:?},\"same_complete_set\":63000}}");
    }
}
