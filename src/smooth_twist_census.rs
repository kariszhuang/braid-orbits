// Twist-constrained exact STAR census. No S_n enumeration.
//
// Setting: (a,b,c,e) with abce = 1, a the unique smallest index s, STAR index
// lists (every pair avoiding a has index sum > n-1). Group convention
// (xy)(i) = x(y(i)); move m in 0..2 is sigma_{m+1}, m+3 its inverse.
//
// Frames. Every orbit has a state with a first; its cusp (12)|(34) forces
// (a,b) additive, i.e. a <=_T p = ab. Such a are noncrossing partitions on each
// p-cycle, so frames are built from noncrossing factorizations of one
// representative p per cycle type, then canonicalized up to S_n.
//
// Lemma T1. sigma3^{2k} conjugates (c,e) by ce = p^{-1}; its cusp (23)|(41)
// forces (p^{-k} e p^k, a) additive. Word: 2^(2k) then the rotation 0,1,2,
// then cusp (12)|(34). Additivity depends only on the two cycle partitions
// (incidence forest), so the cycle partition of e is forest-compatible with
// every p-rotation of the partition of a. We enumerate exactly those set
// partitions (union-find pruning is exact: a cycle in the incidence graph is
// created by the edge that closes it) and then every cyclic order.
//
// Lemma T2. sigma2^{2j} conjugates (b,c) by bc = r^{-1}, r = ea; cusp (12)|(34)
// forces (a, r^{-j} b r^j) additive. Word 1^(2j). Then sigma2^{-1} gives
// (a, r^{-j} c r^j) additive. Word 1^(2j),4.
//
// Every excluded start is therefore excluded by an explicit literal word;
// T1 exclusions are never enumerated, and --brute re-enumerates S_n to check
// both the enumeration and the T1 words literally.
// AtomicUsize comes from smooth_star_census.rs

const TW_COLUMNS: &str = "e_partitions\te_orders\tindex_eligible\ttransitive\tt2_reject\tshort_reject\tscan_reject\tclosed_occurrences\tunresolved";

const TW_N: usize = 32;

// sc_pair with 5-bit labels (n <= 25): canonical form of a pair up to S_n.
fn tw_pair(a: Perm, b: Perm, n: usize) -> ([u128; 2], u128) {
    assert!(n <= 25);
    let mut remaining = (1u64 << n) - 1;
    let mut components = Vec::new();
    while remaining != 0 {
        let first = remaining.trailing_zeros() as usize;
        let mut mask = 0u64;
        let mut best = (usize::MAX, [u128::MAX; 2]);
        let mut aut = 0u128;
        let mut roots = remaining;
        let mut initial = true;
        while roots != 0 {
            let root = if initial { first } else { roots.trailing_zeros() as usize };
            let mut labels = [u8::MAX; MAX_D];
            let mut order = [0usize; MAX_D];
            labels[root] = 0;
            order[0] = root;
            let (mut head, mut len) = (0, 1);
            while head < len {
                let x = order[head];
                head += 1;
                for p in [a, b] {
                    let y = p.img[x] as usize;
                    if labels[y] == u8::MAX { labels[y] = len as u8; order[len] = y; len += 1; }
                }
            }
            if initial {
                for &x in &order[..len] { mask |= 1 << x; }
                roots = mask;
                initial = false;
            }
            roots &= !(1 << root);
            let mut key = [0u128; 2];
            for (j, p) in [a, b].iter().enumerate() {
                for i in 0..len { key[j] |= (labels[p.img[order[i]] as usize] as u128) << (5 * i); }
            }
            let value = (len, key);
            if value < best { best = value; aut = 1; } else if value == best { aut += 1; }
        }
        remaining &= !mask;
        components.push((best.0, best.1, aut));
    }
    components.sort_unstable();
    let mut key = [0u128; 2];
    let (mut offset, mut centralizer, mut multiplicity) = (0, 1u128, 0u128);
    let mut previous = None;
    for (len, local, aut) in components {
        let signature = (len, local);
        if previous == Some(signature) { multiplicity += 1; } else { multiplicity = 1; previous = Some(signature); }
        centralizer *= aut * multiplicity;
        for j in 0..2 {
            for i in 0..len { key[j] |= (((local[j] >> (5 * i)) & 31) + offset as u128) << (5 * (offset + i)); }
        }
        offset += len;
    }
    assert_eq!(offset, n);
    (key, centralizer)
}

fn tw_unpack(k: [u128; 2], n: usize) -> (Perm, Perm) {
    let mut t = [Perm::identity(n); 2];
    for j in 0..2 { for i in 0..n { t[j].img[i] = ((k[j] >> (5 * i)) & 31) as u8; } }
    (t[0], t[1])
}

fn tw_order(p: Perm, n: usize) -> usize {
    let mut seen = [false; MAX_D];
    let mut l = 1usize;
    for i in 0..n {
        if seen[i] { continue; }
        let mut j = i;
        let mut len = 0;
        while !seen[j] { seen[j] = true; j = p.img[j] as usize; len += 1; }
        let g = { let (mut x, mut y) = (l, len); while y != 0 { let t = x % y; x = y; y = t; } x };
        l = l / g * len;
    }
    l
}

fn tw_labels(p: Perm, n: usize) -> [u8; TW_N] {
    let mut lab = [u8::MAX; TW_N];
    for i in 0..n {
        if lab[i] != u8::MAX { continue; }
        let mut j = i;
        while lab[j] == u8::MAX { lab[j] = i as u8; j = p.img[j] as usize; }
    }
    lab
}

fn tw_additive(x: Perm, y: Perm, n: usize) -> bool {
    sm_index(x.compose(y, n), n) == sm_index(x, n) + sm_index(y, n)
}

// Noncrossing partitions of each p-cycle, in p-order; callback gets a.
fn tw_nc(
    i: usize, n: usize, starts: &[bool; TW_N], stack: &mut Vec<usize>, block_of: &mut [usize; TW_N],
    nblocks: usize, emit: &mut dyn FnMut(&[usize; TW_N], usize),
) {
    if i == n { emit(block_of, nblocks); return; }
    let saved: Vec<usize> = if starts[i] { std::mem::take(stack) } else { Vec::new() };
    stack.push(nblocks);
    block_of[i] = nblocks;
    tw_nc(i + 1, n, starts, stack, block_of, nblocks + 1, emit);
    stack.pop();
    for t in 0..stack.len() {
        let above = stack[t + 1..].to_vec();
        stack.truncate(t + 1);
        block_of[i] = stack[t];
        tw_nc(i + 1, n, starts, stack, block_of, nblocks, emit);
        stack.extend(above);
    }
    if starts[i] { *stack = saved; }
}

fn tw_frames(n: usize, patterns: &[[usize; 4]], all_partners: bool) -> Vec<([u128; 2], u128)> {
    let acc = |ia: usize, ib: usize| patterns.iter().any(|q| q[0] == ia && if all_partners { q[1..].contains(&ib) } else { q[3] == ib });
    nc_frames(n, &acc)
}

// Canonical index-additive pairs (a,b), a <=_T ab, with accept(ind a, ind b).
fn nc_frames(n: usize, accept: &dyn Fn(usize, usize) -> bool) -> Vec<([u128; 2], u128)> {
    assert!(n <= 25);
    let mut parts = Vec::new();
    sc_partitions(n, 1, &mut Vec::new(), &mut parts);
    let mut frames = std::collections::BTreeMap::new();
    for lengths in parts {
        let ip = n - lengths.len();
        if !(0..=ip).any(|ia| accept(ia, ip - ia)) { continue; }
        let mut p = Perm::identity(n);
        let mut starts = [false; TW_N];
        let mut offset = 0;
        for &len in &lengths {
            starts[offset] = true;
            for i in 0..len { p.img[offset + i] = (offset + (i + 1) % len) as u8; }
            offset += len;
        }
        let mut raw = 0u128;
        let mut emit = |block_of: &[usize; TW_N], nb: usize| {
            raw += 1;
            let ia = n - nb;
            let ib = ip - ia;
            if !accept(ia, ib) { return; }
            let mut a = Perm::identity(n);
            let mut last = [usize::MAX; TW_N];
            let mut first = [usize::MAX; TW_N];
            for x in 0..n {
                let bl = block_of[x];
                if first[bl] == usize::MAX { first[bl] = x; } else { a.img[last[bl]] = x as u8; }
                last[bl] = x;
            }
            for bl in 0..nb { a.img[last[bl]] = first[bl] as u8; }
            let b = a.inverse(n).compose(p, n);
            assert_eq!(sm_index(a, n) + sm_index(b, n), ip, "noncrossing factor is geodesic");
            let (key, z) = tw_pair(a, b, n);
            if let Some(old) = frames.insert(key, z) { assert_eq!(old, z); }
        };
        tw_nc(0, n, &starts, &mut Vec::new(), &mut [0; TW_N], 0, &mut emit);
        // Catalan numbers by C(k) = C(k-1) 2(2k-1)/(k+1); (2l)! overflows u128 for l >= 18
        let catalan = |l: usize| -> u128 { let mut c = 1u128; for k in 1..=l as u128 { c = c * 2 * (2 * k - 1) / (k + 1); } c };
        let expected: u128 = lengths.iter().map(|&l| catalan(l)).product();
        assert_eq!(raw, expected, "Catalan mass of geodesic factorizations");
    }
    frames.into_iter().collect()
}

struct TwEnum<'a> {
    n: usize,
    m: usize,
    labs: &'a [[u8; TW_N]],
    block: [u8; TW_N],
    first: [u8; TW_N],
    nblocks: usize,
}

fn tw_find(uf: &[u8; TW_N], mut x: u8) -> u8 { while uf[x as usize] != x { x = uf[x as usize]; } x }

impl<'a> TwEnum<'a> {
    fn run(&mut self, i: usize, ufs: &mut Vec<[u8; TW_N]>, emit: &mut dyn FnMut(&[u8; TW_N], usize)) {
        let n = self.n;
        if i == n {
            if self.nblocks == self.m { emit(&self.block, self.m); }
            return;
        }
        let k = self.labs.len();
        // join an existing block
        if self.nblocks + (n - i - 1) >= self.m {
            for j in 0..self.nblocks {
                let f = self.first[j] as usize;
                let mut ok = true;
                for t in 0..k {
                    if tw_find(&ufs[t], self.labs[t][i]) == tw_find(&ufs[t], self.labs[t][f]) { ok = false; break; }
                }
                if !ok { continue; }
                let saved: Vec<[u8; TW_N]> = ufs.clone();
                for t in 0..k {
                    let (x, y) = (tw_find(&ufs[t], self.labs[t][i]), tw_find(&ufs[t], self.labs[t][f]));
                    ufs[t][x as usize] = y;
                }
                self.block[i] = j as u8;
                self.run(i + 1, ufs, emit);
                *ufs = saved;
            }
        }
        // open a new block
        if self.nblocks < self.m {
            self.block[i] = self.nblocks as u8;
            self.first[self.nblocks] = i as u8;
            self.nblocks += 1;
            self.run(i + 1, ufs, emit);
            self.nblocks -= 1;
        }
    }
}

// Every e with the given block partition: each block gets every cyclic order.
fn tw_orders(n: usize, block: &[u8; TW_N], m: usize, emit: &mut dyn FnMut(Perm)) {
    let mut members: Vec<Vec<u8>> = vec![Vec::new(); m];
    for x in 0..n { members[block[x] as usize].push(x as u8); }
    fn rec(bi: usize, members: &mut Vec<Vec<u8>>, e: &mut Perm, emit: &mut dyn FnMut(Perm)) {
        if bi == members.len() { emit(*e); return; }
        let len = members[bi].len();
        // permute positions 1..len, keeping members[bi][0] first (Heap-free: recursive swaps)
        fn perm(bi: usize, pos: usize, members: &mut Vec<Vec<u8>>, e: &mut Perm, emit: &mut dyn FnMut(Perm)) {
            let len = members[bi].len();
            if pos >= len.saturating_sub(1) {
                let cyc = &members[bi];
                for t in 0..len { e.img[cyc[t] as usize] = cyc[(t + 1) % len]; }
                rec(bi + 1, members, e, emit);
                return;
            }
            for j in pos..len {
                members[bi].swap(pos, j);
                perm(bi, pos + 1, members, e, emit);
                members[bi].swap(pos, j);
            }
        }
        if len <= 1 {
            if len == 1 { let x = members[bi][0]; e.img[x as usize] = x; }
            rec(bi + 1, members, e, emit);
        } else {
            perm(bi, 1, members, e, emit);
        }
    }
    let mut e = Perm::identity(n);
    rec(0, &mut members, &mut e, emit);
}


// Atkinson: the smallest block containing {0, y}; primitive iff it is always everything.
fn tw_primitive(gens: &[Perm], n: usize) -> bool {
    for y in 1..n {
        let mut uf: Vec<usize> = (0..n).collect();
        fn find(uf: &mut Vec<usize>, mut x: usize) -> usize { while uf[x] != x { uf[x] = uf[uf[x]]; x = uf[x]; } x }
        let (r0, ry) = (find(&mut uf, 0), find(&mut uf, y));
        uf[r0] = ry;
        let mut queue = vec![(0usize, y)];
        while let Some((u, v)) = queue.pop() {
            for g in gens {
                let (gu, gv) = (g.img[u] as usize, g.img[v] as usize);
                let (ru, rv) = (find(&mut uf, gu), find(&mut uf, gv));
                if ru != rv { uf[ru] = rv; queue.push((gu, gv)); }
            }
        }
        let root = find(&mut uf, 0);
        if (0..n).filter(|&x| find(&mut uf, x) == root).count() < n { return false; }
    }
    true
}


// Cyclic orders of the blocks of e, built point by point. Each placed point x
// fixes c(x) = p^{-1}(e^{-1}(x)); every c-arc is checked against every p-rotation
// of a (lemma T1c: sigma1^{2k+1} forces (p^k a p^-k, c) additive), exactly as in
// TwEnum. The c-cycle count is pruned to [tmin, tmax].
struct TwOrd<'a> {
    n: usize,
    pinv: Perm,
    labs: &'a [[u8; TW_N]],
    tmin: usize,
    tmax: usize,
}

impl<'a> TwOrd<'a> {
    // add arc x -> y of c; false if a cycle appears in some incidence graph
    fn arc(&self, x: u8, y: u8, ufs: &mut [[u8; TW_N]], cuf: &mut [u8; TW_N], closed: &mut usize, open: &mut usize) -> bool {
        let (cx, cy) = (tw_find(cuf, x), tw_find(cuf, y));
        if cx == cy {
            *closed += 1;
            *open -= 1;
        } else {
            for t in 0..self.labs.len() {
                if tw_find(&ufs[t], self.labs[t][x as usize]) == tw_find(&ufs[t], self.labs[t][y as usize]) { return false; }
            }
            for t in 0..self.labs.len() {
                let (u, v) = (tw_find(&ufs[t], self.labs[t][x as usize]), tw_find(&ufs[t], self.labs[t][y as usize]));
                ufs[t][u as usize] = v;
            }
            cuf[cx as usize] = cy;
            *open -= 1;
        }
        *closed <= self.tmax && *closed + *open >= self.tmin
    }

    #[allow(clippy::too_many_arguments)]
    fn place(&self, members: &[Vec<u8>], bi: usize, prev: u8, rest: &mut Vec<u8>, e: &mut Perm,
             ufs: &mut Vec<[u8; TW_N]>, cuf: &mut [u8; TW_N], closed: usize, open: usize, emit: &mut dyn FnMut(Perm)) {
        if bi == members.len() { emit(*e); return; }
        let x0 = members[bi][0];
        if rest.is_empty() {
            // close the cycle: e(prev) = x0, so c(x0) = p^{-1}(prev)
            e.img[prev as usize] = x0;
            let (mut cl, mut op) = (closed, open);
            let (su, sc) = (ufs.clone(), *cuf);
            if self.arc(x0, self.pinv.img[prev as usize], ufs, cuf, &mut cl, &mut op) {
                if bi + 1 == members.len() {
                    if op == 0 && cl >= self.tmin { emit(*e); }
                } else {
                    let mut next: Vec<u8> = members[bi + 1][1..].to_vec();
                    self.place(members, bi + 1, members[bi + 1][0], &mut next, e, ufs, cuf, cl, op, emit);
                }
            }
            *ufs = su; *cuf = sc;
            return;
        }
        for i in 0..rest.len() {
            let y = rest[i];
            e.img[prev as usize] = y;
            let (mut cl, mut op) = (closed, open);
            let (su, sc) = (ufs.clone(), *cuf);
            // e(prev) = y, so c(y) = p^{-1}(prev)
            if self.arc(y, self.pinv.img[prev as usize], ufs, cuf, &mut cl, &mut op) {
                rest.swap_remove(i);
                self.place(members, bi, y, rest, e, ufs, cuf, cl, op, emit);
                rest.push(y);
                let last = rest.len() - 1;
                rest.swap(i, last);
            }
            *ufs = su; *cuf = sc;
        }
    }
}

fn tw_orders_c(n: usize, block: &[u8; TW_N], m: usize, pinv: Perm, labs: &[[u8; TW_N]], tmin: usize, tmax: usize, emit: &mut dyn FnMut(Perm)) {
    let mut members: Vec<Vec<u8>> = vec![Vec::new(); m];
    for x in 0..n { members[block[x] as usize].push(x as u8); }
    let o = TwOrd { n, pinv, labs, tmin, tmax };
    let _ = o.n;
    let mut ufs: Vec<[u8; TW_N]> = labs.iter().map(|_| { let mut u = [0u8; TW_N]; for i in 0..TW_N { u[i] = i as u8; } u }).collect();
    let mut cuf = [0u8; TW_N];
    for i in 0..TW_N { cuf[i] = i as u8; }
    let mut e = Perm::identity(n);
    let mut rest: Vec<u8> = members[0][1..].to_vec();
    o.place(&members, 0, members[0][0], &mut rest, &mut e, &mut ufs, &mut cuf, 0, n, emit);
}

// T1c literal word: 0^(2k+1) then 0,1,2; cusp (12)|(34) then holds (p^k a p^-k, c) | (e, .).
fn tw_t1c_word(t: [Perm; 4], n: usize, p: Perm) -> Option<Vec<usize>> {
    let mut ak = t[0];
    for k in 0..tw_order(p, n) {
        if !tw_additive(ak, t[2], n) {
            let mut w = vec![0usize; 2 * k + 1];
            w.extend([0, 1, 2]);
            return Some(w);
        }
        ak = ak.conjugate_by(p, n);
    }
    None
}


// Block-first search over e. Points of U = supp(a) come first, so the e-blocks
// meeting U are completed (as sets, then with cyclic orders) before any other
// block. At that moment r = ea is known on W' = union of those blocks (outside
// W' it does not matter: r^j a r^-j is supported in W'), so lemma T2 for b,
// (r^j a r^-j, b) additive for every j, is checked once for the whole subtree.
// T1 (partition) and T1c (c-arcs) are pruned exactly as before.
#[derive(Clone)]
struct TwState {
    assigned: [bool; TW_N],
    blocks: Vec<Vec<u8>>,
    ufs: Vec<[u8; TW_N]>,
    e: Perm,
    cufs: Vec<[u8; TW_N]>,
    cuf: [u8; TW_N],
    closed: usize,
    open: usize,
    u_blocks: usize,
    u_done: bool,
}

struct TwSearch<'a> {
    n: usize,
    m: usize,
    labs: &'a [[u8; TW_N]],
    order: Vec<u8>,
    is_u: [bool; TW_N],
    pinv: Perm,
    a: Perm,
    acyc: Vec<Vec<u8>>,
    blab: [u8; TW_N],
    tmin: usize,
    tmax: usize,
}

impl<'a> TwSearch<'a> {
    fn fresh(&self) -> TwState {
        let mut id = [0u8; TW_N];
        for i in 0..TW_N { id[i] = i as u8; }
        TwState {
            assigned: [false; TW_N], blocks: Vec::new(), ufs: vec![id; self.labs.len()], e: Perm::identity(self.n),
            cufs: vec![id; self.labs.len()], cuf: id, closed: 0, open: self.n, u_blocks: 0, u_done: false,
        }
    }

    fn arc(&self, st: &mut TwState, x: u8, y: u8) -> bool {
        let (cx, cy) = (tw_find(&st.cuf, x), tw_find(&st.cuf, y));
        if cx == cy {
            st.closed += 1;
            st.open -= 1;
        } else {
            for t in 0..self.labs.len() {
                if tw_find(&st.cufs[t], self.labs[t][x as usize]) == tw_find(&st.cufs[t], self.labs[t][y as usize]) { return false; }
            }
            for t in 0..self.labs.len() {
                let (u, v) = (tw_find(&st.cufs[t], self.labs[t][x as usize]), tw_find(&st.cufs[t], self.labs[t][y as usize]));
                st.cufs[t][u as usize] = v;
            }
            st.cuf[cx as usize] = cy;
            st.open -= 1;
        }
        st.closed <= self.tmax && st.closed + st.open >= self.tmin
    }

    fn step(&self, st: &mut TwState, stats: &mut [u128; 2], emit: &mut dyn FnMut(Perm)) {
        let next_u = self.order.iter().copied().find(|&x| self.is_u[x as usize] && !st.assigned[x as usize]);
        if let Some(x) = next_u { self.open_block(st, x, stats, emit); return; }
        if !st.u_done {
            st.u_blocks = st.blocks.len();
            let mut s2 = st.clone();
            self.order_from(&mut s2, 0, stats, emit);
            return;
        }
        let next = self.order.iter().copied().find(|&x| !st.assigned[x as usize]);
        match next {
            Some(x) => self.open_block(st, x, stats, emit),
            None => {
                if st.blocks.len() != self.m { return; }
                stats[0] += 1;
                let mut s2 = st.clone();
                let bi = s2.u_blocks;
                self.order_from(&mut s2, bi, stats, emit);
            }
        }
    }

    fn open_block(&self, st: &mut TwState, x: u8, stats: &mut [u128; 2], emit: &mut dyn FnMut(Perm)) {
        if st.blocks.len() >= self.m { return; }
        let mut s2 = st.clone();
        s2.assigned[x as usize] = true;
        s2.blocks.push(vec![x]);
        let pos = self.order.iter().position(|&y| y == x).unwrap();
        self.extend(&mut s2, pos, stats, emit);
    }

    // add members of larger rank to the last block, or close it
    fn extend(&self, st: &mut TwState, last_pos: usize, stats: &mut [u128; 2], emit: &mut dyn FnMut(Perm)) {
        let unassigned = (0..self.n).filter(|&y| !st.assigned[y]).count();
        if st.blocks.len() + unassigned >= self.m {
            let mut s2 = st.clone();
            self.step(&mut s2, stats, emit);
        }
        let bl = st.blocks.len() - 1;
        let first = st.blocks[bl][0] as usize;
        for pos in last_pos + 1..self.order.len() {
            let y = self.order[pos] as usize;
            if st.assigned[y] { continue; }
            if (0..self.labs.len()).any(|t| tw_find(&st.ufs[t], self.labs[t][y]) == tw_find(&st.ufs[t], self.labs[t][first])) { continue; }
            if st.blocks.len() + (unassigned - 1) < self.m { continue; }
            let mut s2 = st.clone();
            for t in 0..self.labs.len() {
                let (u, v) = (tw_find(&s2.ufs[t], self.labs[t][y]), tw_find(&s2.ufs[t], self.labs[t][first]));
                s2.ufs[t][u as usize] = v;
            }
            s2.assigned[y] = true;
            s2.blocks[bl].push(y as u8);
            self.extend(&mut s2, pos, stats, emit);
        }
    }

    // cyclic orders of blocks bi.. (up to u_blocks in the U phase, else to the end)
    fn order_from(&self, st: &mut TwState, bi: usize, stats: &mut [u128; 2], emit: &mut dyn FnMut(Perm)) {
        let end = if st.u_done { st.blocks.len() } else { st.u_blocks };
        if bi == end {
            if !st.u_done {
                if !self.t2b_ok(st) { stats[1] += 1; return; }
                let mut s2 = st.clone();
                s2.u_done = true;
                self.step(&mut s2, stats, emit);
            } else if st.open == 0 && st.closed >= self.tmin && st.closed <= self.tmax {
                emit(st.e);
            }
            return;
        }
        let x0 = st.blocks[bi][0];
        let rest: Vec<u8> = st.blocks[bi][1..].to_vec();
        self.place(st, bi, x0, x0, rest, stats, emit);
    }

    #[allow(clippy::too_many_arguments)]
    fn place(&self, st: &mut TwState, bi: usize, x0: u8, prev: u8, rest: Vec<u8>, stats: &mut [u128; 2], emit: &mut dyn FnMut(Perm)) {
        if rest.is_empty() {
            let mut s2 = st.clone();
            s2.e.img[prev as usize] = x0;
            if self.arc(&mut s2, x0, self.pinv.img[prev as usize]) { self.order_from(&mut s2, bi + 1, stats, emit); }
            return;
        }
        for i in 0..rest.len() {
            let y = rest[i];
            let mut s2 = st.clone();
            s2.e.img[prev as usize] = y;
            if self.arc(&mut s2, y, self.pinv.img[prev as usize]) {
                let mut r2 = rest.clone();
                r2.remove(i);
                self.place(&mut s2, bi, x0, y, r2, stats, emit);
            }
        }
    }

    // r = ea on the union of the ordered U-blocks (r^j a r^-j is supported there).
    // (r^j a r^-j, b) additive <=> the rotated a-cycles form a hyperforest on the
    // b-cycles; track only the positions r^j(u), u in U, until they all return.
    fn t2b_ok(&self, st: &TwState) -> bool {
        let n = self.n;
        let mut r = [0u8; TW_N];
        for x in 0..n { r[x] = x as u8; }
        for bl in &st.blocks[..st.u_blocks] {
            for &x in bl { r[x as usize] = st.e.img[self.a.img[x as usize] as usize]; }
        }
        let mut cur = [0u8; TW_N];
        for x in 0..n { cur[x] = x as u8; }
        loop {
            let mut uf = [0u8; TW_N];
            for i in 0..TW_N { uf[i] = i as u8; }
            for cyc in &self.acyc {
                let mut roots = [0u8; TW_N];
                for (t, &x) in cyc.iter().enumerate() {
                    let root = tw_find(&uf, self.blab[cur[x as usize] as usize]);
                    if roots[..t].contains(&root) { return false; }
                    roots[t] = root;
                }
                for t in 1..cyc.len() { uf[roots[t] as usize] = roots[0]; }
            }
            let mut back = true;
            for cyc in &self.acyc {
                for &x in cyc {
                    let y = r[cur[x as usize] as usize];
                    cur[x as usize] = y;
                    if y != x { back = false; }
                }
            }
            if back { return true; }
        }
    }
}

// T2 (b-version) literal word for a full tuple: 1^(2j), cusp (12)|(34).
fn tw_t2b_word(t: [Perm; 4], n: usize) -> Option<Vec<usize>> {
    let r = t[3].compose(t[0], n);
    let rinv = r.inverse(n);
    let (ia, ib) = (sm_index(t[0], n), sm_index(t[1], n));
    let mut bj = t[1];
    for j in 0..tw_order(r, n) {
        if sm_index(t[0].compose(bj, n), n) != ia + ib { return Some(vec![1usize; 2 * j]); }
        bj = rinv.compose(bj, n).compose(r, n);
    }
    None
}


// Arc-first search over e (v7). Points of U = supp(a) first; each choice
// e(x) = y adds an e-arc x -> y (lemma T1, forest against every p-rotation of a)
// and the c-arc y -> p^{-1}(x) (lemma T1c), both with exact union-find tests and
// cycle-count bounds. When e (resp. e^{-1}) is known on a whole cycle A_i of a,
// r(A_i) = e(A_i) (resp. r^{-1}(A_i) = a^{-1} e^{-1}(A_i)) is added as a
// hyperedge on the b-cycles: lemma T2 at j = 1 and j = -1 (a hyperforest test).
// Union-find changes are undone on backtrack through a log.
struct ArcSearch<'a> {
    n: usize,
    m: usize,
    tmin: usize,
    tmax: usize,
    labs: &'a [[u8; TW_N]],
    pinv: Perm,
    a: Perm,
    ainv: Perm,
    order: Vec<u8>,
    acyc_of: [u8; TW_N],
    acyc: Vec<Vec<u8>>,
    blab: [u8; TW_N],
    // v10: letters for lemma words, the gated lemma list and the symmetry group
    lt: [Perm; 6], // a, a^-1, b, b^-1, p, p^-1
    lem: Vec<TwLemma>,
    sym: Vec<[u8; TW_N]>,
    lazy_depth: usize,
}

// ---------------------------------------------------------------------------
// v10 in-search lemmas. Every word w below ends at a state whose cusp (12)|(34)
// pairs a conjugate g a g^-1 with a conjugate h Y h^-1 (Y in {b, c, e}); the
// other half pairs conjugates of the remaining two entries. When their index
// sum is >= n that half is never additive (the lemma is GATED on this), so the
// word fails iff W(A_i) = h^-1 g (A_i) is not a hyperforest on the Y-cycles.
// W was derived symbolically in the free group and checked numerically
// (evidence/primitive_star_20261103/V10_LEMMAS.md); --brute re-checks every
// prune against the literal words. Partial information is sound: points of
// W(A_i) that are not yet defined are left out (a sub-hyperedge), and e- and
// c-cycles are read as current path components (merging only adds cycles).
// Families: W_j = X^j for every j != 0 up to --jmax in both directions.
// Letters: a A b B p P e E c C (capital = inverse), product read left to right,
// applied right to left.
#[derive(Clone)]
struct TwLemma {
    name: &'static str,
    x: Vec<u8>,
    xi: Vec<u8>,
    pre: Vec<u8>,
    y: u8,          // 0 = b, 1 = e, 2 = c
    other: [u8; 2], // entries of the other cusp half: 1 = b, 2 = c, 3 = e
    fam: u8,        // 0 single word, else the move prefix of the family (+1); 8 = T2b; 9 = T2c
    word: &'static [usize],
    // shift m: the lemma applied to (a, b, p^-m c p^m, p^-m e p^m), which is the state
    // after sigma1^(2m) relabelled by p^-m. sh = p^m, shi = p^-m (identity for m = 0).
    shift: i64,
    sh: [u8; TW_N],
    shi: [u8; TW_N],
}

const TW_LEMMA_SPECS: &[(&str, &str, &str, char, &str, u8, &[usize])] = &[
    // name, X, pre, Y, other, family prefix code, single word
    ("T2c", "ea", "", 'c', "be", 9, &[]),
    ("T2b", "ea", "", 'b', "ce", 8, &[]),  // the walk lemma; used only with shifts m != 0
    ("R3", "ca", "", 'b', "ce", 1, &[]),   // sigma1 sigma2^{2j}
    ("R4", "EB", "", 'b', "ce", 4, &[]),   // sigma1^-1 sigma2^{2j}
    ("R5", "cECB", "", 'b', "ce", 3, &[]), // sigma3 sigma2^{2j} (never pruned in tests; not default)
    ("S131", "", "bCB", 'e', "bc", 0, &[1, 3, 1]),
    ("S404", "", "CBc", 'e', "bc", 0, &[4, 0, 4]),
    ("S0044", "", "Acab", 'b', "ce", 0, &[0, 0, 4, 4]),
    ("S0404", "", "aBAC", 'c', "be", 0, &[0, 4, 0, 4]),
    ("S0444", "", "CAC", 'e', "bc", 0, &[0, 4, 4, 4]),
    ("S1131", "", "bcBCB", 'e', "bc", 0, &[1, 1, 3, 1]),
    ("S1311", "", "BabcbCB", 'c', "be", 0, &[1, 3, 1, 1]),
    ("S1331", "", "CBAbCB", 'b', "ce", 0, &[1, 3, 3, 1]),
    ("S3111", "", "BabcB", 'e', "bc", 0, &[3, 1, 1, 1]),
];

fn tw_id() -> [u8; TW_N] { let mut x = [0u8; TW_N]; for i in 0..TW_N { x[i] = i as u8; } x }

// base lemmas (T2b only for m != 0) times the shifts m in --shifts, distinct mod ord(p)
fn tw_shift_lemmas(base: Vec<TwLemma>, p: Perm, n: usize) -> Vec<TwLemma> {
    let ord = tw_order(p, n) as i64;
    let mut ms: Vec<i64> = vec![0];
    for &m in &tw_opts().shifts {
        if ms.iter().all(|&q| (q - m).rem_euclid(ord) != 0) { ms.push(m); }
    }
    let pinv = p.inverse(n);
    let mut out = Vec::new();
    for &m in &ms {
        let (mut sh, mut shi) = (tw_id(), tw_id());
        for _ in 0..m.unsigned_abs() {
            for i in 0..n { sh[i] = if m > 0 { p.img[sh[i] as usize] } else { pinv.img[sh[i] as usize] }; }
        }
        for i in 0..n { shi[sh[i] as usize] = i as u8; }
        for l in &base {
            if m == 0 && l.fam == 8 { continue; }
            // shifted copies that never pruned on the test frames (they coincide with other members)
            if m != 0 && (matches!(l.name, "S404" | "S0444") || (l.name == "S0044" && m != 1)
                || (l.name == "R3" && m < 0) || (l.name == "R4" && m > 0)) { continue; }
            let mut q = l.clone();
            q.shift = m; q.sh = sh; q.shi = shi;
            out.push(q);
        }
    }
    out
}

fn tw_letters(s: &str) -> Vec<u8> {
    s.chars().map(|ch| match ch { 'a' => 0, 'A' => 1, 'b' => 2, 'B' => 3, 'p' => 4, 'P' => 5, 'e' => 6, 'E' => 7, 'c' => 8, 'C' => 9, _ => panic!("letter") }).collect()
}
fn tw_inverse_word(w: &[u8]) -> Vec<u8> { w.iter().rev().map(|&l| l ^ 1).collect() }

fn tw_lemmas(names: &str) -> Vec<TwLemma> {
    let want: Vec<&str> = names.split(',').filter(|s| !s.is_empty()).collect();
    TW_LEMMA_SPECS.iter().filter(|s| names == "all" || want.contains(&s.0)).map(|&(name, x, pre, y, other, fam, word)| {
        let ent = |ch: char| match ch { 'b' => 1u8, 'c' => 2, 'e' => 3, _ => panic!("entry") };
        let ox: Vec<char> = other.chars().collect();
        TwLemma { name, x: tw_letters(x), xi: tw_inverse_word(&tw_letters(x)), pre: tw_letters(pre),
            y: match y { 'b' => 0, 'e' => 1, _ => 2 }, other: [ent(ox[0]), ent(ox[1])], fam, word, shift: 0, sh: tw_id(), shi: tw_id() }
    }).collect()
}

// literal source word of lemma member j (j != 0 for families)
fn tw_lemma_word(l: &TwLemma, j: i64) -> Vec<usize> {
    let k = j.unsigned_abs() as usize;
    let mut w: Vec<usize> = vec![if l.shift > 0 { 0 } else { 3 }; 2 * l.shift.unsigned_abs() as usize];
    match l.fam {
        0 => w.extend(l.word.iter().copied()),
        9 => if j > 0 { w.extend(vec![1usize; 2 * k]); w.push(4); } else { w.extend(vec![4usize; 2 * k + 1]); },
        8 => w.extend(vec![if j > 0 { 1usize } else { 4 }; 2 * k]),
        f => { w.push((f - 1) as usize); w.extend(std::iter::repeat(if j > 0 { 1 } else { 4 }).take(2 * k)); }
    }
    w
}

// lemmas whose other cusp half has index sum >= n for every allowed ind c
fn tw_gated_lemmas(n: usize, ib: usize, ie: usize, ics: &[usize]) -> Vec<TwLemma> {
    let opts = tw_opts();
    tw_lemmas(&opts.lemmas).into_iter().filter(|l| ics.iter().all(|&ic| {
        let ind = |x: u8| match x { 1 => ib, 2 => ic, _ => ie };
        ind(l.other[0]) + ind(l.other[1]) >= n
    })).collect()
}

// Does any gated lemma's literal word fail on the full tuple t (members j up to jmax)?
fn tw_lemma_literal(t: [Perm; 4], n: usize, lem: &[TwLemma]) -> Option<Vec<usize>> {
    let jmax = tw_opts().jmax as i64;
    for l in lem {
        if l.fam == 0 {
            let w = tw_lemma_word(l, 0);
            if tw_replay_fails(t, n, &w) { return Some(w); }
            continue;
        }
        for j in (1..=jmax).chain((1..=jmax).map(|k| -k)) {
            let w = tw_lemma_word(l, j);
            if tw_replay_fails(t, n, &w) { return Some(w); }
        }
    }
    None
}

// canonical representative of e under conjugation by the group {1} + sym
fn tw_sym_canon(e: Perm, n: usize, sym: &[[u8; TW_N]]) -> [u8; TW_N] {
    let mut best = [0u8; TW_N];
    for i in 0..n { best[i] = e.img[i]; }
    for h in sym {
        let mut hinv = [0u8; TW_N];
        for i in 0..n { hinv[h[i] as usize] = i as u8; }
        let mut v = [0u8; TW_N];
        for i in 0..n { v[i] = h[e.img[hinv[i] as usize] as usize]; }
        if v[..n] < best[..n] { best = v; }
    }
    best
}

const TW_DEFAULT_LEMMAS: &str = "T2c,T2b,R3,R4,S131,S404,S0044,S0404,S0444,S1131,S1311,S3111";
struct TwOpts { sym: bool, lemmas: String, jmax: usize, sym_cap: usize, drive: bool, shifts: Vec<i64>, lazy_pct: usize, b3_depth: usize }
// shifted lemmas are tested in the search down to this depth (percent of n), always at leaves
fn tw_lazy_depth(n: usize) -> usize { n * tw_opts().lazy_pct / 100 }
static TW_OPTS: std::sync::OnceLock<TwOpts> = std::sync::OnceLock::new();
fn tw_opts() -> &'static TwOpts {
    TW_OPTS.get_or_init(|| TwOpts { sym: true, lemmas: TW_DEFAULT_LEMMAS.into(), jmax: 64, sym_cap: 5040, drive: true, shifts: Vec::new(), lazy_pct: 60, b3_depth: 0 })
}

// A subgroup H of the centralizer of <a,b> in S_n, as an explicit element list
// (identity excluded). Generators: automorphisms of each orbit and swaps of
// isomorphic orbits; a generator whose closure would exceed `cap` is dropped,
// so the result is always a group.
fn tw_centralizer(gens: &[Perm], n: usize, cap: usize) -> Vec<[u8; TW_N]> {
    let mut gi: Vec<Perm> = gens.to_vec();
    for g in gens { gi.push(g.inverse(n)); }
    let mut seen = [false; TW_N];
    let mut orbs: Vec<Vec<u8>> = Vec::new();
    for x in 0..n {
        if seen[x] { continue; }
        let mut o = vec![x as u8];
        seen[x] = true;
        let mut k = 0;
        while k < o.len() {
            let y = o[k] as usize; k += 1;
            for g in &gi { let z = g.img[y] as usize; if !seen[z] { seen[z] = true; o.push(z as u8); } }
        }
        orbs.push(o);
    }
    // the map base(o) -> t extended along generators, if consistent and injective
    let ext = |o: u8, t: u8| -> Option<[u8; TW_N]> {
        let mut m = [u8::MAX; TW_N];
        let mut hit = [false; TW_N];
        m[o as usize] = t; hit[t as usize] = true;
        let mut stack = vec![o];
        while let Some(y) = stack.pop() {
            for g in &gi {
                let (z, w) = (g.img[y as usize], g.img[m[y as usize] as usize]);
                if m[z as usize] == u8::MAX {
                    if hit[w as usize] { return None; }
                    m[z as usize] = w; hit[w as usize] = true; stack.push(z);
                } else if m[z as usize] != w { return None; }
            }
        }
        Some(m)
    };
    let mut idp = [0u8; TW_N];
    for i in 0..TW_N { idp[i] = i as u8; }
    let mut gens_h: Vec<[u8; TW_N]> = Vec::new();
    for (i, o) in orbs.iter().enumerate() {
        for &t in o.iter().skip(1) {
            if let Some(m) = ext(o[0], t) {
                let mut h = idp;
                for &x in o { h[x as usize] = m[x as usize]; }
                gens_h.push(h);
            }
        }
        // swap with the next isomorphic orbit
        if let Some(j) = (i + 1..orbs.len()).find(|&j| orbs[j].len() == o.len() && orbs[j].iter().any(|&t| ext(o[0], t).is_some())) {
            let t = orbs[j].iter().copied().find(|&t| ext(o[0], t).is_some()).unwrap();
            let m = ext(o[0], t).unwrap();
            let mut h = idp;
            for &x in o { h[x as usize] = m[x as usize]; h[m[x as usize] as usize] = x; }
            gens_h.push(h);
        }
    }
    let mul = |g: &[u8; TW_N], h: &[u8; TW_N]| { let mut r = idp; for i in 0..n { r[i] = g[h[i] as usize]; } r };
    let mut group: Vec<[u8; TW_N]> = vec![idp];
    let mut set: HashSet<[u8; TW_N]> = group.iter().copied().collect();
    let mut kept: Vec<[u8; TW_N]> = Vec::new();
    for g in gens_h {
        if set.contains(&g) { continue; }
        let mut trial = kept.clone();
        trial.push(g);
        let mut grp: Vec<[u8; TW_N]> = vec![idp];
        let mut s: HashSet<[u8; TW_N]> = grp.iter().copied().collect();
        let mut k = 0;
        let mut over = false;
        while k < grp.len() {
            let x = grp[k]; k += 1;
            for q in &trial {
                let y = mul(q, &x);
                if s.insert(y) { grp.push(y); if grp.len() > cap { over = true; break; } }
            }
            if over { break; }
        }
        if !over { kept = trial; group = grp; set = s; }
    }
    group.into_iter().filter(|h| *h != idp).collect()
}

struct ArcState {
    e: [u8; TW_N],
    used: [bool; TW_N],
    einv: [u8; TW_N],
    ufs: Vec<[u8; TW_N]>, // 0..K: e-forests, K..2K: c-forests
    comp: [[u8; TW_N]; 2], // e-components, c-components
    hyp: [[u8; TW_N]; 3],  // b-cycle hyperforests: j = 1, j = -1, and c(A_i) (T3)
    walked: bool,
    closed: [usize; 2],
    open: [usize; 2],
    img_done: [u8; TW_N],
    pre_done: [u8; TW_N],
    log: Vec<(u16, u8, u8)>,
    kills: Vec<u64>,
    sym_skip: u64,
    nodes: u64,
    need: u16, // first unknown value a lemma met: x+1 = position e(x), 0x100+y = preimage e^-1(y)
    // per lemma: bits of the unknowns its last evaluation stopped at (positions 0..31,
    // preimages 32..63). A lemma is re-evaluated only if the new arc fills one of them
    // or merges the Y-components it reads; otherwise its result cannot change.
    lmask: Vec<u64>,
    mlog: Vec<(u8, u64)>,
    merged: [bool; 2],
    depth: usize, // number of arcs after the current try_arc
}

fn arc_find(uf: &[u8; TW_N], mut x: u8) -> u8 { while uf[x as usize] != x { x = uf[x as usize]; } x }

impl<'a> ArcSearch<'a> {
    // array ids: 0..2K forests, 2K+0/1 comp, 2K+2/3 hyp
    fn set(&self, st: &mut ArcState, id: usize, idx: u8, val: u8) {
        let k2 = 2 * self.labs.len();
        let arr = if id < k2 { &mut st.ufs[id] } else if id < k2 + 2 { &mut st.comp[id - k2] } else { &mut st.hyp[id - k2 - 2] };
        st.log.push((id as u16, idx, arr[idx as usize]));
        arr[idx as usize] = val;
    }
    fn undo(&self, st: &mut ArcState, mark: usize) {
        let k2 = 2 * self.labs.len();
        while st.log.len() > mark {
            let (id, idx, old) = st.log.pop().unwrap();
            let id = id as usize;
            let arr = if id < k2 { &mut st.ufs[id] } else if id < k2 + 2 { &mut st.comp[id - k2] } else { &mut st.hyp[id - k2 - 2] };
            arr[idx as usize] = old;
        }
    }
    // side 0: e-arc, side 1: c-arc
    fn arc(&self, st: &mut ArcState, side: usize, x: u8, y: u8) -> bool {
        let k = self.labs.len();
        let k2 = 2 * k;
        let (cx, cy) = (arc_find(&st.comp[side], x), arc_find(&st.comp[side], y));
        st.merged[side] = cx != cy;
        if cx == cy {
            st.closed[side] += 1;
            st.open[side] -= 1;
        } else {
            for t in 0..k {
                let uf = &st.ufs[side * k + t];
                if arc_find(uf, self.labs[t][x as usize]) == arc_find(uf, self.labs[t][y as usize]) { return false; }
            }
            for t in 0..k {
                let id = side * k + t;
                let (u, v) = (arc_find(&st.ufs[id], self.labs[t][x as usize]), arc_find(&st.ufs[id], self.labs[t][y as usize]));
                self.set(st, id, u, v);
            }
            self.set(st, k2 + side, cx, cy);
            st.open[side] -= 1;
        }
        let (lo, hi) = if side == 0 { (self.m, self.m) } else { (self.tmin, self.tmax) };
        st.closed[side] <= hi && st.closed[side] + st.open[side] >= lo
    }
    // hyperedge of b-cycle labels of `pts` in hyperforest h
    fn hyper(&self, st: &mut ArcState, h: usize, pts: &[u8]) -> bool {
        let k2 = 2 * self.labs.len();
        let mut roots = [0u8; TW_N];
        for (t, &x) in pts.iter().enumerate() {
            let rt = arc_find(&st.hyp[h], self.blab[x as usize]);
            if roots[..t].contains(&rt) { return false; }
            roots[t] = rt;
        }
        for t in 1..pts.len() { self.set(st, k2 + 2 + h, roots[t], roots[0]); }
        true
    }
    fn go(&self, st: &mut ArcState, t: usize, emit: &mut dyn FnMut(Perm)) {
        let all: Vec<u16> = (0..self.sym.len() as u16).collect();
        self.go_sym(st, t, &all, emit);
    }
    // `sym` lists the elements of H fixing every point already used by an arc
    // (a subgroup). At point x only one value y per orbit of Stab_sym(x) is
    // tried: y and h(y) give subtrees exchanged by h, hence the same orbits.
    fn go_sym(&self, st: &mut ArcState, t: usize, sym: &[u16], emit: &mut dyn FnMut(Perm)) {
        let n = self.n;
        st.nodes += 1;
        let walked = st.walked;
        let x = match self.next_point(st, t) {
            TwNext::Fail => { st.walked = walked; return; }
            TwNext::Leaf => {
                if st.open == [0, 0] && st.closed[0] == self.m && st.closed[1] >= self.tmin && st.closed[1] <= self.tmax && self.leaf_ok(st) {
                    let mut e = Perm::identity(n);
                    for i in 0..n { e.img[i] = st.e[i]; }
                    emit(e);
                }
                st.walked = walked;
                return;
            }
            TwNext::Point(x) => x,
            TwNext::Value(y) => {
                // branch on the preimage x of y
                let ty: Vec<u16> = sym.iter().copied().filter(|&h| self.sym[h as usize][y as usize] == y).collect();
                for x in 0..n as u8 {
                    if st.e[x as usize] != u8::MAX { continue; }
                    if ty.iter().any(|&h| self.sym[h as usize][x as usize] < x) { st.sym_skip += 1; continue; }
                    let child: Vec<u16> = ty.iter().copied().filter(|&h| self.sym[h as usize][x as usize] == x).collect();
                    let mark = st.log.len(); let mm = st.mlog.len();
                    let (sc, so, py) = (st.closed, st.open, st.pre_done);
                    st.depth = t + 1; if self.try_arc(st, x, y) { self.go_sym(st, t + 1, &child, emit); }
                    self.undo(st, mark); while st.mlog.len() > mm { let (li, old) = st.mlog.pop().unwrap(); st.lmask[li as usize] = old; }
                    st.closed = sc; st.open = so; st.pre_done = py;
                    st.used[y as usize] = false;
                    st.e[x as usize] = u8::MAX;
                }
                st.walked = walked;
                return;
            }
        };
        let tx: Vec<u16> = sym.iter().copied().filter(|&h| self.sym[h as usize][x as usize] == x).collect();
        for y in 0..n as u8 {
            if st.used[y as usize] { continue; }
            if tx.iter().any(|&h| self.sym[h as usize][y as usize] < y) { st.sym_skip += 1; continue; }
            let child: Vec<u16> = tx.iter().copied().filter(|&h| self.sym[h as usize][y as usize] == y).collect();
            let mark = st.log.len(); let mm = st.mlog.len();
            let (sc, so, py) = (st.closed, st.open, st.pre_done);
            st.depth = t + 1; if self.try_arc(st, x, y) { self.go_sym(st, t + 1, &child, emit); }
            self.undo(st, mark); while st.mlog.len() > mm { let (li, old) = st.mlog.pop().unwrap(); st.lmask[li as usize] = old; }
            st.closed = sc; st.open = so; st.pre_done = py;
            st.used[y as usize] = false;
            st.e[x as usize] = u8::MAX;
        }
        st.walked = walked;
    }
    // letter l applied to point x; u8::MAX if e is not yet known there
    #[inline]
    fn ev1(&self, st: &ArcState, l: u8, x: u8) -> u8 {
        if x == u8::MAX { return x; }
        match l {
            0..=5 => self.lt[l as usize].img[x as usize],
            6 => st.e[x as usize],
            7 => if st.used[x as usize] { st.einv[x as usize] } else { u8::MAX },
            8 => if st.used[x as usize] { self.pinv.img[st.einv[x as usize] as usize] } else { u8::MAX },
            _ => st.e[self.lt[4].img[x as usize] as usize],
        }
    }
    #[inline]
    fn ev(&self, st: &ArcState, w: &[u8], mut x: u8) -> u8 {
        for &l in w.iter().rev() { x = self.ev1(st, l, x); if x == u8::MAX { break; } }
        x
    }
    // like ev, and records in `need` what was missing first
    #[inline]
    fn evn(&self, st: &ArcState, w: &[u8], mut x: u8, mask: &mut u64) -> u8 {
        if x == u8::MAX { return x; }
        for &l in w.iter().rev() {
            let y = self.ev1(st, l, x);
            if y == u8::MAX {
                *mask |= match l { 6 => 1u64 << x, 9 => 1u64 << self.lt[4].img[x as usize], _ => 1u64 << (32 + x as u32) };
                return y;
            }
            x = y;
        }
        x
    }
    // word w of lemma l (with its shift) at x, recording missing unknowns in mask
    #[inline]
    fn evl(&self, st: &ArcState, l: &TwLemma, w: &[u8], mut x: u8, mask: &mut u64) -> u8 {
        if l.shift == 0 { return self.evn(st, w, x, mask); }
        if x == u8::MAX { return x; }
        for &q in w.iter().rev() {
            if q < 6 { x = self.lt[q as usize].img[x as usize]; continue; }
            let z = l.sh[x as usize];
            let y = self.ev1(st, q, z);
            if y == u8::MAX {
                *mask |= match q { 6 => 1u64 << z, 9 => 1u64 << self.lt[4].img[z as usize], _ => 1u64 << (32 + z as u32) };
                return y;
            }
            x = l.shi[y as usize];
        }
        x
    }
    // hyperforest test of the known points cur[u] (u over the a-cycles) on Y-cycles
    fn lem_forest(&self, st: &ArcState, y: u8, cur: &[u8; TW_N], sh: &[u8; TW_N]) -> bool {
        if self.acyc.len() == 1 {
            // one cycle of a: the known points must lie in distinct Y-cycles
            let mut seen = 0u32;
            for &u in &self.acyc[0] {
                let v = cur[u as usize];
                if v == u8::MAX { continue; }
                let lab = match y { 0 => self.blab[v as usize], 1 => arc_find(&st.comp[0], sh[v as usize]), _ => arc_find(&st.comp[1], sh[v as usize]) };
                let bit = 1u32 << lab;
                if seen & bit != 0 { return false; }
                seen |= bit;
            }
            return true;
        }
        let mut uf = [0u8; TW_N];
        for i in 0..TW_N { uf[i] = i as u8; }
        for cyc in &self.acyc {
            let mut roots = [0u8; TW_N];
            let mut k = 0;
            for &u in cyc {
                let v = cur[u as usize];
                if v == u8::MAX { continue; }
                let lab = match y { 0 => self.blab[v as usize], 1 => arc_find(&st.comp[0], sh[v as usize]), _ => arc_find(&st.comp[1], sh[v as usize]) };
                let rt = arc_find(&uf, lab);
                if roots[..k].contains(&rt) { return false; }
                roots[k] = rt; k += 1;
            }
            for i in 1..k { uf[roots[i] as usize] = roots[0]; }
        }
        true
    }
    // one lemma on the current partial e: (passes, unknowns it stopped at)
    fn eval_lemma(&self, st: &ArcState, l: &TwLemma) -> (bool, u64) {
        let jmax = tw_opts().jmax;
        let mut mask = 0u64;
        let mut cur = [u8::MAX; TW_N];
        for cyc in &self.acyc { for &u in cyc { cur[u as usize] = self.evl(st, l, &l.pre, u, &mut mask); } }
        if l.x.is_empty() { return (self.lem_forest(st, l.y, &cur, &l.sh), mask); }
        let start = cur;
        for (di, dir) in [&l.x, &l.xi].into_iter().enumerate() {
            let mut c = start;
            for _ in 0..jmax {
                let mut any = false;
                let mut back = true;
                for cyc in &self.acyc { for &u in cyc {
                    let v = self.evl(st, l, dir, c[u as usize], &mut mask);
                    c[u as usize] = v;
                    if v != u8::MAX { any = true; }
                    if v != start[u as usize] { back = false; }
                } }
                // a full period forward already covers every backward member
                if back { if di == 0 { return (true, mask); } break; }
                if !any { break; }
                if !self.lem_forest(st, l.y, &c, &l.sh) { return (false, mask); }
            }
        }
        (true, mask)
    }
    // all gated lemmas on the current partial e after the arc x -> y; counts the first that fails.
    // Shifted (lazy) lemmas run only at depth <= lazy_depth here, and on every leaf (leaf_ok).
    fn lemmas_ok(&self, st: &mut ArcState, x: u8, y: u8) -> bool {
        let hitbits = (1u64 << x) | (1u64 << (32 + y as u32));
        for (li, l) in self.lem.iter().enumerate() {
            if l.shift != 0 {
                if st.depth > self.lazy_depth { continue; }
                if !self.eval_lemma(st, l).0 { st.kills[li] += 1; return false; }
                continue;
            }
            let comp_change = (l.y == 1 && st.merged[0]) || (l.y == 2 && st.merged[1]);
            if st.lmask[li] & hitbits == 0 && !comp_change { continue; }
            let (ok, mask) = self.eval_lemma(st, l);
            if !ok { st.kills[li] += 1; return false; }
            if mask != st.lmask[li] { st.mlog.push((li as u8, st.lmask[li])); st.lmask[li] = mask; }
        }
        true
    }
    // complete e: the lazy lemmas that were skipped below lazy_depth
    fn leaf_ok(&self, st: &mut ArcState) -> bool {
        if st.depth <= self.lazy_depth { return true; }
        for (li, l) in self.lem.iter().enumerate() {
            if l.shift != 0 && !self.eval_lemma(st, l).0 { st.kills[li] += 1; return false; }
        }
        true
    }
    // next unknown wanted by the first lemma that has one
    fn lemma_need(&self, st: &ArcState) -> u16 {
        for &m in &st.lmask {
            if m != 0 {
                let b = m.trailing_zeros() as u16;
                return if b < 32 { b + 1 } else { 0x100 + (b - 32) };
            }
        }
        0
    }
}


enum TwNext { Fail, Point(u8), Value(u8), Leaf }

impl<'a> ArcSearch<'a> {
    // T2 walk: positions r^j(u), u in U, r(x) = e(a(x)). Each fully known rotation
    // r^j(A_i) must be a hyperforest on the b-cycles (lemma T2 for every j). Returns
    // the first point whose e-value the walk needs, or that the walk has closed.
    fn walk(&self, st: &ArcState) -> Option<Option<u8>> {
        let n = self.n;
        let mut cur = [0u8; TW_N];
        for x in 0..n { cur[x] = x as u8; }
        loop {
            let mut next = cur;
            for cyc in &self.acyc {
                for &u in cyc {
                    let y = self.a.img[cur[u as usize] as usize];
                    if st.e[y as usize] == u8::MAX { return Some(Some(y)); }
                    next[u as usize] = st.e[y as usize];
                }
            }
            let mut uf = [0u8; TW_N];
            for i in 0..TW_N { uf[i] = i as u8; }
            let mut back = true;
            for cyc in &self.acyc {
                let mut roots = [0u8; TW_N];
                for (t, &u) in cyc.iter().enumerate() {
                    let rt = arc_find(&uf, self.blab[next[u as usize] as usize]);
                    if roots[..t].contains(&rt) { return None; }
                    roots[t] = rt;
                    if next[u as usize] != u { back = false; }
                }
                for t in 1..cyc.len() { uf[roots[t] as usize] = roots[0]; }
            }
            if back { return Some(None); }
            cur = next;
        }
    }
    fn next_point(&self, st: &mut ArcState, t: usize) -> TwNext {
        // the walk must finish (and pass) before any leaf, even when it consumed every point
        if !st.walked {
            match self.walk(st) {
                None => return TwNext::Fail,
                Some(Some(x)) => return TwNext::Point(x),
                Some(None) => st.walked = true,
            }
        }
        if t == self.n { return TwNext::Leaf; }
        let need = if tw_opts().drive { self.lemma_need(st) } else { 0 };
        if need != 0 {
            if need < 0x100 { let x = (need - 1) as u8; if st.e[x as usize] == u8::MAX { return TwNext::Point(x); } }
            else { let y = (need - 0x100) as u8; if !st.used[y as usize] { return TwNext::Value(y); } }
        }
        match self.order.iter().copied().find(|&x| st.e[x as usize] == u8::MAX) {
            Some(x) => TwNext::Point(x),
            None => TwNext::Leaf,
        }
    }
    // one branch e(x) = y with every incremental test; false = prune
    fn try_arc(&self, st: &mut ArcState, x: u8, y: u8) -> bool {
        st.e[x as usize] = y;
        st.used[y as usize] = true;
        st.einv[y as usize] = x;
        let mut ok = self.arc(st, 0, x, y) && self.arc(st, 1, y, self.pinv.img[x as usize]);
        let cj = self.acyc_of[y as usize];
        if ok && cj != u8::MAX {
            st.pre_done[cj as usize] += 1;
            let cyc = &self.acyc[cj as usize];
            if st.pre_done[cj as usize] as usize == cyc.len() {
                // j = -1: r^-1(A_i) = a^-1 e^-1(A_i); T3: c(A_i) = p^-1 e^-1(A_i)
                let pts: Vec<u8> = cyc.iter().map(|&u| self.ainv.img[st.einv[u as usize] as usize]).collect();
                let cpts: Vec<u8> = cyc.iter().map(|&u| self.pinv.img[st.einv[u as usize] as usize]).collect();
                ok = self.hyper(st, 1, &pts) && self.hyper(st, 2, &cpts);
            }
        }
        ok && (self.lem.is_empty() || self.lemmas_ok(st, x, y))
    }
}


// T3: sigma2^-1 then sigma1 then the rotation; cusp holds (a, c^-1 b c) | (e, .).
fn tw_t3_word(t: [Perm; 4], n: usize) -> Option<Vec<usize>> {
    let (a, b, c) = (t[0], t[1], t[2]);
    let x = c.inverse(n).compose(b, n).compose(c, n);
    if sm_index(a.compose(x, n), n) != sm_index(a, n) + sm_index(b, n) { Some(vec![4, 0, 0, 1, 2]) } else { None }
}


fn tw_smooth_here(t: &[Perm; 4], n: usize) -> bool {
    let ip = sm_index(t[0].compose(t[1], n), n);
    ip == sm_index(t[0], n) + sm_index(t[1], n) || ip == sm_index(t[2], n) + sm_index(t[3], n)
}

// Literal depth-limited search over the six moves: the first word (length <= depth)
// reaching a state whose cusp (12)|(34) is not smooth. No hashing or canonical keys.
fn tw_short_dfs(t: [Perm; 4], n: usize, depth: usize) -> Option<Vec<usize>> {
    fn rec(t: [Perm; 4], n: usize, depth: usize, word: &mut Vec<usize>) -> bool {
        if !tw_smooth_here(&t, n) { return true; }
        if depth == 0 { return false; }
        for mv in 0..6 {
            // skip a move that undoes the previous one
            if let Some(&last) = word.last() { if (last + 3) % 6 == mv { continue; } }
            word.push(mv);
            if rec(hurwitz_move_dynamic(n, t, mv), n, depth - 1, word) { return true; }
            word.pop();
        }
        false
    }
    let mut w = Vec::new();
    if rec(t, n, depth, &mut w) { Some(w) } else { None }
}

// Literal search over words in sigma1^+-1, sigma2^+-1 only (B3 on positions 1..3, e fixed),
// length <= depth: the first word reaching a non-smooth cusp (12)|(34).
fn tw_b3_dfs(t: [Perm; 4], n: usize, depth: usize) -> Option<Vec<usize>> {
    fn rec(t: [Perm; 4], n: usize, depth: usize, word: &mut Vec<usize>) -> bool {
        if !tw_smooth_here(&t, n) { return true; }
        if depth == 0 { return false; }
        for mv in [0usize, 1, 3, 4] {
            if let Some(&last) = word.last() { if (last + 3) % 6 == mv { continue; } }
            word.push(mv);
            if rec(hurwitz_move_dynamic(n, t, mv), n, depth - 1, word) { return true; }
            word.pop();
        }
        false
    }
    let mut w = Vec::new();
    if rec(t, n, depth, &mut w) { Some(w) } else { None }
}

// T2 at j = 1 and j = -1 for a full tuple (the predicate pruned by ArcSearch).
fn tw_t2b_pm1_word(t: [Perm; 4], n: usize) -> Option<Vec<usize>> {
    let r = t[3].compose(t[0], n);
    let rinv = r.inverse(n);
    let ia_ib = sm_index(t[0], n) + sm_index(t[1], n);
    let ord = tw_order(r, n);
    for (j, h) in [(1usize, r), (ord - 1, rinv)] {
        let aj = h.compose(t[0], n).compose(h.inverse(n), n);
        if sm_index(aj.compose(t[1], n), n) != ia_ib { return Some(vec![1usize; 2 * j]); }
    }
    None
}

// T1 literal word for the excluded e: 2^(2k) then 0,1,2, cusp (12)|(34) fails.
fn tw_t1_word(t: [Perm; 4], n: usize, p: Perm) -> Option<Vec<usize>> {
    let mut ak = t[0];
    for k in 0..tw_order(p, n) {
        if !tw_additive(t[3], ak, n) {
            let mut w = vec![2usize; 2 * k];
            w.extend([0, 1, 2]);
            return Some(w);
        }
        ak = ak.conjugate_by(p, n); // ak = p^k a p^-k, so the failing state is sigma3^{2k}
    }
    None
}

fn tw_replay_fails(t: [Perm; 4], n: usize, word: &[usize]) -> bool {
    let mut u = t;
    for &m in word { u = hurwitz_move_dynamic(n, u, m); }
    !sm_values(&u, n).0
}

#[allow(clippy::too_many_arguments)]
fn tw_frame(
    n: usize, fi: usize, key: [u128; 2], patterns: &[[usize; 4]], out: &Path, proofs: &Mutex<StarProofs>,
    cache_limit: usize, cap: usize, closure_seconds: u64, brute: Option<&[(Perm, usize, Perm)]>, all_partners: bool,
    primitive_only: bool,
) -> Result<([u128; 9], String, String), String> {
    let mut imprimitive = 0u128;
    let mut t2_pruned = 0u128;
    let mut brute_note = String::from("-");
    let (a, b) = tw_unpack(key, n);
    let p = a.compose(b, n);
    let pinv = p.inverse(n);
    let (ia, ib) = (sm_index(a, n), sm_index(b, n));
    assert_eq!(ia + ib, sm_index(p, n));
    // allowed (ie, ic): remaining two entries of a matching pattern
    let mut allowed: Vec<(usize, usize)> = Vec::new();
    for q in patterns.iter().filter(|q| q[0] == ia) {
        // Normal form (default): ind b = q3, ind e = q2, ind c = q1. Every orbit has
        // such a state: sigma2^-1 or sigma3^-1 sigma2^-1 brings the largest entry to
        // position 2, then sigma3 puts the larger remaining index in position 4.
        if !all_partners {
            if q[3] == ib && !allowed.contains(&(q[2], q[1])) { allowed.push((q[2], q[1])); }
            continue;
        }
        let mut rest = vec![q[1], q[2], q[3]];
        let Some(pos) = rest.iter().position(|&x| x == ib) else { continue };
        rest.remove(pos);
        for (ie, ic) in [(rest[0], rest[1]), (rest[1], rest[0])] {
            if !allowed.contains(&(ie, ic)) { allowed.push((ie, ic)); }
        }
    }
    // distinct p-rotations of the cycle partition of a
    let mut labs: Vec<[u8; TW_N]> = Vec::new();
    let mut ak = a;
    for _ in 0..tw_order(p, n) {
        let l = tw_labels(ak, n);
        if !labs.contains(&l) { labs.push(l); }
        ak = ak.conjugate_by(p, n);
    }
    let mut count = [0u128; 9];
    let mut words: HashMap<Vec<usize>, u64> = HashMap::new();
    let mut enumerated: HashSet<Perm> = HashSet::new();
    let mut ies: Vec<usize> = allowed.iter().map(|x| x.0).collect();
    ies.sort_unstable(); ies.dedup();
    let mut err: Option<String> = None;
    let symg = tw_centralizer(&[a, b], n, tw_opts().sym_cap);
    let (mut nodes, mut sym_skip) = (0u64, 0u64);
    let (mut b3_kills, mut w_survivors) = (0u64, 0u64);
    let mut kills: Vec<(String, u64)> = Vec::new();
    for &ie in &ies {
        let m = n - ie;
        let ics: Vec<usize> = allowed.iter().filter(|x| x.0 == ie).map(|x| x.1).collect();
        let (tmin, tmax) = (n - ics.iter().max().unwrap(), n - ics.iter().min().unwrap());
        let mut order: Vec<u8> = (0..n as u8).filter(|&x| a.img[x as usize] != x).collect();
        order.extend((0..n as u8).filter(|&x| a.img[x as usize] == x));
        let mut acyc: Vec<Vec<u8>> = Vec::new();
        let mut acyc_of = [u8::MAX; TW_N];
        {
            let mut seen = [false; TW_N];
            for x in 0..n {
                if seen[x] || a.img[x] as usize == x { continue; }
                let mut cyc = Vec::new();
                let mut y = x;
                while !seen[y] { seen[y] = true; acyc_of[y] = acyc.len() as u8; cyc.push(y as u8); y = a.img[y] as usize; }
                acyc.push(cyc);
            }
        }
        let lem = tw_shift_lemmas(tw_gated_lemmas(n, ib, ie, &ics), p, n);
        let sym = if tw_opts().sym { symg.clone() } else { Vec::new() };
        let srch = ArcSearch { n, m, tmin, tmax, labs: &labs, pinv, a, ainv: a.inverse(n), order, acyc_of, acyc, blab: tw_labels(b, n),
            lt: [a, a.inverse(n), b, b.inverse(n), p, pinv], lem, sym, lazy_depth: tw_lazy_depth(n) };
        let mut id = [0u8; TW_N];
        for i in 0..TW_N { id[i] = i as u8; }
        let mut st = ArcState { e: [u8::MAX; TW_N], used: [false; TW_N], einv: [u8::MAX; TW_N], ufs: vec![id; 2 * labs.len()],
            comp: [id; 2], hyp: [id; 3], walked: false, closed: [0; 2], open: [n; 2], img_done: [0; TW_N], pre_done: [0; TW_N], log: Vec::new(),
            kills: vec![0; srch.lem.len()], sym_skip: 0, nodes: 0, need: 0,
            lmask: srch.lem.iter().map(|l| if l.shift == 0 { u64::MAX } else { 0 }).collect(), mlog: Vec::new(), merged: [false; 2], depth: 0 };
        let sstats = [0u128; 2];
        {
            srch.go(&mut st, 0, &mut |e: Perm| {
                if err.is_some() { return; }
                count[1] += 1;
                let c = pinv.compose(e.inverse(n), n);
                let ic = sm_index(c, n);
                if !allowed.contains(&(ie, ic)) || !st_matches(patterns, ia, ib, ic, ie) { return; }
                count[2] += 1;
                if brute.is_some() { enumerated.insert(e); }
                let t = [a, b, c, e];
                if !generators_transitive(n, &t[..3]) { return; }
                count[3] += 1;
                // T2 and its c-version: r = ea, conjugates r^{-j} b r^j and r^{-j} c r^j
                let r = e.compose(a, n);
                let rinv = r.inverse(n);
                let mut word: Option<Vec<usize>> = None;
                let (mut bj, mut cj) = (b, c);
                for j in 0..tw_order(r, n) {
                    if sm_index(a.compose(bj, n), n) != ia + ib { word = Some(vec![1usize; 2 * j]); break; }
                    if sm_index(a.compose(cj, n), n) != ia + ic { let mut w = vec![1usize; 2 * j]; w.push(4); word = Some(w); break; }
                    bj = rinv.compose(bj, n).compose(r, n);
                    cj = rinv.compose(cj, n).compose(r, n);
                }
                if let Some(w) = word {
                    if !tw_replay_fails(t, n, &w) { err = Some(format!("T2 word did not fail: frame {fi} {w:?}")); return; }
                    count[4] += 1;
                    *words.entry(w).or_insert(0) += 1;
                    return;
                }
                if let Some(w) = sc_word(t, n).or_else(|| tw_short_dfs(t, n, 3)) {
                    count[5] += 1;
                    *words.entry(w).or_insert(0) += 1;
                    return;
                }
                let b3d = tw_opts().b3_depth;
                if b3d > 0 {
                    // iterative deepening: the recorded word has the minimal B3 length
                    if let Some(w) = (1..=b3d).find_map(|k| tw_b3_dfs(t, n, k)) {
                        count[5] += 1;
                        b3_kills += 1;
                        *words.entry(w).or_insert(0) += 1;
                        return;
                    }
                }
                // a start that survives every listed family (statement W) reaches this point
                w_survivors += 1;
                if primitive_only && !tw_primitive(&t[..3], n) { imprimitive += 1; return; }
                let sk = sm_key(&t, n);
                let known = match proofs.lock() {
                    Ok(cache) => if cache.closed.contains_key(&sk) { Some(None) } else { cache.bad.get(&sk).cloned().map(Some) },
                    Err(_) => { err = Some("poisoned".into()); return; }
                };
                let status = match known {
                    Some(None) => "closed_smooth",
                    Some(Some(w)) => { *words.entry(w).or_insert(0) += 1; "nonsmooth" }
                    None => {
                        let scan = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(closure_seconds));
                        let mut cache = match proofs.lock() { Ok(c) => c, Err(_) => { err = Some("poisoned".into()); return; } };
                        match scan.status {
                            "nonsmooth" => {
                                *words.entry(scan.word.clone()).or_insert(0) += 1;
                                if cache.bad.len() < cache_limit { cache.bad.insert(sk, scan.word); }
                            }
                            "closed_smooth" => {
                                if !cache.closed.contains_key(&sk) {
                                    let id = cache.components;
                                    if let Err(x) = sm_save(&t, n, &scan, &out.join(format!("component_f{fi}_{id:04}"))) { err = Some(x); return; }
                                    cache.components += 1;
                                    for k in scan.states { cache.closed.insert(k, id); }
                                }
                            }
                            _ => {
                                if let Err(x) = sm_save(&t, n, &scan, &out.join(format!("unresolved_{fi}_{}", count[1]))) { err = Some(x); return; }
                            }
                        }
                        scan.status
                    }
                };
                match status {
                    "nonsmooth" => count[6] += 1,
                    "closed_smooth" => count[7] += 1,
                    _ => count[8] += 1,
                }
            });
        }
        count[0] += sstats[0];
        t2_pruned += sstats[1];
        nodes += st.nodes;
        sym_skip += st.sym_skip;
        for (l, &k) in srch.lem.iter().zip(st.kills.iter()) {
            let key = if l.shift == 0 { l.name.to_string() } else { format!("{}@{}", l.name, l.shift) };
            match kills.iter_mut().find(|x| x.0 == key) { Some(x) => x.1 += k, None => kills.push((key, k)) }
        }
    }
    if let Some(x) = err { return Err(x); }
    // Optional literal cross-check against all of S_n. Eligible = passes T1, T1c,
    // full T2b, T3 and every gated v10 lemma word (literal replay). The
    // enumeration must lie inside the eligible set and meet every symmetry class.
    if let Some(perms) = brute {
        let mut eligible = 0u128;
        let mut replayed = 0u64;
        let mut lem_cache: HashMap<usize, Vec<TwLemma>> = HashMap::new();
        let mut want: HashSet<[u8; TW_N]> = HashSet::new();
        for &(c, ic, cinv) in perms {
            let e = cinv.compose(pinv, n);
            let ie = sm_index(e, n);
            if !allowed.contains(&(ie, ic)) || !st_matches(patterns, ia, ib, ic, ie) { continue; }
            let lem = lem_cache.entry(ie).or_insert_with(|| {
                let ics: Vec<usize> = allowed.iter().filter(|x| x.0 == ie).map(|x| x.1).collect();
                tw_shift_lemmas(tw_gated_lemmas(n, ib, ie, &ics), p, n)
            });
            let t = [a, b, c, e];
            let w = tw_t1_word(t, n, p).or_else(|| tw_t1c_word(t, n, p)).or_else(|| tw_t2b_word(t, n)).or_else(|| tw_t3_word(t, n))
                .or_else(|| tw_lemma_literal(t, n, lem));
            match w {
                None => {
                    eligible += 1;
                    want.insert(tw_sym_canon(e, n, if tw_opts().sym { &symg } else { &[] }));
                }
                Some(w) => {
                    if enumerated.contains(&e) { return Err(format!("frame {fi}: enumerated e fails a literal word {w:?}")); }
                    if replayed < 64 || c.img[0] == 0 {
                        if !tw_replay_fails(t, n, &w) { return Err(format!("frame {fi}: word {w:?} does not fail")); }
                        replayed += 1;
                    }
                }
            }
        }
        let got: HashSet<[u8; TW_N]> = enumerated.iter().map(|&e| tw_sym_canon(e, n, if tw_opts().sym { &symg } else { &[] })).collect();
        if got != want { return Err(format!("frame {fi}: symmetry classes differ: brute {} vs enumerated {}", want.len(), got.len())); }
        if !tw_opts().sym && eligible != count[2] { return Err(format!("frame {fi}: brute eligible {eligible} vs enumerated {}", count[2])); }
        brute_note = format!("ok:{eligible}:{}:{replayed}", want.len());
    }
    let mut ws: Vec<_> = words.into_iter().collect();
    ws.sort();
    let mut wordlines = String::new();
    for (word, k) in ws {
        wordlines.push_str(&format!("{fi}\t{}\t{k}\n", word.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",")));
    }
    let kill_text = if kills.is_empty() { "-".to_string() } else { kills.iter().map(|(k, v)| format!("{k}:{v}")).collect::<Vec<_>>().join(",") };
    let record = format!("{fi}\t{}\t{}\t{imprimitive}\t{brute_note}\t{t2_pruned}\t{nodes}\t{sym_skip}\t{}\t{kill_text}\t{b3_kills}\t{w_survivors}\n", labs.len(),
        count.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"), symg.len() + 1);
    Ok((count, record, wordlines))
}

pub fn twist_census_main(args: &[String]) -> Result<(), String> {
    let _ = TW_OPTS.set(TwOpts {
        sym: !args.iter().any(|x| x == "--no-sym"),
        lemmas: args.windows(2).find(|w| w[0] == "--lemmas").map(|w| w[1].clone()).unwrap_or(TW_DEFAULT_LEMMAS.into()),
        jmax: args.windows(2).find(|w| w[0] == "--jmax").map(|w| w[1].parse().unwrap()).unwrap_or(64),
        sym_cap: args.windows(2).find(|w| w[0] == "--sym-cap").map(|w| w[1].parse().unwrap()).unwrap_or(5040),
        drive: !args.iter().any(|x| x == "--no-drive"),
        b3_depth: args.windows(2).find(|w| w[0] == "--b3-depth").map(|w| w[1].parse().unwrap()).unwrap_or(0),
        lazy_pct: args.windows(2).find(|w| w[0] == "--lazy-pct").map(|w| w[1].parse().unwrap()).unwrap_or(60),
        shifts: args.windows(2).find(|w| w[0] == "--shifts").map(|w| w[1].split(',').filter(|s| !s.is_empty()).map(|x| x.parse().unwrap()).collect()).unwrap_or_default(),
    });
    if args.iter().any(|x| x == "--sample") { return twist_sampler_main(args); }
    let arg = |name: &str, default: &str| args.windows(2).find(|w| w[0] == name)
        .map(|w| w[1].clone()).unwrap_or(default.into());
    let n: usize = arg("--d", "9").parse().map_err(|_| "degree")?;
    let workers: usize = arg("--workers", "6").parse().map_err(|_| "workers")?;
    let cap: usize = arg("--cap", "1000000").parse().map_err(|_| "cap")?;
    let cache_limit: usize = arg("--cache-entries", "2000000").parse().map_err(|_| "cache")?;
    let closure_seconds: u64 = arg("--closure-seconds", "120").parse().map_err(|_| "closure seconds")?;
    let brute = args.iter().any(|x| x == "--brute");
    let all_partners = args.iter().any(|x| x == "--all-partners");
    let primitive_only = args.iter().any(|x| x == "--primitive-only");
    let stride: usize = arg("--brute-stride", "1").parse().map_err(|_| "stride")?;
    let only: Option<usize> = args.windows(2).find(|w| w[0] == "--frame").map(|w| w[1].parse().unwrap());
    let triangle = args.iter().any(|x| x == "--triangle");
    let allow_genus0 = args.iter().any(|x| x == "--allow-genus0");
    let patterns = if triangle { tri_parse(&arg("--patterns", ""), n, allow_genus0)? } else { st_parse_patterns(&arg("--patterns", ""))? };
    if !(3..=25).contains(&n) || !(1..=10).contains(&workers) || (brute && n > 11) {
        return Err("degree3..25, workers1..10, --brute only for degree<=11".into());
    }
    for q in patterns.iter().filter(|_| !triangle) {
        if q.iter().any(|&x| x == 0 || x >= n) { return Err("indices must lie in 1..n-1".into()); }
        let s: usize = q.iter().sum();
        if s % 2 != 0 || s < 2 * n { return Err("pattern must have even sum >= 2n (positive genus)".into()); }
        if q[0] + q[3] > n - 1 { return Err("not a STAR list: smallest plus largest exceeds n-1".into()); }
    }
    let resume = args.iter().any(|x| x == "--resume");
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() { return Err("--out required".into()); }
    if resume && out.exists() {
        if !out.join("PATTERNS.tsv").exists() { return Err("--resume target has no PATTERNS.tsv".into()); }
    } else {
        fs::create_dir(&out).map_err(|e| e.to_string())?;
    }
    let pattern_text: String = patterns.iter()
        .map(|q| format!("{}\t{}\t{}\t{}\tgenus={}\n", q[0], q[1], q[2], q[3], (q.iter().sum::<usize>() + 2 - 2 * n) / 2))
        .collect::<String>();
    if resume && out.join("PATTERNS.tsv").exists() {
        if fs::read_to_string(out.join("PATTERNS.tsv")).map_err(|e| e.to_string())? != pattern_text { return Err("--resume patterns differ".into()); }
    } else {
        fs::write(out.join("PATTERNS.tsv"), &pattern_text).map_err(|e| e.to_string())?;
    }
    // frames already recorded by earlier sessions (one flushed line per finished frame)
    let mut prior: HashMap<usize, [u128; 9]> = HashMap::new();
    if resume {
        for entry in fs::read_dir(&out).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            let name = path.file_name().and_then(|x| x.to_str()).unwrap_or("").to_string();
            if !name.starts_with("COUNTS_") { continue; }
            for line in fs::read_to_string(&path).map_err(|e| e.to_string())?.lines() {
                let v: Vec<&str> = line.split('\t').collect();
                if v.len() < 13 || v[0] == "frame" { continue; }
                let mut c = [0u128; 9];
                for j in 0..9 { c[j] = v[2 + j].parse().map_err(|_| "bad COUNTS line")?; }
                prior.insert(v[0].parse().map_err(|_| "bad frame id")?, c);
            }
        }
    }
    let session = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let start = Instant::now();
    let frames = if triangle {
        let acc = |ia: usize, ib: usize| patterns.iter().any(|q| q[1] == ia && q[2] == ib);
        nc_frames(n, &acc)
    } else { tw_frames(n, &patterns, all_partners) };
    let mut table = BufWriter::new(File::create(out.join("FRAMES.tsv")).map_err(|e| e.to_string())?);
    writeln!(table, "frame\tpair_orbit\ta\tb").map_err(|e| e.to_string())?;
    for (fi, &(key, z)) in frames.iter().enumerate() {
        let (a, b) = tw_unpack(key, n);
        writeln!(table, "{fi}\t{}\t{}\t{}", factorial_u128(n) / z, a.cycle_string(n), b.cycle_string(n)).map_err(|e| e.to_string())?;
    }
    table.flush().map_err(|e| e.to_string())?;
    let perms: Option<Vec<(Perm, usize, Perm)>> = if brute {
        Some(sc_perms(n).into_iter().map(|(p, i)| (p, i, p.inverse(n))).collect())
    } else { None };
    println!("degree={n} workers={workers} triangle={triangle} frames={} normal_form={} primitive_only={primitive_only} setup_seconds={:.3}", frames.len(), !all_partners, start.elapsed().as_secs_f64());
    let todo: Vec<usize> = match only { Some(f) => vec![f], None => (0..frames.len()).filter(|fi| !prior.contains_key(fi)).collect() };
    println!("previously_recorded={} remaining={}", prior.len(), todo.len());
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let proofs = Mutex::new(StarProofs::default());
    let results = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers).map(|w| { let out = &out; let todo = &todo; let next = &next; let done = &done; let proofs = &proofs; let frames = &frames; let patterns = &patterns; let perms = &perms;
          scope.spawn(move || -> Result<Vec<[u128; 9]>, String> {
            let mut result = Vec::new();
            let open = |name: String| OpenOptions::new().create(true).append(true).open(out.join(name)).map_err(|e| e.to_string());
            let mut cf = open(format!("COUNTS_s{session}_w{w}.tsv"))?;
            let mut wf = open(format!("WORDS_s{session}_w{w}.tsv"))?;
            writeln!(cf, "frame\tdistinct_rotations\t{TW_COLUMNS}\timprimitive_skipped\tbrute\tt2_pruned_prefixes\tnodes\tsym_skip\tsym_order\tlemma_kills\tb3_kills\tw_survivors").map_err(|e| e.to_string())?;
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= todo.len() { break; }
                let fi = todo[i];
                let check = if fi % stride == 0 { perms.as_deref() } else { None };
                let (c, record, wordlines) = if triangle {
                    tri_frame(n, fi, frames[fi].0, patterns, out, proofs, cache_limit, cap, closure_seconds, check, primitive_only)?
                } else {
                    tw_frame(n, fi, frames[fi].0, patterns, out, proofs, cache_limit, cap, closure_seconds, check, all_partners, primitive_only)?
                };
                wf.write_all(wordlines.as_bytes()).map_err(|e| e.to_string())?;
                cf.write_all(record.as_bytes()).map_err(|e| e.to_string())?;
                cf.flush().map_err(|e| e.to_string())?;
                result.push(c);
                let completed = done.fetch_add(1, Ordering::Relaxed) + 1;
                if completed % 50 == 0 || completed == todo.len() {
                    println!("frames_processed={completed}/{} elapsed_seconds={:.3}", todo.len(), start.elapsed().as_secs_f64());
                }
            }
            Ok(result)
        })}).collect();
        handles.into_iter().map(|h| h.join().map_err(|_| "worker panicked".to_string()).and_then(|r| r))
            .collect::<Result<Vec<_>, _>>()
    })?;
    let mut total = [0u128; 9];
    let mut completed = 0;
    for c in results.into_iter().flatten().chain(prior.values().copied()) { completed += 1; for j in 0..9 { total[j] += c[j]; } }
    let complete = completed == frames.len() && only.is_none();
    let _ = proofs.into_inner().map_err(|_| "poisoned proof cache")?;
    let components = fs::read_dir(&out).map_err(|e| e.to_string())?.filter(|x| x.as_ref().map(|y| y.file_name().to_string_lossy().starts_with("component_")).unwrap_or(false)).count();
    fs::write(out.join("COVERAGE.tsv"), format!("degree\tframes\tcompleted_frames\t{TW_COLUMNS}\tclosed_components\tenumeration_complete\tfully_resolved\tbrute_checked\telapsed_seconds\n{n}\t{}\t{completed}\t{}\t{components}\t{complete}\t{}\t{brute}\t{:.3}\n",
        frames.len(), total.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"), complete && total[8] == 0, start.elapsed().as_secs_f64()))
        .map_err(|e| e.to_string())?;
    println!("complete={complete} degree={n} counts={total:?} components={components} seconds={:.3}", start.elapsed().as_secs_f64());
    Ok(())
}

// ---------------------------------------------------------------------------
// TRIANGLE regime (--triangle). Sorted list x <= y <= z < h with y + z <= n-1:
// every pair containing the high entry e is non-additive by index.
// Normal form: e in position 4 (rotation) and, by sigma1/sigma2 on the first
// three, ind a = y, ind b = z, ind c = x. Then (a,b) is additive (cusp (12)|(34)).
// Lemma T1-tri: sigma1^{2k} gives (p^k a p^-k, p^k b p^-k, c, e) and its cusp
// (23)|(41) forces (p^k b p^-k, c) additive; sigma1^{2k+1} gives
// (., p^k a p^-k, c, e) and forces (p^k a p^-k, c) additive. Words:
// 0^(2k),0,1,2 and 0^(2k+1),0,1,2. So the cycle partition of c is
// forest-compatible with every p-rotation of a AND of b. c is built arc by arc;
// each c-arc x -> y also fixes the e-arc p(y) -> x (e^{-1} = pc), which is used
// only for the cycle count of e.
fn tri_parse(raw: &str, n: usize, allow_genus0: bool) -> Result<Vec<[usize; 4]>, String> {
    let mut out = Vec::new();
    for part in raw.split(';').filter(|s| !s.trim().is_empty()) {
        let v: Vec<usize> = part.split(',').map(|x| x.trim().parse().map_err(|_| "pattern"))
            .collect::<Result<_, _>>()?;
        if v.len() != 4 { return Err("pattern needs four indices".into()); }
        let mut q = [v[0], v[1], v[2], v[3]];
        q.sort_unstable();
        if q[2] == q[3] { return Err("TRIANGLE list needs a unique largest index".into()); }
        if q[1] + q[2] > n - 1 { return Err("not a TRIANGLE list: y + z > n-1".into()); }
        if q[0] == 0 || q[3] >= n { return Err("indices must lie in 1..n-1".into()); }
        let s: usize = q.iter().sum();
        if s % 2 != 0 || s < 2 * n - 2 || (!allow_genus0 && s < 2 * n) { return Err("bad genus".into()); }
        out.push(q);
    }
    if out.is_empty() { return Err("no patterns".into()); }
    Ok(out)
}

struct TriSearch<'a> {
    n: usize,
    labs: &'a [[u8; TW_N]],
    p: Perm,
    bounds: [(usize, usize); 2], // c cycles, e cycles
}

struct TriState {
    c: [u8; TW_N],
    used: [bool; TW_N],
    ufs: Vec<[u8; TW_N]>,
    comp: [[u8; TW_N]; 2],
    closed: [usize; 2],
    open: [usize; 2],
    log: Vec<(u16, u8, u8)>,
}

impl<'a> TriSearch<'a> {
    fn set(&self, st: &mut TriState, id: usize, idx: u8, val: u8) {
        let k = self.labs.len();
        let arr = if id < k { &mut st.ufs[id] } else { &mut st.comp[id - k] };
        st.log.push((id as u16, idx, arr[idx as usize]));
        arr[idx as usize] = val;
    }
    fn undo(&self, st: &mut TriState, mark: usize) {
        let k = self.labs.len();
        while st.log.len() > mark {
            let (id, idx, old) = st.log.pop().unwrap();
            let id = id as usize;
            let arr = if id < k { &mut st.ufs[id] } else { &mut st.comp[id - k] };
            arr[idx as usize] = old;
        }
    }
    fn arc(&self, st: &mut TriState, side: usize, x: u8, y: u8) -> bool {
        let k = self.labs.len();
        let (cx, cy) = (arc_find(&st.comp[side], x), arc_find(&st.comp[side], y));
        if cx == cy {
            st.closed[side] += 1;
            st.open[side] -= 1;
        } else {
            if side == 0 {
                for t in 0..k {
                    if arc_find(&st.ufs[t], self.labs[t][x as usize]) == arc_find(&st.ufs[t], self.labs[t][y as usize]) { return false; }
                }
                for t in 0..k {
                    let (u, v) = (arc_find(&st.ufs[t], self.labs[t][x as usize]), arc_find(&st.ufs[t], self.labs[t][y as usize]));
                    self.set(st, t, u, v);
                }
            }
            self.set(st, k + side, cx, cy);
            st.open[side] -= 1;
        }
        let (lo, hi) = self.bounds[side];
        st.closed[side] <= hi && st.closed[side] + st.open[side] >= lo
    }
    fn go(&self, st: &mut TriState, x: usize, emit: &mut dyn FnMut(Perm)) {
        let n = self.n;
        if x == n {
            let (lc, hc) = self.bounds[0];
            let (le, he) = self.bounds[1];
            if st.open == [0, 0] && st.closed[0] >= lc && st.closed[0] <= hc && st.closed[1] >= le && st.closed[1] <= he {
                let mut c = Perm::identity(n);
                for i in 0..n { c.img[i] = st.c[i]; }
                emit(c);
            }
            return;
        }
        for y in 0..n as u8 {
            if st.used[y as usize] { continue; }
            let mark = st.log.len();
            let (sc, so) = (st.closed, st.open);
            st.c[x] = y;
            st.used[y as usize] = true;
            if self.arc(st, 0, x as u8, y) && self.arc(st, 1, self.p.img[y as usize], x as u8) { self.go(st, x + 1, emit); }
            self.undo(st, mark);
            st.closed = sc; st.open = so;
            st.used[y as usize] = false;
        }
    }
}

fn tri_t1_word(t: [Perm; 4], n: usize, p: Perm) -> Option<Vec<usize>> {
    let (mut ak, mut bk) = (t[0], t[1]);
    for k in 0..tw_order(p, n) {
        if !tw_additive(bk, t[2], n) { let mut w = vec![0usize; 2 * k]; w.extend([0, 1, 2]); return Some(w); }
        if !tw_additive(ak, t[2], n) { let mut w = vec![0usize; 2 * k + 1]; w.extend([0, 1, 2]); return Some(w); }
        ak = ak.conjugate_by(p, n);
        bk = bk.conjugate_by(p, n);
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn tri_frame(
    n: usize, fi: usize, key: [u128; 2], patterns: &[[usize; 4]], out: &Path, proofs: &Mutex<StarProofs>,
    cache_limit: usize, cap: usize, closure_seconds: u64, brute: Option<&[(Perm, usize, Perm)]>, primitive_only: bool,
) -> Result<([u128; 9], String, String), String> {
    let (a, b) = tw_unpack(key, n);
    let p = a.compose(b, n);
    let (ia, ib) = (sm_index(a, n), sm_index(b, n));
    assert_eq!(ia + ib, sm_index(p, n));
    let allowed: Vec<(usize, usize)> = patterns.iter().filter(|q| q[1] == ia && q[2] == ib).map(|q| (q[0], q[3])).collect();
    let mut labs: Vec<[u8; TW_N]> = Vec::new();
    let (mut ak, mut bk) = (a, b);
    for _ in 0..tw_order(p, n) {
        for l in [tw_labels(ak, n), tw_labels(bk, n)] { if !labs.contains(&l) { labs.push(l); } }
        ak = ak.conjugate_by(p, n);
        bk = bk.conjugate_by(p, n);
    }
    let cyc = |i: usize| n - i;
    let bounds = [
        (allowed.iter().map(|x| cyc(x.0)).min().unwrap(), allowed.iter().map(|x| cyc(x.0)).max().unwrap()),
        (allowed.iter().map(|x| cyc(x.1)).min().unwrap(), allowed.iter().map(|x| cyc(x.1)).max().unwrap()),
    ];
    let srch = TriSearch { n, labs: &labs, p, bounds };
    let mut id = [0u8; TW_N];
    for i in 0..TW_N { id[i] = i as u8; }
    let mut st = TriState { c: [u8::MAX; TW_N], used: [false; TW_N], ufs: vec![id; labs.len()], comp: [id; 2], closed: [0; 2], open: [n; 2], log: Vec::new() };
    let mut count = [0u128; 9];
    let mut imprimitive = 0u128;
    let mut words: HashMap<Vec<usize>, u64> = HashMap::new();
    let mut enumerated: HashSet<Perm> = HashSet::new();
    let mut err: Option<String> = None;
    srch.go(&mut st, 0, &mut |c: Perm| {
        if err.is_some() { return; }
        count[1] += 1;
        let e = p.compose(c, n).inverse(n);
        let (ic, ie) = (sm_index(c, n), sm_index(e, n));
        if !allowed.contains(&(ic, ie)) { return; }
        count[2] += 1;
        if brute.is_some() { enumerated.insert(c); }
        let t = [a, b, c, e];
        if !generators_transitive(n, &t[..3]) { return; }
        count[3] += 1;
        // sigma2^{2j}: (a, B) additive; sigma2^{2j+1}: (a, B C B^-1) additive, B,C = (bc)^j (b,c) (bc)^-j
        let h = b.compose(c, n);
        let hinv = h.inverse(n);
        let (mut bj, mut cj) = (b, c);
        let mut word = None;
        for j in 0..tw_order(h, n) {
            if sm_index(a.compose(bj, n), n) != ia + ib { word = Some(vec![1usize; 2 * j]); break; }
            let x2 = bj.compose(cj, n).compose(bj.inverse(n), n);
            if sm_index(a.compose(x2, n), n) != ia + ic { word = Some(vec![1usize; 2 * j + 1]); break; }
            bj = h.compose(bj, n).compose(hinv, n);
            cj = h.compose(cj, n).compose(hinv, n);
        }
        if let Some(w) = word {
            if !tw_replay_fails(t, n, &w) { err = Some(format!("triangle twist word did not fail: frame {fi} {w:?}")); return; }
            count[4] += 1;
            *words.entry(w).or_insert(0) += 1;
            return;
        }
        if let Some(w) = sc_word(t, n).or_else(|| tw_short_dfs(t, n, 3)) { count[5] += 1; *words.entry(w).or_insert(0) += 1; return; }
        if primitive_only && !tw_primitive(&t[..3], n) { imprimitive += 1; return; }
        let sk = sm_key(&t, n);
        let known = match proofs.lock() {
            Ok(cache) => if cache.closed.contains_key(&sk) { Some(None) } else { cache.bad.get(&sk).cloned().map(Some) },
            Err(_) => { err = Some("poisoned".into()); return; }
        };
        let status = match known {
            Some(None) => "closed_smooth",
            Some(Some(w)) => { *words.entry(w).or_insert(0) += 1; "nonsmooth" }
            None => {
                let scan = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(closure_seconds));
                let mut cache = match proofs.lock() { Ok(x) => x, Err(_) => { err = Some("poisoned".into()); return; } };
                match scan.status {
                    "nonsmooth" => {
                        *words.entry(scan.word.clone()).or_insert(0) += 1;
                        if cache.bad.len() < cache_limit { cache.bad.insert(sk, scan.word); }
                    }
                    "closed_smooth" => {
                        if !cache.closed.contains_key(&sk) {
                            let id = cache.components;
                            if let Err(x) = sm_save(&t, n, &scan, &out.join(format!("component_f{fi}_{id:04}"))) { err = Some(x); return; }
                            cache.components += 1;
                            for k in scan.states { cache.closed.insert(k, id); }
                        }
                    }
                    _ => { if let Err(x) = sm_save(&t, n, &scan, &out.join(format!("unresolved_{fi}_{}", count[1]))) { err = Some(x); return; } }
                }
                scan.status
            }
        };
        match status { "nonsmooth" => count[6] += 1, "closed_smooth" => count[7] += 1, _ => count[8] += 1 }
    });
    if let Some(x) = err { return Err(x); }
    let mut brute_note = String::from("-");
    if let Some(perms) = brute {
        let (mut eligible, mut replayed) = (0u128, 0u64);
        for &(c, ic, _) in perms {
            let e = p.compose(c, n).inverse(n);
            if !allowed.contains(&(ic, sm_index(e, n))) { continue; }
            let t = [a, b, c, e];
            match tri_t1_word(t, n, p) {
                None => {
                    eligible += 1;
                    if !enumerated.contains(&c) { return Err(format!("frame {fi}: T1-tri-compatible c missed")); }
                }
                Some(w) => {
                    if enumerated.contains(&c) { return Err(format!("frame {fi}: enumerated c fails T1-tri")); }
                    if replayed < 64 || c.img[0] == 0 {
                        if !tw_replay_fails(t, n, &w) { return Err(format!("frame {fi}: T1-tri word {w:?} does not fail")); }
                        replayed += 1;
                    }
                }
            }
        }
        if eligible != count[2] { return Err(format!("frame {fi}: brute {eligible} vs enumerated {}", count[2])); }
        brute_note = format!("ok:{eligible}:{replayed}");
    }
    let mut ws: Vec<_> = words.into_iter().collect();
    ws.sort();
    let mut wordlines = String::new();
    for (word, k) in ws { wordlines.push_str(&format!("{fi}\t{}\t{k}\n", word.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(","))); }
    let record = format!("{fi}\t{}\t{}\t{imprimitive}\t{brute_note}\t0\n", labs.len(), count.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"));
    Ok((count, record, wordlines))
}

// ---------------------------------------------------------------------------
// Constraint-first STAR sampler (--sample). NOT coverage: a sampled negative is
// not an exclusion. Each proposal draws a normal-form frame (ind b = largest
// index; a random with index s; b built from ind b cycle-merging transpositions so
// (a,b) is additive), then runs ArcSearch in RANDOM branch order with a node
// budget. So every emitted start already passes T1, T1c and T2 at j = +-1.
// Leaves are kept only with primitive <a,b,c> (exact Atkinson test), then full
// T2/T2c words, short words and closure. Closed and capped seeds are saved.
impl<'a> ArcSearch<'a> {
    fn go_random(&self, st: &mut ArcState, t: usize, rng: &mut SplitMix64, budget: &mut u64, leaves: &mut Vec<Perm>, want: usize) {
        if *budget == 0 || leaves.len() >= want { return; }
        *budget -= 1;
        let n = self.n;
        let walked = st.walked;
        let x = match self.next_point(st, t) {
            TwNext::Fail => { st.walked = walked; return; }
            TwNext::Leaf => {
                if st.open == [0, 0] && st.closed[0] == self.m && st.closed[1] >= self.tmin && st.closed[1] <= self.tmax && self.leaf_ok(st) {
                    let mut e = Perm::identity(n);
                    for i in 0..n { e.img[i] = st.e[i]; }
                    leaves.push(e);
                }
                st.walked = walked;
                return;
            }
            TwNext::Point(x) => x,
            TwNext::Value(y) => {
                let mut xs: Vec<u8> = (0..n as u8).filter(|&x| st.e[x as usize] == u8::MAX).collect();
                for i in (1..xs.len()).rev() { let j = rng.below(i + 1); xs.swap(i, j); }
                for x in xs {
                    if *budget == 0 || leaves.len() >= want { break; }
                    let mark = st.log.len(); let mm = st.mlog.len();
                    let (sc, so, py) = (st.closed, st.open, st.pre_done);
                    st.depth = t + 1; if self.try_arc(st, x, y) { self.go_random(st, t + 1, rng, budget, leaves, want); }
                    self.undo(st, mark); while st.mlog.len() > mm { let (li, old) = st.mlog.pop().unwrap(); st.lmask[li as usize] = old; }
                    st.closed = sc; st.open = so; st.pre_done = py;
                    st.used[y as usize] = false;
                    st.e[x as usize] = u8::MAX;
                }
                st.walked = walked;
                return;
            }
        };
        let mut ys: Vec<u8> = (0..n as u8).filter(|&y| !st.used[y as usize]).collect();
        for i in (1..ys.len()).rev() { let j = rng.below(i + 1); ys.swap(i, j); }
        for y in ys {
            if *budget == 0 || leaves.len() >= want { break; }
            let mark = st.log.len(); let mm = st.mlog.len();
            let (sc, so, py) = (st.closed, st.open, st.pre_done);
            st.depth = t + 1; if self.try_arc(st, x, y) { self.go_random(st, t + 1, rng, budget, leaves, want); }
            self.undo(st, mark); while st.mlog.len() > mm { let (li, old) = st.mlog.pop().unwrap(); st.lmask[li as usize] = old; }
            st.closed = sc; st.open = so; st.pre_done = py;
            st.used[y as usize] = false;
            st.e[x as usize] = u8::MAX;
        }
        st.walked = walked;
    }
}

// Uniform permutation of 0..n with exactly k cycles (Feller coupling, then relabel).
fn tws_with_cycles(n: usize, k: usize, st: &[Vec<f64>], rng: &mut SplitMix64) -> Perm {
    let mut fresh = vec![false; n];
    let mut kk = k;
    for m in (1..=n).rev() {
        let p_new = if kk == 0 { 0.0 } else { st[m - 1][kk - 1] / st[m][kk] };
        let u = (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        if u < p_new { fresh[m - 1] = true; kk -= 1; }
    }
    let mut p = Perm::identity(n);
    for i in 0..n { if !fresh[i] { let j = rng.below(i); p.img[i] = p.img[j]; p.img[j] = i as u8; } }
    let mut lab: Vec<u8> = (0..n as u8).collect();
    for i in (1..n).rev() { let j = rng.below(i + 1); lab.swap(i, j); }
    let mut q = Perm::identity(n);
    for i in 0..n { q.img[lab[i] as usize] = lab[p.img[i] as usize]; }
    q
}

fn tws_star_patterns(d: usize) -> Vec<[usize; 4]> {
    let mut out = Vec::new();
    for g in 1..d {
        let s = 2 * d + 2 * g - 2;
        for a in 2..d { for b in a + 1..d { for c in b..d {
            if s < a + b + c { continue; }
            let e = s - a - b - c;
            if c <= e && e < d && a + e <= d - 1 { out.push([a, b, c, e]); }
        }}}
    }
    out
}

pub fn twist_sampler_main(args: &[String]) -> Result<(), String> {
    let arg = |name: &str, default: &str| args.windows(2).find(|w| w[0] == name).map(|w| w[1].clone()).unwrap_or(default.into());
    let degrees: Vec<usize> = arg("--degrees", "16,17,18,19,20,21,22,23,24,25").split(',').map(|x| x.parse().map_err(|_| "degrees")).collect::<Result<_, _>>()?;
    let seed: u64 = arg("--seed", "1").parse().map_err(|_| "seed")?;
    let seconds: u64 = arg("--seconds", "31536000").parse().map_err(|_| "seconds")?;
    let budget0: u64 = arg("--node-budget", "20000").parse().map_err(|_| "budget")?;
    let want: usize = arg("--leaves", "8").parse().map_err(|_| "leaves")?;
    let cap: usize = arg("--cap", "300000").parse().map_err(|_| "cap")?;
    let closure_seconds: u64 = arg("--closure-seconds", "60").parse().map_err(|_| "closure")?;
    let any_group = args.iter().any(|a| a == "--any-group");
    let out = PathBuf::from(arg("--out", ""));
    if out.as_os_str().is_empty() { return Err("--out required".into()); }
    if degrees.iter().any(|&d| d < 9 || d > 25) { return Err("degrees 9..25".into()); }
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let mut stirling = vec![vec![0.0f64; 27]; 27];
    stirling[0][0] = 1.0;
    for n in 1..=26 { for k in 1..=n { stirling[n][k] = stirling[n - 1][k - 1] + (n - 1) as f64 * stirling[n - 1][k]; } }
    let pats: Vec<Vec<[usize; 4]>> = degrees.iter().map(|&d| tws_star_patterns(d)).collect();
    let mut rng = SplitMix64::new(seed.wrapping_mul(0x9E3779B97F4A7C15) ^ 0x7A157);
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    // frames, leaves, transitive, primitive, t2_reject, short_reject, scan_reject, closed, unresolved
    let mut stats = vec![[0u64; 9]; degrees.len()];
    let mut last = Instant::now();
    let mut seen_closed: HashSet<SmoothKey> = HashSet::new();
    let mut log = BufWriter::new(OpenOptions::new().create(true).append(true).open(out.join("PROGRESS.jsonl")).map_err(|e| e.to_string())?);
    let mut found = 0usize;
    while Instant::now() < deadline {
        let di = rng.below(degrees.len());
        let n = degrees[di];
        let q = pats[di][rng.below(pats[di].len())];
        let (ia, ib, ic, ie) = (q[0], q[3], q[1], q[2]);
        let a = tws_with_cycles(n, n - ia, &stirling, &mut rng);
        let mut cur = a;
        let mut b = Perm::identity(n);
        let mut ok = true;
        for _ in 0..ib {
            let lab = tw_labels(cur, n);
            let x = rng.below(n);
            let mut y = rng.below(n);
            let mut tries = 0;
            while lab[y] == lab[x] && tries < 1000 { y = rng.below(n); tries += 1; }
            if lab[y] == lab[x] { ok = false; break; }
            let mut t = Perm::identity(n); t.img[x] = y as u8; t.img[y] = x as u8;
            cur = cur.compose(t, n);
            b = b.compose(t, n);
        }
        if !ok || sm_index(b, n) != ib { continue; }
        stats[di][0] += 1;
        let p = a.compose(b, n);
        let pinv = p.inverse(n);
        let mut labs: Vec<[u8; TW_N]> = Vec::new();
        let mut ak = a;
        for _ in 0..tw_order(p, n) { let l = tw_labels(ak, n); if !labs.contains(&l) { labs.push(l); } ak = ak.conjugate_by(p, n); }
        let mut order: Vec<u8> = (0..n as u8).filter(|&x| a.img[x as usize] != x).collect();
        order.extend((0..n as u8).filter(|&x| a.img[x as usize] == x));
        let mut acyc: Vec<Vec<u8>> = Vec::new();
        let mut acyc_of = [u8::MAX; TW_N];
        let mut seen = [false; TW_N];
        for x in 0..n {
            if seen[x] || a.img[x] as usize == x { continue; }
            let mut cyc = Vec::new();
            let mut y = x;
            while !seen[y] { seen[y] = true; acyc_of[y] = acyc.len() as u8; cyc.push(y as u8); y = a.img[y] as usize; }
            acyc.push(cyc);
        }
        let lem = tw_shift_lemmas(tw_gated_lemmas(n, sm_index(b, n), ie, &[ic]), p, n);
        let srch = ArcSearch { n, m: n - ie, tmin: n - ic, tmax: n - ic, labs: &labs, pinv, a, ainv: a.inverse(n), order, acyc_of, acyc, blab: tw_labels(b, n),
            lt: [a, a.inverse(n), b, b.inverse(n), p, pinv], lem, sym: Vec::new(), lazy_depth: tw_lazy_depth(n) };
        let mut id = [0u8; TW_N];
        for i in 0..TW_N { id[i] = i as u8; }
        let mut st = ArcState { e: [u8::MAX; TW_N], used: [false; TW_N], einv: [u8::MAX; TW_N], ufs: vec![id; 2 * labs.len()],
            comp: [id; 2], hyp: [id; 3], walked: false, closed: [0; 2], open: [n; 2], img_done: [0; TW_N], pre_done: [0; TW_N], log: Vec::new(),
            kills: vec![0; srch.lem.len()], sym_skip: 0, nodes: 0, need: 0,
            lmask: srch.lem.iter().map(|l| if l.shift == 0 { u64::MAX } else { 0 }).collect(), mlog: Vec::new(), merged: [false; 2], depth: 0 };
        let mut budget = budget0;
        let mut leaves = Vec::new();
        srch.go_random(&mut st, 0, &mut rng, &mut budget, &mut leaves, want);
        for e in leaves {
            stats[di][1] += 1;
            let c = pinv.compose(e.inverse(n), n);
            if sm_index(c, n) != ic { continue; }
            let t = [a, b, c, e];
            if !generators_transitive(n, &t[..3]) { continue; }
            stats[di][2] += 1;
            if !any_group && !tw_primitive(&t[..3], n) { continue; }
            stats[di][3] += 1;
            let r = e.compose(a, n);
            let rinv = r.inverse(n);
            let (mut bj, mut cj) = (b, c);
            let mut word = None;
            for j in 0..tw_order(r, n) {
                if sm_index(a.compose(bj, n), n) != ia + ib { word = Some(vec![1usize; 2 * j]); break; }
                if sm_index(a.compose(cj, n), n) != ia + ic { let mut w = vec![1usize; 2 * j]; w.push(4); word = Some(w); break; }
                bj = rinv.compose(bj, n).compose(r, n);
                cj = rinv.compose(cj, n).compose(r, n);
            }
            if word.is_some() { stats[di][4] += 1; continue; }
            if sc_word(t, n).or_else(|| tw_short_dfs(t, n, 3)).is_some() { stats[di][5] += 1; continue; }
            let scan = sm_scan(&t, n, cap, Instant::now() + Duration::from_secs(closure_seconds));
            match scan.status {
                "nonsmooth" => stats[di][6] += 1,
                "closed_smooth" => {
                    stats[di][7] += 1;
                    let k = sm_key(&t, n);
                    if !seen_closed.contains(&k) {
                        sm_save(&t, n, &scan, &out.join(format!("CLOSED_d{n}_{found:04}")))?;
                        for s in &scan.states { seen_closed.insert(*s); }
                        found += 1;
                        writeln!(log, "{{\"event\":\"CLOSED\",\"degree\":{n},\"pattern\":{:?},\"states\":{},\"tuple\":\"{}\"}}", q, scan.states.len(), tuple_string_dynamic(n, &t)).map_err(|e| e.to_string())?;
                        log.flush().map_err(|e| e.to_string())?;
                    }
                }
                _ => {
                    stats[di][8] += 1;
                    sm_save(&t, n, &scan, &out.join(format!("UNRESOLVED_d{n}_{}", stats[di][8])))?;
                }
            }
        }
        if last.elapsed() > Duration::from_secs(60) {
            last = Instant::now();
            let rows: Vec<String> = degrees.iter().zip(&stats).map(|(d, s)| format!("\"{d}\":{:?}", s)).collect();
            writeln!(log, "{{\"event\":\"progress\",\"secs\":{:.0},\"cols\":\"frames,leaves,transitive,primitive,t2_reject,short_reject,scan_reject,closed,unresolved\",\"by_degree\":{{{}}}}}", start.elapsed().as_secs_f64(), rows.join(",")).map_err(|e| e.to_string())?;
            log.flush().map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}
