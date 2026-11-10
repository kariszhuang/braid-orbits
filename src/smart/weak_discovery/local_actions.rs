//! Independent local markings, not cosets of one common marked group.
use super::*;

pub(super) struct Audit {
    pub models: usize, pub combinations: usize, pub retained: usize,
    pub common_marked_groups: usize, pub complete: bool, pub seed: u64,
    pub block_max: usize, pub block_limit: usize, pub combination_cap: usize,
    pub frame_cap: usize,
}

// Ordered Cayley graph rooted at 1. This agrees with SmallPairGroup's
// marking_key, but also works for actual products/subdirect products.
fn marked_group_key(d: usize, a: Perm, b: Perm, cap: usize) -> Option<Vec<u16>> {
    let mut elements = vec![Perm::identity(d)];
    let mut labels = HashMap::from([(elements[0], 0usize)]);
    let mut key = Vec::new();
    let mut head = 0;
    while head < elements.len() {
        let p = elements[head]; head += 1;
        for g in [a, b] {
            let q = p.compose(g, d);
            let next = if let Some(&i) = labels.get(&q) { i } else {
                if elements.len() == cap { return None; }
                let i = elements.len(); elements.push(q); labels.insert(q, i); i
            };
            key.push(next as u16);
        }
    }
    Some(key)
}

type Catalog = (Vec<ActionModel>, HashSet<Vec<u16>>);
fn catalog(block_max: usize) -> Result<Catalog, String> {
    static CACHE: OnceLock<Mutex<HashMap<usize, Catalog>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(value) = cache.lock().unwrap().get(&block_max) { return Ok(value.clone()); }
    let mut models = Vec::new(); let mut common = HashSet::new(); let mut seen = HashSet::new();
    let mut add = |n: usize, a: Perm, b: Perm| {
        if n > block_max { return; }
        // Each catalog action is transitive, so tuple canonicalization is valid.
        let key = canonical(&[a, b, Perm::identity(n), a.compose(b, n).inverse(n)], n);
        if seen.insert((n, key)) { models.push((n, a, b)); }
    };
    let names = ["S3", "A4", "S4", "A5", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13"];
    for name in names {
        let group = SmallPairGroup::new(name)?;
        for marking in group.markings(false) {
            common.insert(group.marking_key(marking).unwrap().into_iter().map(u16::from).collect());
            for (n, a, b) in group.actions(marking) { add(n, a, b); }
        }
    }
    for (n, x, y) in fano_models() {
        for (a, b) in [(x, y), (x.compose(y, n), y), (x, x.compose(y, n))] {
            if let Some(key) = marked_group_key(n, a, b, 168) { common.insert(key); }
            add(n, a, b);
        }
    }
    models.sort_unstable_by_key(|(n, a, b)| (*n, a.img, b.img));
    let value = (models, common); cache.lock().unwrap().insert(block_max, value.clone()); Ok(value)
}

pub(super) fn frames(args: &[String], d: usize, ambient: Ambient, post_cut: bool) -> Result<(Vec<Frame>, Audit), String> {
    let block_max: usize = parse_num(get_arg(args, "--local-block-max", Some("8"))?, "--local-block-max")?;
    let block_limit: usize = parse_num(get_arg(args, "--point-orbit-max", Some("8"))?, "--point-orbit-max")?;
    let frame_cap: usize = parse_num(get_arg(args, "--local-frame-cap", Some("2048"))?, "--local-frame-cap")?;
    let combination_cap: usize = parse_num(get_arg(args, "--local-combination-cap", Some("50000"))?, "--local-combination-cap")?;
    let seed: u64 = parse_num(get_arg(args, "--seed", Some("20261004"))?, "--seed")?;
    if !(2..=MAX_D).contains(&block_max) || !(2..=MAX_D).contains(&block_limit)
        || !(1..=100000).contains(&frame_cap) || !(1..=10000000).contains(&combination_cap) {
        return Err("local action bounds: block-max 2..34, orbit-max 2..34, frames 1..100000, combinations 1..10000000".into());
    }
    let (mut models, common) = catalog(block_max.min(d))?;
    let mut rng = SplitMix64::new(seed);
    // A reproducible randomized prefix avoids always privileging the same
    // marking in capped runs. It does not make the resulting seeds uniform.
    for i in (1..models.len()).rev() { let j = rng.below(i + 1); models.swap(i, j); }
    struct Builder<'a> { d: usize, ambient: Ambient, post_cut: bool, models: &'a [ActionModel],
        common: &'a HashSet<Vec<u16>>, audit: Audit, out: Vec<Frame>, seen: HashSet<String> }
    impl Builder<'_> {
        fn visit(&mut self, begin: usize, used: usize, blocks: usize, a: Perm, b: Perm) {
            if !self.audit.complete { return; }
            if used == self.d {
                if self.audit.combinations == self.audit.combination_cap || self.out.len() == self.audit.frame_cap {
                    self.audit.complete = false; return;
                }
                self.audit.combinations += 1;
                if blocks < 2 || (!self.post_cut && (index(a,self.d)==0 || index(b,self.d)==0
                    || index(a,self.d)>self.d-3 || index(b,self.d)>self.d-3
                    || (self.ambient==Ambient::Alternating && (!a.is_even(self.d)||!b.is_even(self.d))))) { return; }
                if marked_group_key(self.d,a,b,168).is_some_and(|k|self.common.contains(&k)) { return; }
                let key = colored_representation_key(self.d, [a,b,Perm::identity(self.d)]);
                if !self.seen.insert(key) { return; }
                let h=a.compose(b,self.d);
                let left_genus=defect_of(&[a,b,Perm::identity(self.d),h.inverse(self.d)],self.d).1;
                if left_genus>self.d as i64-5 { return; }
                self.out.push(Frame {h,cycles:permutation_cycles(self.d,h),pair:Some((a,b)),left_genus}); return;
            }
            if blocks == self.audit.block_limit { return; }
            for i in begin..self.models.len() {
                let (n,x,y)=self.models[i]; if used+n>self.d {continue;}
                let mut aa=a; let mut bb=b;
                for j in 0..n {aa.img[used+j]=(used+x.img[j] as usize) as u8;bb.img[used+j]=(used+y.img[j] as usize) as u8;}
                self.visit(i,used+n,blocks+1,aa,bb); if !self.audit.complete {break;}
            }
        }
    }
    let audit=Audit {models:models.len(),combinations:0,retained:0,common_marked_groups:common.len(),complete:true,
        seed,block_max,block_limit,combination_cap,frame_cap};
    let mut builder=Builder {d,ambient,post_cut,models:&models,common:&common,audit,out:Vec::new(),seen:HashSet::new()};
    builder.visit(0,0,0,Perm::identity(d),Perm::identity(d));
    builder.audit.retained=builder.out.len(); Ok((builder.out,builder.audit))
}

impl Audit {
    pub(super) fn write(&self, out: &mut File) -> io::Result<()> {
        writeln!(out,"independent local action models = {}\nlocal combinations inspected = {}\nnovel marked pair frames retained = {}\ncommon marked groups excluded = {}\nlocal combination enumeration complete within configured library = {}\nlocal seed = {}\nlocal block maximum = {}\nlocal point-orbit maximum = {}\nlocal combination cap = {}\nlocal frame cap = {}",
            self.models,self.combinations,self.retained,self.common_marked_groups,self.complete,self.seed,self.block_max,self.block_limit,self.combination_cap,self.frame_cap)
    }
}

#[cfg(test)] mod tests {
    use super::*;
    #[test] fn actual_marked_keys_agree_with_abstract_keys() {
        for name in ["S3","A4","S4","A5","D5"] {
            let g=SmallPairGroup::new(name).unwrap();
            for pair in g.markings(false) {let actions=g.actions(pair);
                let &(n,a,b)=actions.iter().max_by_key(|x|x.0).unwrap();
                if subgroup_order(n,&[a,b])==g.table.len() as u128 {
                    assert_eq!(marked_group_key(n,a,b,168).unwrap(),g.marking_key(pair).unwrap().into_iter().map(u16::from).collect::<Vec<_>>());
                }
            }
        }
    }
    #[test] fn independent_blocks_can_escape_a_common_small_group() {
        let mut a=Perm::identity(9); let mut b=a;
        let left=pure_models("S3").unwrap().into_iter().find(|x|x.0==3).unwrap();
        let right=SmallPairGroup::new("A4").unwrap().actions(SmallPairGroup::new("A4").unwrap().base).into_iter().find(|x|x.0==6).unwrap();
        for (offset,(n,x,y)) in [(0,left),(3,right)] {for j in 0..n {a.img[offset+j]=(offset+x.img[j] as usize) as u8;b.img[offset+j]=(offset+y.img[j] as usize) as u8;}}
        let (_,common)=catalog(6).unwrap(); assert!(!common.contains(&marked_group_key(9,a,b,168).unwrap()));
        assert_eq!(point_orbits(9,&[a,b]).len(),2);
    }
}
