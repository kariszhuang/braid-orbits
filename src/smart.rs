//! Sparse, exact Hurwitz exploration; no dense |G|^2 tables.
use super::*;
use std::collections::{BTreeMap, BTreeSet};
#[path = "hybrid.rs"]
pub(super) mod hybrid;
mod defect_mission;
mod weak_discovery;
mod pair_fiber;
mod braid_frames;
mod positive_witness;
mod pure_stability;

// A Sims filter is only used to compress a generating set, never to guess its
// order. Each reduction multiplies by an already retained generator, so the
// generated subgroup is unchanged. At most d(d-1)/2 generators survive.
struct GeneratorFilter {
    inverse: [[u16; MAX_D]; MAX_D],
    inverses: Vec<Perm>,
    generators: Vec<Perm>,
}

impl GeneratorFilter {
    fn new(d: usize) -> Self {
        let capacity=d*(d-1)/2;
        Self { inverse: [[u16::MAX; MAX_D]; MAX_D], inverses: Vec::with_capacity(capacity),
               generators: Vec::with_capacity(capacity) }
    }
    fn insert(&mut self, mut g: Perm, d: usize, start: usize) {
        for i in start..d {
            let x = g.img[i] as usize;
            if x == i { continue; }
            let slot=self.inverse[i][x];
            if slot!=u16::MAX {
                g = self.inverses[slot as usize].compose(g, d);
            } else {
                self.inverse[i][x] = self.inverses.len() as u16;
                self.inverses.push(g.inverse(d));
                self.generators.push(g);
                return;
            }
        }
    }
}

pub(super) fn subgroup_order(d: usize, generators: &[Perm]) -> u128 {
    let id = Perm::identity(d);
    let mut filter = GeneratorFilter::new(d);
    for &g in generators { filter.insert(g, d, 0); }
    let mut current = filter.generators;
    let mut order = 1;
    for base in 0..d {
        if current.is_empty() { break; }
        let mut trans = [None; MAX_D];
        let mut orbit = Vec::with_capacity(d);
        trans[base] = Some(id);
        orbit.push(base);
        let mut head = 0;
        // Positive generators suffice: in a finite group their inverses are
        // positive powers. Schreier's lemma does not require a symmetric set.
        while head < orbit.len() {
            let x = orbit[head];
            head += 1;
            for &g in &current {
                let y = g.img[x] as usize;
                if trans[y].is_none() {
                    trans[y] = Some(g.compose(trans[x].unwrap(), d));
                    orbit.push(y);
                }
            }
        }
        order *= orbit.len() as u128;
        let inverses = trans.map(|x| x.map(|p| p.inverse(d)));
        let mut next = GeneratorFilter::new(d);
        for x in orbit {
            for &g in &current {
                let y = g.img[x] as usize;
                let h = inverses[y].unwrap().compose(g.compose(trans[x].unwrap(), d), d);
                next.insert(h, d, base + 1);
            }
        }
        current = next.generators;
    }
    order
}

// Three u128 words store the first three permutations.  d<=16 uses the
// historical 4-bit encoding; 17..25 uses 5 bits. Larger degrees use Lehmer
// permutation ranks (34! fits u128), preserving every existing certificate.
type Key = [u128; 3];

fn pair_transitive(d: usize, a: Perm, b: Perm) -> bool {
    generators_transitive(d,&[a,b])
}

fn has_transitive_pair(d: usize, t: &[Perm;4]) -> bool {
    (0..4).any(|i|(i+1..4).any(|j|pair_transitive(d,t[i],t[j])))
}

fn unpack(key: Key, d: usize) -> [Perm; 4] {
    let mut t = [Perm::identity(d); 4];
    let bits = if d <= 16 { 4 } else { 5 };
    let mask = (1u128 << bits) - 1;
    for j in 0..3 {
        if d<=25 {
            for i in 0..d { t[j].img[i] = ((key[j] >> (bits*i)) & mask) as u8; }
        } else {
            let mut rank=key[j];
            let mut digits=[0usize;MAX_D];
            for i in (0..d).rev() {digits[i]=(rank%(d-i) as u128) as usize; rank/=(d-i) as u128;}
            let mut available:Vec<_>=(0..d as u8).collect();
            for i in 0..d {t[j].img[i]=available.remove(digits[i]);}
            // Invalid out-of-range input decodes modulo d!; the verifier's
            // range/canonical checks reject it without panicking on input.
        }
    }
    t[3] = product_one_fourth(d, t[0], t[1], t[2]);
    t
}

// Canonical labeling of a TRANSITIVE tuple under all simultaneous S_d
// conjugations. For each possible root, traverse edges in generator order,
// label vertices by discovery order, then take the least encoded triple.
// This needs only d roots, not d! conjugating permutations. The fourth entry
// is determined by the first three. Relabeling commutes with every braid move
// and preserves both tested properties; quotient exhaustion is thus sound.
fn canonical(t: &[Perm; 4], d: usize) -> Key {
    let mut best = [u128::MAX; 3];
    let bits = if d <= 16 { 4 } else { 5 };
    for root in 0..d {
        let mut labels = [u8::MAX; MAX_D];
        let mut vertices = [0usize; MAX_D];
        labels[root] = 0;
        vertices[0] = root;
        let mut len = 1;
        let mut head = 0;
        while head < len {
            let x = vertices[head];
            head += 1;
            for g in &t[..3] {
                let y = g.img[x] as usize;
                if labels[y] == u8::MAX {
                    labels[y] = len as u8;
                    vertices[len] = y;
                    len += 1;
                }
            }
        }
        assert_eq!(len, d, "canonical labeling requires a transitive tuple");
        let mut key = [0; 3];
        for j in 0..3 {
            let mut available=(1u64<<d)-1;
            for i in 0..d {
                let value=labels[t[j].img[vertices[i]] as usize];
                if d<=25 { key[j] |= (value as u128) << (bits*i); }
                else {
                    let digit=(available & ((1u64<<value)-1)).count_ones();
                    key[j]=key[j]*(d-i) as u128+digit as u128;
                    available &= !(1u64<<value);
                }
            }
        }
        best = best.min(key);
    }
    best
}

#[derive(Debug)]
struct Exploration {
    status: &'static str,
    weak_witness: Option<[Perm; 4]>,
    weak_word: Option<Vec<u8>>,
    states: Vec<Key>,
    orders: BTreeMap<u128, usize>,
}

fn explore(d: usize, seed: [Perm; 4], cap: usize, deadline: Instant) -> Exploration {
    explore_goal(d, seed, cap, deadline, false)
}

fn explore_goal(d: usize, seed: [Perm; 4], cap: usize, deadline: Instant, weak_only: bool) -> Exploration {
    let root = canonical(&seed, d);
    let mut out = Exploration { status: "capped", weak_witness: None, weak_word: None, states: vec![root], orders: BTreeMap::new() };
    let mut parents = if weak_only { vec![(0usize, 0u8)] } else { Vec::new() };
    let mut seen = HashSet::from([root]);
    let mut head = 0;
    while head < out.states.len() {
        if Instant::now() >= deadline { out.status = "timeout"; return out; }
        let t = unpack(out.states[head], d);
        head += 1;
        let connected = pair_transitive(d, t[0], t[1]);
        if connected && out.weak_witness.is_none() { out.weak_witness = Some(t); }
        if connected && weak_only {
            let mut word = Vec::new(); let mut cursor = head - 1;
            while cursor != 0 { let (parent, mv) = parents[cursor]; word.push(mv); cursor = parent; }
            word.reverse(); out.weak_word = Some(word);
            out.status = "weak_found"; return out;
        }
        // Intransitive pairs cannot contain A_d. Count their orders too only
        // for an exhausted certificate, after traversal, to keep search cheap.
        if connected {
            let order = subgroup_order(d, &t[..2]);
            if order >= factorial_u128(d)/2 { out.status = "strong_found"; return out; }
            *out.orders.entry(order).or_default() += 1;
        }
        // Each braid generator is a permutation of this finite state space.
        // A finite set closed under the three positive generators is therefore
        // closed under their inverses too. The verifier still checks all six.
        for mv in 0..3 {
            let key = canonical(&hurwitz_move_dynamic(d, t, mv), d);
            if seen.contains(&key) { continue; }
            if seen.len() >= cap { return out; }
            seen.insert(key);
            out.states.push(key);
            if weak_only { parents.push((head - 1, mv as u8)); }
        }
    }
    out.status = "certified_strong_failure";
    if out.weak_witness.is_none() { out.status = "certified_weak_and_strong_failure"; }
    for &key in &out.states {
        let t=unpack(key,d);
        if !pair_transitive(d,t[0],t[1]) {
            *out.orders.entry(subgroup_order(d,&t[..2])).or_default() += 1;
        }
    }
    out
}

fn parse_tuple(text: &str, d: usize) -> Result<[Perm; 4], String> {
    let text = text.trim().strip_prefix('[').and_then(|s| s.strip_suffix(']')).ok_or("expected [four permutations]")?;
    let parts: Vec<_> = text.split(',').collect();
    if parts.len() != 4 { return Err("expected exactly four permutations".into()); }
    let mut t = [Perm::identity(d); 4];
    for (i, part) in parts.iter().enumerate() {
        let mut rest = part.trim();
        if rest.is_empty() { return Err("empty permutation; use () for identity".into()); }
        let mut used = HashSet::new();
        while !rest.is_empty() {
            rest = rest.strip_prefix('(').ok_or("expected cycle opening parenthesis")?;
            let (cycle, after) = rest.split_once(')').ok_or("unclosed cycle")?;
            let points: Vec<usize> = cycle.split_whitespace().map(|x| x.parse::<usize>().map_err(|_| "invalid cycle point".to_string())).collect::<Result<_, _>>()?;
            for &x in &points {
                if x == 0 || x > d || !used.insert(x) { return Err("cycle points must be distinct and in 1..d".into()); }
            }
            for j in 0..points.len() { t[i].img[points[j]-1] = (points[(j+1)%points.len()]-1) as u8; }
            rest = after.trim();
        }
    }
    if !tuple_product_is_identity(d, &t) { return Err("tuple product is not identity (right-to-left composition)".into()); }
    Ok(t)
}

fn validate_tuple(d: usize, ambient: Ambient, t: &[Perm;4]) -> Result<(), String> {
    if !tuple_product_is_identity(d, t) { return Err("tuple product is not identity".into()); }
    if ambient == Ambient::Alternating && t.iter().any(|p| !p.is_even(d)) { return Err("odd entry in A_d input".into()); }
    let order = subgroup_order(d, &t[..3]);
    if order != ambient.target_order(d) { return Err(format!("tuple generates order {order}, expected {}", ambient.target_order(d))); }
    Ok(())
}

fn write_certificate(path: &Path, d: usize, ambient: Ambient, seed: &[Perm;4], result: &Exploration) -> io::Result<()> {
    let mut w = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?);
    writeln!(w, "HURWITZ_SPARSE_CERTIFICATE_V1\ngroup={}\nd={}\nseed={}\nstatus={}\nquotient=S_d simultaneous conjugation\nstates={}", ambient.short(), d, tuple_string_dynamic(d, seed), result.status, result.states.len())?;
    if d>25 {writeln!(w,"encoding=Lehmer permutation ranks")?;}
    let representative=result.states.iter().min().ok_or_else(||io::Error::new(io::ErrorKind::InvalidData,"empty certificate"))?;
    writeln!(w,"representative={:032x}/{:032x}/{:032x}",representative[0],representative[1],representative[2])?;
    if let Some(t) = result.weak_witness { writeln!(w, "weak_witness={}", tuple_string_dynamic(d, &t))?; }
    for (order, count) in &result.orders { writeln!(w, "ORDER\t{order}\t{count}")?; }
    for key in &result.states { writeln!(w, "STATE\t{:032x}\t{:032x}\t{:032x}", key[0], key[1], key[2])?; }
    w.flush()?;
    w.get_ref().sync_all()
}

fn verify(path: &Path) -> Result<(), String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let value = |name: &str| -> Result<&str,String> {
        text.lines().find_map(|s| s.strip_prefix(name)).ok_or_else(||format!("missing {name}"))
    };
    if text.lines().next() != Some("HURWITZ_SPARSE_CERTIFICATE_V1") { return Err("unsupported certificate".into()); }
    let d: usize = value("d=")?.parse().map_err(|_| "invalid degree")?;
    if !(4..=MAX_D).contains(&d) { return Err("invalid degree".into()); }
    let ambient = Ambient::parse(value("group=")?)?;
    let seed = parse_tuple(value("seed=")?,d)?;
    validate_tuple(d,ambient,&seed)?;
    let mut states = HashSet::new();
    for line in text.lines().filter_map(|s| s.strip_prefix("STATE\t")) {
        let fields:Vec<_>=line.split('\t').collect();
        if fields.len()!=3 { return Err("invalid state".into()); }
        let mut key=[0;3];
        for i in 0..3 { key[i]=u128::from_str_radix(fields[i],16).map_err(|_|"invalid state encoding")?; }
        let t=unpack(key,d);
        for (i,p) in t[..3].iter().enumerate() {
            let mut images=p.img[..d].to_vec(); images.sort_unstable();
            if images != (0..d as u8).collect::<Vec<_>>() || (d<16 && key[i]>>(4*d)!=0) { return Err("state is not a permutation".into()); }
        }
        if !states.insert(key) { return Err("duplicate certificate state".into()); }
    }
    let expected:usize=value("states=")?.parse().map_err(|_|"invalid state count")?;
    if states.len()!=expected || !states.contains(&canonical(&seed,d)) { return Err("missing states or seed".into()); }
    if let Ok(recorded)=value("representative=") {
        if parse_hex_key(&recorded.split('/').collect::<Vec<_>>())?!=*states.iter().min().unwrap() {
            return Err("incorrect component representative".into());
        }
    }
    let mut weak=false;
    let mut histogram=BTreeMap::new();
    for &key in &states {
        let t=unpack(key,d);
        let bits = if d <= 16 { 4 } else { 5 };
        // Check transitivity first, so invalid input never reaches canonical's
        // transitive-only precondition. This also catches malformed vertices.
        validate_tuple(d,ambient,&t)?;
        if canonical(&t,d)!=key || key.iter().any(|word| if d<=25 {(word >> (bits*d)) != 0} else {*word>=factorial_u128(d)}) { return Err("noncanonical state".into()); }
        let order=subgroup_order(d,&t[..2]);
        if order>=factorial_u128(d)/2 { return Err("certificate has a strong pair".into()); }
        *histogram.entry(order).or_insert(0usize)+=1;
        weak |= point_orbit_count_perms(d,t[0],t[1])==1;
        for mv in 0..6 {
            if !states.contains(&canonical(&hurwitz_move_dynamic(d,t,mv),d)) { return Err("certificate set is not closed under braid moves".into()); }
        }
    }
    let mut recorded=BTreeMap::new();
    for line in text.lines().filter_map(|s|s.strip_prefix("ORDER\t")) {
        let (a,b)=line.split_once('\t').ok_or("invalid order histogram")?;
        let a:u128=a.parse().map_err(|_|"invalid order")?;
        let b:usize=b.parse().map_err(|_|"invalid order count")?;
        if recorded.insert(a,b).is_some() { return Err("duplicate order".into()); }
    }
    if recorded!=histogram { return Err("incorrect order histogram".into()); }
    let root=canonical(&seed,d);
    let mut reached=HashSet::from([root]);
    let mut queue=vec![root];
    let mut head=0;
    while head<queue.len() {
        let t=unpack(queue[head],d); head+=1;
        for mv in 0..3 {
            let key=canonical(&hurwitz_move_dynamic(d,t,mv),d);
            if reached.insert(key) { queue.push(key); }
        }
    }
    if reached!=states { return Err("certificate is not one connected quotient component".into()); }
    if weak {
        let witness=parse_tuple(value("weak_witness=")?,d)?;
        validate_tuple(d,ambient,&witness)?;
        if !states.contains(&canonical(&witness,d)) || point_orbit_count_perms(d,witness[0],witness[1])!=1 {
            return Err("invalid weak witness".into());
        }
    }
    let expected_status=if weak {"certified_strong_failure"} else {"certified_weak_and_strong_failure"};
    if value("status=")? != expected_status { return Err("incorrect certificate status".into()); }
    println!("VERIFIED {}_{}: {} quotient states, strong failure, weak={weak}",ambient.short(),d,states.len());
    Ok(())
}

fn index(p: Perm, d: usize) -> usize {
    let mut seen=[false;MAX_D];
    let mut cycles=0;
    for i in 0..d {
        if seen[i] {continue;}
        cycles+=1;
        let mut j=i;
        while !seen[j] {seen[j]=true; j=p.img[j] as usize;}
    }
    d-cycles
}

// ---------- Quotient Cayley--Schreier spectra ----------

#[derive(Clone)]
struct SpectrumCase {
    id: String,
    kind: String,
    source: String,
    ambient: Ambient,
    d: usize,
    seed: [Perm; 4],
    states: Option<Vec<Key>>,
    input_states: Option<usize>,
    aliases: Vec<String>,
}

#[derive(Clone)]
struct SpectrumGraph {
    states: Vec<Key>,
    transitions: Vec<[u32; 6]>,
}

#[derive(Clone)]
struct EigenResult {
    method: String,
    values: Vec<f64>,
    lambda_max: f64,
    lambda2: f64,
    lambda_min: f64,
}

#[derive(Clone)]
struct BaselineResult {
    kind: String,
    replicate: usize,
    components: usize,
    eigen: EigenResult,
    diagnostics: GraphDiagnostics,
}

#[derive(Clone)]
struct ComputedSpectrum {
    case: SpectrumCase,
    graph: SpectrumGraph,
    vertices: usize,
    eigen: EigenResult,
    components: usize,
    loops: usize,
    parallel_excess: usize,
    diagnostics: GraphDiagnostics,
    convergence: Vec<(usize, f64, f64)>,
    baselines: Vec<BaselineResult>,
}

#[derive(Clone, Default)]
struct GraphDiagnostics {
    generator_profiles: [String; 3],
    generator_fixed: [usize; 3],
    generator_cycles: [usize; 3],
    braid_12_defects: usize,
    braid_23_defects: usize,
    commute_13_defects: usize,
    bipartite: bool,
    maxcut_conflict_fraction: f64,
    wl_hash: u64,
    canonical_labeled_code: String,
    moments: Vec<(usize, f64, f64, String)>,
}

struct SpectrumOptions {
    baseline_reps: usize,
    baseline_kinds: HashSet<String>,
    trace_order: usize,
    trace_probes: usize,
    maxcut_restarts: usize,
    canonical_max: usize,
    switch_attempts_per_edge: usize,
    lanczos_checks: Vec<usize>,
}

fn optional_arg(args: &[String], name: &str) -> Option<String> {
    args.iter().position(|x| x == name).and_then(|i| args.get(i + 1)).cloned()
}

fn csv_fields(line: &str) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let ch = chars[i];
        if quoted {
            if ch == '"' {
                if i + 1 < chars.len() && chars[i + 1] == '"' {
                    field.push('"');
                    i += 1;
                } else {
                    quoted = false;
                }
            } else {
                field.push(ch);
            }
        } else {
            match ch {
                '"' if field.trim().is_empty() => quoted = true,
                ',' => {
                    out.push(field.trim().to_string());
                    field.clear();
                }
                _ => field.push(ch),
            }
        }
        i += 1;
    }
    if quoted { return Err("unterminated CSV quote".into()); }
    out.push(field.trim().to_string());
    Ok(out)
}

fn infer_group_degree(path: &Path) -> Result<(Ambient, usize), String> {
    for parent in path.ancestors().take(8) {
        let Some(name) = parent.file_name().and_then(|x| x.to_str()) else { continue; };
        let group_name = name.split('_').next().unwrap_or(name);
        let bytes = group_name.as_bytes();
        if bytes.len() < 2 || (bytes[0] != b'A' && bytes[0] != b'S') || !bytes[1..].iter().all(u8::is_ascii_digit) {
            continue;
        }
        let d: usize = group_name[1..].parse().map_err(|_| format!("invalid degree in {name}"))?;
        if (4..=MAX_D).contains(&d) {
            return Ok((if bytes[0] == b'A' { Ambient::Alternating } else { Ambient::Symmetric }, d));
        }
    }
    Err(format!("cannot infer A_d/S_d and d from component path {}", path.display()))
}

fn parse_components_file(path: &Path) -> Result<Vec<SpectrumCase>, String> {
    let (ambient, d) = infer_group_degree(path)?;
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut lines = text.lines();
    let header = lines.next().ok_or_else(|| format!("empty component file {}", path.display()))?;
    let names = csv_fields(header)?;
    let seed_col = names.iter().position(|x| x == "seed").ok_or("component CSV has no seed column")?;
    let component_col = names.iter().position(|x| x == "component");
    let states_col = names.iter().position(|x| x == "quotient_states");
    let mut out = Vec::new();
    for (line_no, line) in lines.enumerate() {
        if line.trim().is_empty() { continue; }
        let fields = csv_fields(line).map_err(|e| format!("{}:{}: {e}", path.display(), line_no + 2))?;
        if fields.len() != names.len() { return Err(format!("{}:{}: CSV field count mismatch", path.display(), line_no + 2)); }
        let seed = parse_tuple(&fields[seed_col], d).map_err(|e| format!("{}:{}: {e}", path.display(), line_no + 2))?;
        validate_tuple(d, ambient, &seed).map_err(|e| format!("{}:{}: {e}", path.display(), line_no + 2))?;
        let component = component_col.map(|i| fields[i].clone()).unwrap_or_else(|| (line_no + 1).to_string());
        let input_states = states_col.map(|i| fields[i].parse::<usize>().map_err(|_| format!("{}:{}: invalid quotient_states", path.display(), line_no + 2))).transpose()?;
        let shard = path.parent().and_then(Path::file_name).and_then(|x| x.to_str()).filter(|x| x.starts_with("shard_"));
        let scope = shard.map(|x| format!(":{x}")).unwrap_or_default();
        out.push(SpectrumCase {
            id: format!("{}{}{}:component-{}", ambient.short(), d, scope, component),
            kind: "exact_component".into(),
            source: path.display().to_string(),
            ambient,
            d,
            seed,
            states: None,
            input_states,
            aliases: vec![component],
        });
    }
    Ok(out)
}

fn parse_hex_key(fields: &[&str]) -> Result<Key, String> {
    if fields.len() != 3 { return Err("invalid state encoding".into()); }
    Ok([
        u128::from_str_radix(fields[0], 16).map_err(|_| "invalid state encoding")?,
        u128::from_str_radix(fields[1], 16).map_err(|_| "invalid state encoding")?,
        u128::from_str_radix(fields[2], 16).map_err(|_| "invalid state encoding")?,
    ])
}

fn parse_certificate_case(path: &Path) -> Result<SpectrumCase, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let value = |name: &str| -> Result<&str, String> {
        text.lines().find_map(|s| s.strip_prefix(name)).ok_or_else(|| format!("{}: missing {name}", path.display()))
    };
    if text.lines().next() != Some("HURWITZ_SPARSE_CERTIFICATE_V1") {
        return Err(format!("{}: unsupported certificate", path.display()));
    }
    let d: usize = value("d=")?.parse().map_err(|_| format!("{}: invalid degree", path.display()))?;
    let ambient = Ambient::parse(value("group=")?)?;
    let seed = parse_tuple(value("seed=")?, d)?;
    validate_tuple(d, ambient, &seed)?;
    let mut states = Vec::new();
    for line in text.lines().filter_map(|s| s.strip_prefix("STATE\t")) {
        states.push(parse_hex_key(&line.split('\t').collect::<Vec<_>>())?);
    }
    let expected: usize = value("states=")?.parse().map_err(|_| "invalid state count")?;
    if states.len() != expected || states.is_empty() { return Err(format!("{}: state count mismatch", path.display())); }
    Ok(SpectrumCase {
        id: format!("{}:{}", ambient.short(), path.file_stem().and_then(|x| x.to_str()).unwrap_or("certificate")),
        kind: "certificate".into(),
        source: path.display().to_string(),
        ambient,
        d,
        seed,
        states: Some(states),
        input_states: Some(expected),
        aliases: Vec::new(),
    })
}

fn collect_certificate_paths(path: &Path, out: &mut Vec<PathBuf>) -> Result<(), String> {
    let meta = fs::metadata(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if meta.is_file() {
        if path.file_name().and_then(|x| x.to_str()).unwrap_or("").starts_with("CERTIFIED_") { out.push(path.to_path_buf()); }
        return Ok(());
    }
    let mut entries: Vec<_> = fs::read_dir(path).map_err(|e| format!("{}: {e}", path.display()))?.collect::<Result<Vec<_>, _>>().map_err(|e| e.to_string())?;
    entries.sort_by_key(|e| e.path());
    for entry in entries { collect_certificate_paths(&entry.path(), out)?; }
    Ok(())
}

fn add_spectrum_case(cases: &mut Vec<SpectrumCase>, roots: &mut HashMap<(Ambient, usize, Key), usize>, mut case: SpectrumCase) {
    let root = canonical(&case.seed, case.d);
    let key = (case.ambient, case.d, root);
    if let Some(&i) = roots.get(&key) {
        cases[i].aliases.extend(case.aliases.drain(..));
        return;
    }
    roots.insert(key, cases.len());
    cases.push(case);
}

fn graph_from_states(d: usize, mut states: Vec<Key>) -> Result<SpectrumGraph, String> {
    states.sort_unstable();
    states.dedup();
    if states.is_empty() { return Err("empty quotient graph".into()); }
    if states.len() > u32::MAX as usize { return Err("quotient graph exceeds u32 vertex capacity".into()); }
    let mut index = HashMap::with_capacity(states.len() * 2);
    for (i, &key) in states.iter().enumerate() { index.insert(key, i as u32); }
    let mut transitions = vec![[0u32; 6]; states.len()];
    for (i, &key) in states.iter().enumerate() {
        let t = unpack(key, d);
        for mv in 0..6 {
            let next = canonical(&hurwitz_move_dynamic(d, t, mv), d);
            transitions[i][mv] = *index.get(&next).ok_or_else(|| format!("missing braid neighbor for state {i}, move {mv}"))?;
        }
    }
    Ok(SpectrumGraph { states, transitions })
}

fn enumerate_spectrum_graph(d: usize, seed: [Perm; 4]) -> Result<SpectrumGraph, String> {
    let root = canonical(&seed, d);
    let mut states = vec![root];
    let mut seen = HashSet::from([root]);
    let mut head = 0;
    while head < states.len() {
        let t = unpack(states[head], d);
        head += 1;
        for mv in 0..6 {
            let next = canonical(&hurwitz_move_dynamic(d, t, mv), d);
            if seen.insert(next) { states.push(next); }
        }
    }
    drop(seen);
    graph_from_states(d, states)
}

fn graph_components(graph: &SpectrumGraph) -> usize {
    let n = graph.transitions.len();
    let mut seen = vec![false; n];
    let mut count = 0;
    for start in 0..n {
        if seen[start] { continue; }
        count += 1;
        let mut q = vec![start];
        seen[start] = true;
        let mut head = 0;
        while head < q.len() {
            let i = q[head];
            head += 1;
            for &j in &graph.transitions[i] {
                let j = j as usize;
                if !seen[j] { seen[j] = true; q.push(j); }
            }
        }
    }
    count
}

fn graph_shape(graph: &SpectrumGraph) -> (usize, usize) {
    let mut loops = 0;
    let mut parallel_excess = 0;
    for (i, row) in graph.transitions.iter().enumerate() {
        loops += row.iter().filter(|&&j| j as usize == i).count();
        let mut counts = HashMap::<u32, usize>::new();
        for &j in row { *counts.entry(j).or_default() += 1; }
        parallel_excess += counts.values().map(|&n| n.saturating_sub(1)).sum::<usize>();
    }
    (loops, parallel_excess)
}

fn permutation_profile(p: &[usize]) -> (String, usize, usize) {
    let mut image = vec![false; p.len()];
    if p.iter().any(|&value| value >= p.len() || std::mem::replace(&mut image[value], true)) {
        return ("not-permutation".into(), p.iter().enumerate().filter(|&(i, value)| i == *value).count(), 0);
    }
    let mut seen = vec![false; p.len()];
    let mut counts = BTreeMap::<usize, usize>::new();
    for start in 0..p.len() {
        if seen[start] { continue; }
        let (mut i, mut length) = (start, 0usize);
        while !seen[i] { seen[i] = true; length += 1; i = p[i]; }
        *counts.entry(length).or_default() += 1;
    }
    let fixed = counts.get(&1).copied().unwrap_or(0);
    let cycles = counts.values().sum();
    let profile = counts.into_iter().map(|(length, count)| format!("{length}^{count}")).collect::<Vec<_>>().join(",");
    (profile, fixed, cycles)
}

fn positive_generators(graph: &SpectrumGraph) -> [Vec<usize>; 3] {
    std::array::from_fn(|slot| graph.transitions.iter().map(|row| row[slot] as usize).collect())
}

fn relation_defects(generators: &[Vec<usize>; 3]) -> (usize, usize, usize) {
    let apply = |start: usize, word: &[usize]| word.iter().fold(start, |i, &g| generators[g][i]);
    let mut braid_12 = 0;
    let mut braid_23 = 0;
    let mut commute_13 = 0;
    for i in 0..generators[0].len() {
        braid_12 += usize::from(apply(i, &[0, 1, 0]) != apply(i, &[1, 0, 1]));
        braid_23 += usize::from(apply(i, &[1, 2, 1]) != apply(i, &[2, 1, 2]));
        commute_13 += usize::from(apply(i, &[0, 2]) != apply(i, &[2, 0]));
    }
    (braid_12, braid_23, commute_13)
}

fn graph_is_bipartite(graph: &SpectrumGraph) -> bool {
    let mut color = vec![2u8; graph.transitions.len()];
    for start in 0..color.len() {
        if color[start] != 2 { continue; }
        color[start] = 0;
        let mut queue = vec![start];
        let mut head = 0;
        while head < queue.len() {
            let i = queue[head]; head += 1;
            for &raw in &graph.transitions[i] {
                let j = raw as usize;
                if i == j { return false; }
                if color[j] == 2 { color[j] = 1 - color[i]; queue.push(j); }
                else if color[j] == color[i] { return false; }
            }
        }
    }
    true
}

fn maxcut_conflict_fraction(graph: &SpectrumGraph, restarts: usize, seed: u64) -> f64 {
    if graph.transitions.is_empty() || restarts == 0 { return f64::NAN; }
    let mut rng = SplitMix64::new(seed);
    let mut best = usize::MAX;
    for restart in 0..restarts {
        let mut color = (0..graph.transitions.len()).map(|i| if restart == 0 { (i & 1) as u8 } else { rng.below(2) as u8 }).collect::<Vec<_>>();
        for _ in 0..32 {
            let mut changed = false;
            for i in 0..color.len() {
                let mut same = 0;
                let mut different = 0;
                for &raw in &graph.transitions[i] {
                    let j = raw as usize;
                    if i == j { continue; }
                    if color[i] == color[j] { same += 1; } else { different += 1; }
                }
                if same > different { color[i] ^= 1; changed = true; }
            }
            if !changed { break; }
        }
        let conflicts = graph.transitions.iter().enumerate().map(|(i, row)| row.iter().filter(|&&j| color[i] == color[j as usize]).count()).sum();
        best = best.min(conflicts);
    }
    best as f64 / (6 * graph.transitions.len()) as f64
}

fn fnv_mix(mut hash: u64, value: u64) -> u64 {
    for byte in value.to_le_bytes() { hash = (hash ^ byte as u64).wrapping_mul(0x100000001b3); }
    hash
}

fn wl_graph_hash(graph: &SpectrumGraph) -> u64 {
    if graph.transitions.is_empty() { return 0; }
    let mut signatures = Vec::with_capacity(graph.transitions.len());
    for (i, row) in graph.transitions.iter().enumerate() {
        let mut signature = Vec::with_capacity(12);
        for &j in row { signature.push(usize::from(j as usize == i)); }
        for a in 0..6 { for b in 0..a { signature.push(usize::from(row[a] == row[b])); } }
        signatures.push(signature);
    }
    let mut unique = signatures.clone(); unique.sort(); unique.dedup();
    let mut colors = signatures.iter().map(|x| unique.binary_search(x).unwrap()).collect::<Vec<_>>();
    for _ in 0..32 {
        let next_signatures = graph.transitions.iter().enumerate().map(|(i, row)| {
            let mut signature = Vec::with_capacity(7);
            signature.push(colors[i]);
            signature.extend(row.iter().map(|&j| colors[j as usize]));
            signature
        }).collect::<Vec<_>>();
        let mut unique = next_signatures.clone(); unique.sort(); unique.dedup();
        let next = next_signatures.iter().map(|x| unique.binary_search(x).unwrap()).collect::<Vec<_>>();
        if next == colors { break; }
        colors = next;
    }
    let mut final_signatures = graph.transitions.iter().enumerate().map(|(i, row)| {
        let mut signature = Vec::with_capacity(7);
        signature.push(colors[i]);
        signature.extend(row.iter().map(|&j| colors[j as usize]));
        signature
    }).collect::<Vec<_>>();
    final_signatures.sort();
    let mut hash = fnv_mix(0xcbf29ce484222325, graph.transitions.len() as u64);
    for signature in final_signatures { for value in signature { hash = fnv_mix(hash, value as u64); } }
    hash
}

// For a connected graph with fixed move labels, an isomorphism is determined
// by the image of one vertex. BFS from every possible image therefore gives
// a collision-free (unlike WL) canonical code under vertex relabeling.
fn canonical_labeled_code(graph: &SpectrumGraph, maximum: usize) -> String {
    let n = graph.transitions.len();
    if n == 0 || n > maximum { return String::new(); }
    let mut best: Option<Vec<u32>> = None;
    for root in 0..n {
        let mut old_to_new = vec![usize::MAX; n];
        let mut queue = Vec::with_capacity(n);
        old_to_new[root] = 0;
        queue.push(root);
        let mut code = Vec::with_capacity(6 * n);
        let mut head = 0;
        while head < queue.len() {
            let vertex = queue[head];
            head += 1;
            for &next in &graph.transitions[vertex] {
                let next = next as usize;
                if old_to_new[next] == usize::MAX {
                    old_to_new[next] = queue.len();
                    queue.push(next);
                }
                code.push(old_to_new[next] as u32);
            }
        }
        if queue.len() != n { return String::new(); }
        if best.as_ref().is_none_or(|previous| code < *previous) { best = Some(code); }
    }
    best.unwrap().iter().map(u32::to_string).collect::<Vec<_>>().join(",")
}

fn stochastic_trace_moments(graph: &SpectrumGraph, eigen: &EigenResult, max_order: usize, probes: usize, seed: u64) -> Vec<(usize, f64, f64, String)> {
    if max_order == 0 || graph.transitions.is_empty() { return Vec::new(); }
    if !eigen.values.is_empty() {
        return (1..=max_order).map(|k| {
            let estimate = eigen.values.iter().map(|x| x.powi(k as i32)).sum::<f64>() / (graph.transitions.len() as f64 * 6f64.powi(k as i32));
            (k, estimate, 0.0, "eigen-exact".into())
        }).collect();
    }
    let probes = probes.max(1);
    let mut rng = SplitMix64::new(seed);
    let mut samples = vec![Vec::<f64>::with_capacity(probes); max_order];
    for _ in 0..probes {
        let z = (0..graph.transitions.len()).map(|_| if rng.below(2) == 0 { -1.0 } else { 1.0 }).collect::<Vec<_>>();
        let mut v = z.clone();
        for k in 1..=max_order {
            v = graph_matvec(graph, &v);
            samples[k - 1].push(dot(&z, &v) / (graph.transitions.len() as f64 * 6f64.powi(k as i32)));
        }
    }
    samples.into_iter().enumerate().map(|(i, values)| {
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance = if values.len() > 1 { values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64 } else { 0.0 };
        (i + 1, mean, (variance / values.len() as f64).sqrt(), format!("hutchinson-{probes}"))
    }).collect()
}

fn graph_diagnostics(graph: &SpectrumGraph, eigen: &EigenResult, trace_order: usize, trace_probes: usize, maxcut_restarts: usize, canonical_max: usize, seed: u64) -> GraphDiagnostics {
    let generators = positive_generators(graph);
    let profiles = generators.iter().map(|p| permutation_profile(p)).collect::<Vec<_>>();
    let (braid_12_defects, braid_23_defects, commute_13_defects) = relation_defects(&generators);
    GraphDiagnostics {
        generator_profiles: std::array::from_fn(|i| profiles[i].0.clone()),
        generator_fixed: std::array::from_fn(|i| profiles[i].1),
        generator_cycles: std::array::from_fn(|i| profiles[i].2),
        braid_12_defects, braid_23_defects, commute_13_defects,
        bipartite: graph_is_bipartite(graph),
        maxcut_conflict_fraction: maxcut_conflict_fraction(graph, maxcut_restarts, seed ^ 0xb1b1),
        wl_hash: if graph.transitions.len() <= 4096 { wl_graph_hash(graph) } else { 0 },
        canonical_labeled_code: canonical_labeled_code(graph, canonical_max),
        moments: stochastic_trace_moments(graph, eigen, trace_order, trace_probes, seed ^ 0x7ace),
    }
}

fn jacobi_eigenvalues(mut a: Vec<f64>, n: usize) -> Vec<f64> {
    if n == 0 { return Vec::new(); }
    let tol = 1e-11;
    let max_sweeps = 80 * n.max(1);
    for _ in 0..max_sweeps {
        let mut largest: f64 = 0.0;
        for p in 0..n.saturating_sub(1) {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                largest = largest.max(apq.abs());
                if apq.abs() <= tol { continue; }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                let tau = (aqq - app) / (2.0 * apq);
                let t = if tau >= 0.0 { 1.0 / (tau + (1.0 + tau * tau).sqrt()) } else { -1.0 / (-tau + (1.0 + tau * tau).sqrt()) };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                for k in 0..n {
                    if k == p || k == q { continue; }
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    let new_kp = c * akp - s * akq;
                    let new_kq = s * akp + c * akq;
                    a[k * n + p] = new_kp;
                    a[p * n + k] = new_kp;
                    a[k * n + q] = new_kq;
                    a[q * n + k] = new_kq;
                }
                a[p * n + p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                a[q * n + q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;
            }
        }
        if largest < tol { break; }
    }
    let mut values = (0..n).map(|i| a[i * n + i]).collect::<Vec<_>>();
    values.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    values
}

fn graph_matrix(graph: &SpectrumGraph) -> Vec<f64> {
    let n = graph.transitions.len();
    let mut a = vec![0.0; n * n];
    for (i, row) in graph.transitions.iter().enumerate() {
        for &j in row { a[i * n + j as usize] += 1.0; }
    }
    a
}

fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() }
fn norm(a: &[f64]) -> f64 { dot(a, a).sqrt() }

fn graph_matvec(graph: &SpectrumGraph, v: &[f64]) -> Vec<f64> {
    graph.transitions.iter().map(|row| row.iter().map(|&j| v[j as usize]).sum()).collect()
}

fn project_constant(v: &mut [f64]) {
    if v.is_empty() { return; }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    for x in v { *x -= mean; }
}

fn lanczos_extremes(graph: &SpectrumGraph, steps: usize, seed: u64) -> (f64, f64) {
    let n = graph.transitions.len();
    if n <= 1 { return (f64::NAN, f64::NAN); }
    let m = steps.min(n - 1).max(1);
    let mut rng = SplitMix64::new(seed);
    let mut q = (0..n).map(|_| (rng.next_u64() as f64 / u64::MAX as f64) - 0.5).collect::<Vec<_>>();
    project_constant(&mut q);
    if norm(&q) < 1e-14 {
        for (i, x) in q.iter_mut().enumerate() { *x = if i % 2 == 0 { 1.0 } else { -1.0 }; }
        project_constant(&mut q);
    }
    let qn = norm(&q);
    for x in &mut q { *x /= qn; }
    let mut previous = vec![0.0; n];
    let mut diagonal = Vec::with_capacity(m);
    let mut off_diagonal = Vec::with_capacity(m.saturating_sub(1));
    let mut beta_previous = 0.0;
    for step in 0..m {
        let mut z = graph_matvec(graph, &q);
        if step > 0 { for (z, p) in z.iter_mut().zip(&previous) { *z -= beta_previous * p; } }
        project_constant(&mut z);
        let alpha = dot(&q, &z);
        for (z, qx) in z.iter_mut().zip(&q) { *z -= alpha * qx; }
        // Three-term Lanczos keeps memory O(n), which is essential for A8
        // components with millions of quotient states. Reprojecting the
        // constant mode prevents the trivial eigenvalue 6 from leaking back.
        project_constant(&mut z);
        let beta = norm(&z);
        diagonal.push(alpha);
        if step + 1 < m {
            off_diagonal.push(beta);
            if beta < 1e-10 { break; }
            previous = q;
            q = z.into_iter().map(|x| x / beta).collect();
            beta_previous = beta;
        }
    }
    let k = diagonal.len();
    let mut tri = vec![0.0; k * k];
    for i in 0..k { tri[i * k + i] = diagonal[i]; }
    for i in 0..k.saturating_sub(1) { tri[i * k + i + 1] = off_diagonal[i]; tri[(i + 1) * k + i] = off_diagonal[i]; }
    let values = jacobi_eigenvalues(tri, k);
    (*values.first().unwrap_or(&f64::NAN), *values.last().unwrap_or(&f64::NAN))
}

fn eigen_for_graph(graph: &SpectrumGraph, full_max: usize, lanczos_steps: usize, seed: u64) -> EigenResult {
    let n = graph.transitions.len();
    if n <= full_max {
        let values = jacobi_eigenvalues(graph_matrix(graph), n);
        let lambda_max = *values.first().unwrap_or(&f64::NAN);
        let lambda2 = *values.get(1).unwrap_or(&f64::NAN);
        let lambda_min = *values.last().unwrap_or(&f64::NAN);
        EigenResult { method: "jacobi-exact".into(), values, lambda_max, lambda2, lambda_min }
    } else {
        let (lambda2, lambda_min) = lanczos_extremes(graph, lanczos_steps, seed);
        EigenResult { method: format!("lanczos-{lanczos_steps}"), values: Vec::new(), lambda_max: 6.0, lambda2, lambda_min }
    }
}

fn random_permutation(n: usize, rng: &mut SplitMix64) -> Vec<usize> {
    let mut p: Vec<usize> = (0..n).collect();
    rng.shuffle(&mut p);
    p
}

fn random_schreier_graph(n: usize, rng: &mut SplitMix64) -> SpectrumGraph {
    assert!(n <= u32::MAX as usize);
    let mut transitions = vec![[0u32; 6]; n];
    for slot in 0..3 {
        let p = random_permutation(n, rng);
        let mut inv = vec![0usize; n];
        for (i, &j) in p.iter().enumerate() { inv[j] = i; }
        for i in 0..n { transitions[i][slot] = p[i] as u32; transitions[i][slot + 3] = inv[i] as u32; }
    }
    SpectrumGraph { states: Vec::new(), transitions }
}

fn cycle_matched_schreier_graph(source: &SpectrumGraph, rng: &mut SplitMix64) -> SpectrumGraph {
    let n = source.transitions.len();
    let generators = positive_generators(source);
    let mut transitions = vec![[0u32; 6]; n];
    for slot in 0..3 {
        let labels = random_permutation(n, rng);
        let mut p = vec![0usize; n];
        for i in 0..n { p[labels[i]] = labels[generators[slot][i]]; }
        let mut inverse = vec![0usize; n];
        for (i, &j) in p.iter().enumerate() { inverse[j] = i; }
        for i in 0..n {
            transitions[i][slot] = p[i] as u32;
            transitions[i][slot + 3] = inverse[i] as u32;
        }
    }
    SpectrumGraph { states: Vec::new(), transitions }
}

// Independent conjugation within each generator's nonfixed support preserves
// its cycle type AND the exact set of loops at every vertex.
fn loop_matched_schreier_graph(source: &SpectrumGraph, rng: &mut SplitMix64) -> SpectrumGraph {
    let n = source.transitions.len();
    let generators = positive_generators(source);
    let mut transitions = vec![[0u32; 6]; n];
    for slot in 0..3 {
        let mut labels: Vec<usize> = (0..n).collect();
        let support: Vec<usize> = (0..n).filter(|&i| generators[slot][i] != i).collect();
        let mut shuffled = support.clone();
        rng.shuffle(&mut shuffled);
        for (from, to) in support.into_iter().zip(shuffled) { labels[from] = to; }
        let mut inverse = vec![0usize; n];
        for i in 0..n {
            let from = labels[i];
            let to = labels[generators[slot][i]];
            transitions[from][slot] = to as u32;
            inverse[to] = from;
        }
        for i in 0..n { transitions[i][slot + 3] = inverse[i] as u32; }
    }
    SpectrumGraph { states: Vec::new(), transitions }
}

fn random_regular_graph(n: usize, degree: usize, switch_attempts_per_edge: usize, rng: &mut SplitMix64) -> Option<SpectrumGraph> {
    if n <= degree || degree != 6 { return None; }
    if n == degree + 1 {
        let transitions = (0..n).map(|i| {
            let mut row = [0u32; 6];
            let mut k = 0;
            for j in 0..n { if i != j { row[k] = j as u32; k += 1; } }
            row
        }).collect();
        return Some(SpectrumGraph { states: Vec::new(), transitions });
    }
    // A configuration-model rejection sampler almost never succeeds for
    // millions of vertices. Start from a randomly relabeled circulant graph
    // and use degree-preserving double-edge switches. This samples the usual
    // simple d-regular switch chain without loops or parallel edges.
    let labels = random_permutation(n, rng);
    let mut transitions = vec![[0u32; 6]; n];
    for i in 0..n {
        for offset in 1..=3 {
            transitions[labels[i]][2 * (offset - 1)] = labels[(i + offset) % n] as u32;
            transitions[labels[i]][2 * (offset - 1) + 1] = labels[(i + n - offset) % n] as u32;
        }
    }
    for _ in 0..(3 * switch_attempts_per_edge * n) {
        let a = rng.below(n); let ai = rng.below(6); let x = transitions[a][ai] as usize;
        let b = rng.below(n); let bi = rng.below(6); let y = transitions[b][bi] as usize;
        if a == b || a == y || b == x || x == y || transitions[a].contains(&(y as u32)) || transitions[b].contains(&(x as u32)) { continue; }
        let Some(xi) = transitions[x].iter().position(|&v| v as usize == a) else { continue; };
        let Some(yi) = transitions[y].iter().position(|&v| v as usize == b) else { continue; };
        transitions[a][ai] = y as u32; transitions[y][yi] = a as u32;
        transitions[b][bi] = x as u32; transitions[x][xi] = b as u32;
    }
    Some(SpectrumGraph { states: Vec::new(), transitions })
}

fn case_seed(case: &SpectrumCase, salt: u64) -> u64 {
    let mut h = 0xcbf29ce484222325u64 ^ salt;
    for byte in case.id.as_bytes() { h = (h ^ *byte as u64).wrapping_mul(0x100000001b3); }
    h
}

fn passport_string(t: &[Perm; 4], d: usize) -> String {
    t.iter().map(|p| cycle_type(*p, d).iter().map(|n| n.to_string()).collect::<Vec<_>>().join(".")).collect::<Vec<_>>().join("/")
}

fn genus_of(t: &[Perm; 4], d: usize) -> i64 {
    let total: i64 = t.iter().map(|p| index(*p, d) as i64).sum();
    (2 - 2 * d as i64 + total) / 2
}

fn defect_of(t: &[Perm; 4], d: usize) -> (i64, i64, i64) {
    // Fixed-size cycle/DSU counters: no heap allocations in a walk's hot loop.
    let product_index = index(t[0].compose(t[1], d), d);
    let cycles = (d - product_index) as i64;
    let left = point_orbit_count_perms(d, t[0], t[1]) as i64;
    let right = point_orbit_count_perms(d, t[2], t[3]) as i64;
    let defect = cycles - left - right + 1;
    let left_genus = (index(t[0], d) + index(t[1], d)
        + product_index) as i64 / 2 - d as i64 + left;
    let right_genus = (index(t[2], d) + index(t[3], d)
        + product_index) as i64 / 2 - d as i64 + right;
    (defect, left_genus, right_genus)
}

fn defect_census(args: &[String]) -> Result<(), String> {
    let mut cases = Vec::new();
    let mut roots = HashMap::new();
    if let Some(raw) = optional_arg(args, "--components") {
        for item in raw.split(',') {
            for case in parse_components_file(Path::new(item))? {
                add_spectrum_case(&mut cases, &mut roots, case);
            }
        }
    }
    if let Some(raw) = optional_arg(args, "--certificates") {
        let mut paths = Vec::new();
        for item in raw.split(',') { collect_certificate_paths(Path::new(item), &mut paths)?; }
        for path in paths {
            verify(&path)?;
            add_spectrum_case(&mut cases, &mut roots, parse_certificate_case(&path)?);
        }
    }
    if has_flag(args, "--exclude-identity") {
        cases.retain(|case| !case.seed.iter().any(|p| p.is_identity(case.d)));
    }
    if cases.is_empty() { return Err("defect-census needs --components and/or --certificates".into()); }
    let out = PathBuf::from(get_arg(args, "--out", None)?);
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let mut rows = BufWriter::new(File::create(out.join("DEFECT_DISTRIBUTIONS.tsv")).map_err(|e| e.to_string())?);
    writeln!(rows, "case_id\tgroup\tdegree\tgenus\tstates\tdefect\tcount\tleft_genus_min\tleft_genus_max\tright_genus_min\tright_genus_max\tsource").map_err(|e| e.to_string())?;
    for case in cases {
        let states = match case.states {
            Some(states) => graph_from_states(case.d, states)?.states,
            None => enumerate_spectrum_graph(case.d, case.seed)?.states,
        };
        if case.input_states.is_some_and(|n| n != states.len()) { return Err(format!("{}: state count mismatch", case.id)); }
        let genus = genus_of(&case.seed, case.d);
        let mut counts: BTreeMap<i64, (usize, i64, i64, i64, i64)> = BTreeMap::new();
        for key in &states {
            let (defect, left, right) = defect_of(&unpack(*key, case.d), case.d);
            if defect < 0 || left < 0 || right < 0 || defect + left + right != genus {
                return Err(format!("{}: defect/genus identity failed", case.id));
            }
            let entry = counts.entry(defect).or_insert((0, left, left, right, right));
            entry.0 += 1;
            entry.1 = entry.1.min(left); entry.2 = entry.2.max(left);
            entry.3 = entry.3.min(right); entry.4 = entry.4.max(right);
        }
        for (defect, (count, lmin, lmax, rmin, rmax)) in counts {
            writeln!(rows, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", case.id, case.ambient.short(), case.d, genus, states.len(), defect, count, lmin, lmax, rmin, rmax, case.source).map_err(|e| e.to_string())?;
        }
        eprintln!("defect census {}: {} states", case.id, states.len());
    }
    Ok(())
}

fn compute_spectrum_case(case: SpectrumCase, full_max: usize, lanczos_steps: usize, options: &SpectrumOptions, summary_only: bool) -> Result<ComputedSpectrum, String> {
    let mut graph = match case.states.clone() {
        Some(states) => graph_from_states(case.d, states)?,
        None => enumerate_spectrum_graph(case.d, case.seed)?,
    };
    let components = graph_components(&graph);
    let (loops, parallel_excess) = graph_shape(&graph);
    if summary_only { graph.states.clear(); graph.states.shrink_to_fit(); }
    let eigen = eigen_for_graph(&graph, full_max, lanczos_steps, case_seed(&case, 0));
    let diagnostics = graph_diagnostics(&graph, &eigen, options.trace_order, options.trace_probes, options.maxcut_restarts, options.canonical_max, case_seed(&case, 0xd1a6));
    let convergence = if graph.transitions.len() > full_max {
        options.lanczos_checks.iter().map(|&steps| {
            let (lambda2, lambda_min) = lanczos_extremes(&graph, steps, case_seed(&case, 0xc011));
            (steps, lambda2, lambda_min)
        }).collect()
    } else { Vec::new() };
    let mut baselines = Vec::new();
    if case.kind == "exact_component" || case.kind == "certificate" {
        for replicate in 0..options.baseline_reps {
            let mut rng = SplitMix64::new(case_seed(&case, replicate as u64 + 1));
            if options.baseline_kinds.contains("random_schreier") {
                let baseline = random_schreier_graph(graph.transitions.len(), &mut rng);
                let components = graph_components(&baseline);
                let eigen = eigen_for_graph(&baseline, full_max, lanczos_steps, case_seed(&case, 0x51 + replicate as u64));
                let diagnostics = graph_diagnostics(&baseline, &eigen, options.trace_order, options.trace_probes, options.maxcut_restarts, options.canonical_max, case_seed(&case, 0x5100 + replicate as u64));
                baselines.push(BaselineResult { kind: "random_schreier".into(), replicate, components, eigen, diagnostics });
            }
            if options.baseline_kinds.contains("cycle_matched_schreier") {
                let baseline = cycle_matched_schreier_graph(&graph, &mut rng);
                let components = graph_components(&baseline);
                let eigen = eigen_for_graph(&baseline, full_max, lanczos_steps, case_seed(&case, 0x71 + replicate as u64));
                let diagnostics = graph_diagnostics(&baseline, &eigen, options.trace_order, options.trace_probes, options.maxcut_restarts, options.canonical_max, case_seed(&case, 0x7100 + replicate as u64));
                baselines.push(BaselineResult { kind: "cycle_matched_schreier".into(), replicate, components, eigen, diagnostics });
            }
            if options.baseline_kinds.contains("loop_matched_schreier") {
                let baseline = loop_matched_schreier_graph(&graph, &mut rng);
                let components = graph_components(&baseline);
                let eigen = eigen_for_graph(&baseline, full_max, lanczos_steps, case_seed(&case, 0x81 + replicate as u64));
                let diagnostics = graph_diagnostics(&baseline, &eigen, options.trace_order, options.trace_probes, options.maxcut_restarts, options.canonical_max, case_seed(&case, 0x8100 + replicate as u64));
                baselines.push(BaselineResult { kind: "loop_matched_schreier".into(), replicate, components, eigen, diagnostics });
            }
            if options.baseline_kinds.contains("random_6_regular") {
                if let Some(baseline) = random_regular_graph(graph.transitions.len(), 6, options.switch_attempts_per_edge, &mut rng) {
                    let components = graph_components(&baseline);
                    let eigen = eigen_for_graph(&baseline, full_max, lanczos_steps, case_seed(&case, 0x61 + replicate as u64));
                    let diagnostics = graph_diagnostics(&baseline, &eigen, options.trace_order, options.trace_probes, options.maxcut_restarts, options.canonical_max, case_seed(&case, 0x6100 + replicate as u64));
                    baselines.push(BaselineResult { kind: "random_6_regular".into(), replicate, components, eigen, diagnostics });
                }
            }
        }
    }
    let vertices = graph.transitions.len();
    Ok(ComputedSpectrum { case, graph, vertices, eigen, components, loops, parallel_excess, diagnostics, convergence, baselines })
}

fn format_float(x: f64) -> String {
    if x.is_finite() { format!("{x:.10}") } else { "nan".into() }
}

fn finite_stats(values: impl IntoIterator<Item = f64>) -> (usize, f64, f64, f64) {
    let (mut n, mut sum, mut min, mut max) = (0, 0.0, f64::INFINITY, f64::NEG_INFINITY);
    for value in values.into_iter().filter(|x| x.is_finite()) {
        n += 1;
        sum += value;
        min = min.min(value);
        max = max.max(value);
    }
    if n == 0 { (0, f64::NAN, f64::NAN, f64::NAN) } else { (n, sum / n as f64, min, max) }
}

fn spectral_wasserstein(actual: &[f64], baseline: &[f64]) -> f64 {
    if actual.len() != baseline.len() || actual.len() <= 1 { return f64::NAN; }
    actual[1..].iter().zip(&baseline[1..]).map(|(a, b)| (a - b).abs()).sum::<f64>() / (actual.len() - 1) as f64
}

fn write_graph_diagnostic(writer: &mut impl Write, case_id: &str, source: &str, replicate: Option<usize>, vertices: usize, diagnostic: &GraphDiagnostics) -> Result<(), String> {
    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:016x}\t{}",
        case_id, source, replicate.map(|x| x.to_string()).unwrap_or_default(), vertices,
        diagnostic.generator_profiles[0], diagnostic.generator_profiles[1], diagnostic.generator_profiles[2],
        diagnostic.generator_fixed[0], diagnostic.generator_fixed[1], diagnostic.generator_fixed[2],
        diagnostic.generator_cycles[0], diagnostic.generator_cycles[1], diagnostic.generator_cycles[2],
        diagnostic.braid_12_defects, diagnostic.braid_23_defects, diagnostic.commute_13_defects,
        diagnostic.bipartite, format_float(diagnostic.maxcut_conflict_fraction), diagnostic.wl_hash, diagnostic.canonical_labeled_code).map_err(|e| e.to_string())
}

fn write_spectrum_outputs(out: &Path, computed: &[ComputedSpectrum], input_rows: usize, options: &SpectrumOptions, full_max: usize, degree5_full_max: usize, lanczos_steps: usize, summary_only: bool) -> Result<(), String> {
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let mut cases = BufWriter::new(File::create(out.join("SPECTRUM.tsv")).map_err(|e| e.to_string())?);
    writeln!(cases, "case_id\tkind\tgroup\tdegree\tsource_rows\tinput_states\tquotient_states\tpassport\tindices\tgenus\tconnected_components\tloop_moves\tparallel_excess\tlambda_max\tlambda2\tlambda_min\tabs_nontrivial_radius\tspectral_gap\tnormalized_gap\tmethod\tsource\taliases").map_err(|e| e.to_string())?;
    let mut values = BufWriter::new(File::create(out.join("EIGENVALUES.tsv")).map_err(|e| e.to_string())?);
    writeln!(values, "case_id\tsource\teigen_index\teigenvalue").map_err(|e| e.to_string())?;
    let mut states = BufWriter::new(File::create(out.join("STATES.tsv")).map_err(|e| e.to_string())?);
    writeln!(states, "case_id\tstate_index\tkey0\tkey1\tkey2").map_err(|e| e.to_string())?;
    let mut edges = BufWriter::new(File::create(out.join("EDGES.tsv")).map_err(|e| e.to_string())?);
    writeln!(edges, "case_id\tfrom\tmove\tto").map_err(|e| e.to_string())?;
    let mut baselines = BufWriter::new(File::create(out.join("BASELINES.tsv")).map_err(|e| e.to_string())?);
    writeln!(baselines, "case_id\tbaseline\treplicate\tvertices\tconnected_components\tlambda_max\tlambda2\tlambda_min\tabs_nontrivial_radius\tspectral_gap\tnormalized_gap\tmethod").map_err(|e| e.to_string())?;
    let mut matched = BufWriter::new(File::create(out.join("MATCHED_BASELINES.tsv")).map_err(|e| e.to_string())?);
    writeln!(matched, "case_id\tgroup\tdegree\tgenus\tvertices\tbaseline\treplicates\tdisconnected_replicates\tactual_lambda2\tmean_baseline_lambda2\tactual_minus_baseline_lambda2\tactual_normalized_gap\tmean_baseline_normalized_gap\tactual_minus_baseline_normalized_gap\tmean_actual_to_baseline_wasserstein\tmean_baseline_pair_wasserstein\twasserstein_ratio").map_err(|e| e.to_string())?;
    let mut diagnostics = BufWriter::new(File::create(out.join("GRAPH_DIAGNOSTICS.tsv")).map_err(|e| e.to_string())?);
    writeln!(diagnostics, "case_id\tsource\treplicate\tvertices\tgenerator_1_cycles\tgenerator_2_cycles\tgenerator_3_cycles\tgenerator_1_fixed\tgenerator_2_fixed\tgenerator_3_fixed\tgenerator_1_cycle_count\tgenerator_2_cycle_count\tgenerator_3_cycle_count\tbraid_12_defects\tbraid_23_defects\tcommute_13_defects\texact_bipartite\tmaxcut_conflict_fraction\twl_hash\tcanonical_labeled_code").map_err(|e| e.to_string())?;
    let mut moments = BufWriter::new(File::create(out.join("TRACE_MOMENTS.tsv")).map_err(|e| e.to_string())?);
    writeln!(moments, "case_id\tsource\treplicate\torder\tnormalized_trace\tstandard_error\tmethod").map_err(|e| e.to_string())?;
    let mut convergence = BufWriter::new(File::create(out.join("LANCZOS_CONVERGENCE.tsv")).map_err(|e| e.to_string())?);
    writeln!(convergence, "case_id\tvertices\tsteps\tlambda2\tlambda_min\tnormalized_gap\tabsolute_gap").map_err(|e| e.to_string())?;
    let mut full_cases = 0;
    let mut iterative_cases = 0;
    let mut total_states = 0usize;
    for item in computed {
        let t = item.case.seed;
        let indices = t.iter().map(|p| index(*p, item.case.d).to_string()).collect::<Vec<_>>().join(";");
        let lambda2 = item.eigen.lambda2;
        let radius = lambda2.abs().max(item.eigen.lambda_min.abs());
        let gap = 6.0 - lambda2;
        let normalized = 1.0 - lambda2 / 6.0;
        let alias_text = item.case.aliases.join("|");
        writeln!(cases, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            item.case.id, item.case.kind, item.case.ambient.short(), item.case.d, item.case.aliases.len().max(1), item.case.input_states.map(|x| x.to_string()).unwrap_or_default(), item.vertices, passport_string(&t, item.case.d), indices, genus_of(&t, item.case.d), item.components, item.loops, item.parallel_excess, format_float(item.eigen.lambda_max), format_float(lambda2), format_float(item.eigen.lambda_min), format_float(radius), format_float(gap), format_float(normalized), item.eigen.method, item.case.source, alias_text).map_err(|e| e.to_string())?;
        write_graph_diagnostic(&mut diagnostics, &item.case.id, "case", None, item.vertices, &item.diagnostics)?;
        for &(order, estimate, standard_error, ref method) in &item.diagnostics.moments {
            writeln!(moments, "{}\tcase\t\t{}\t{}\t{}\t{}", item.case.id, order, format_float(estimate), format_float(standard_error), method).map_err(|e| e.to_string())?;
        }
        for &(steps, check_lambda2, check_lambda_min) in &item.convergence {
            let check_radius = check_lambda2.abs().max(check_lambda_min.abs());
            writeln!(convergence, "{}\t{}\t{}\t{}\t{}\t{}\t{}", item.case.id, item.vertices, steps, format_float(check_lambda2), format_float(check_lambda_min), format_float(1.0 - check_lambda2 / 6.0), format_float(1.0 - check_radius / 6.0)).map_err(|e| e.to_string())?;
        }
        total_states += item.vertices;
        if item.eigen.values.is_empty() { iterative_cases += 1; } else { full_cases += 1; }
        if !summary_only {
            for (i, &key) in item.graph.states.iter().enumerate() {
                writeln!(states, "{}\t{}\t{:032x}\t{:032x}\t{:032x}", item.case.id, i, key[0], key[1], key[2]).map_err(|e| e.to_string())?;
                for mv in 0..6 { writeln!(edges, "{}\t{}\t{}\t{}", item.case.id, i, mv, item.graph.transitions[i][mv]).map_err(|e| e.to_string())?; }
            }
        }
        for (i, &value) in item.eigen.values.iter().enumerate() { writeln!(values, "{}\tcase\t{}\t{}", item.case.id, i, format_float(value)).map_err(|e| e.to_string())?; }
        for baseline in &item.baselines {
            let b = &baseline.eigen;
            let radius = b.lambda2.abs().max(b.lambda_min.abs());
            writeln!(baselines, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", item.case.id, baseline.kind, baseline.replicate, item.vertices, baseline.components, format_float(b.lambda_max), format_float(b.lambda2), format_float(b.lambda_min), format_float(radius), format_float(6.0 - b.lambda2), format_float(1.0 - b.lambda2 / 6.0), b.method).map_err(|e| e.to_string())?;
            write_graph_diagnostic(&mut diagnostics, &item.case.id, &baseline.kind, Some(baseline.replicate), item.vertices, &baseline.diagnostics)?;
            for &(order, estimate, standard_error, ref method) in &baseline.diagnostics.moments {
                writeln!(moments, "{}\t{}\t{}\t{}\t{}\t{}\t{}", item.case.id, baseline.kind, baseline.replicate, order, format_float(estimate), format_float(standard_error), method).map_err(|e| e.to_string())?;
            }
            for (i, &value) in b.values.iter().enumerate() { writeln!(values, "{}\t{}-{}\t{}\t{}", item.case.id, baseline.kind, baseline.replicate, i, format_float(value)).map_err(|e| e.to_string())?; }
        }
        let mut baseline_groups = BTreeMap::<&str, Vec<&BaselineResult>>::new();
        for baseline in &item.baselines { baseline_groups.entry(&baseline.kind).or_default().push(baseline); }
        for (kind, samples) in baseline_groups {
            let (_, mean_lambda2, _, _) = finite_stats(samples.iter().map(|x| x.eigen.lambda2));
            let (_, mean_wasserstein, _, _) = finite_stats(samples.iter().map(|x| spectral_wasserstein(&item.eigen.values, &x.eigen.values)));
            let mut pair_distances = Vec::new();
            for i in 0..samples.len() { for j in i + 1..samples.len() { pair_distances.push(spectral_wasserstein(&samples[i].eigen.values, &samples[j].eigen.values)); } }
            let (_, mean_pair_wasserstein, _, _) = finite_stats(pair_distances);
            let actual_gap = 1.0 - lambda2 / 6.0;
            let baseline_gap = 1.0 - mean_lambda2 / 6.0;
            let ratio = if mean_pair_wasserstein > 0.0 { mean_wasserstein / mean_pair_wasserstein } else { f64::NAN };
            writeln!(matched, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                item.case.id, item.case.ambient.short(), item.case.d, genus_of(&t, item.case.d), item.vertices, kind, samples.len(), samples.iter().filter(|x| x.components != 1).count(), format_float(lambda2), format_float(mean_lambda2), format_float(lambda2 - mean_lambda2), format_float(actual_gap), format_float(baseline_gap), format_float(actual_gap - baseline_gap), format_float(mean_wasserstein), format_float(mean_pair_wasserstein), format_float(ratio)).map_err(|e| e.to_string())?;
        }
    }
    let mut genus_groups = BTreeMap::<(String, usize, i64, String), Vec<&ComputedSpectrum>>::new();
    for item in computed {
        genus_groups.entry((item.case.ambient.short().into(), item.case.d, genus_of(&item.case.seed, item.case.d), item.case.kind.clone())).or_default().push(item);
    }
    let mut genus = BufWriter::new(File::create(out.join("GENUS_SUMMARY.tsv")).map_err(|e| e.to_string())?);
    writeln!(genus, "group\tdegree\tgenus\tkind\tcases\ttotal_vertices\tfull_spectra\titerative_cases\tmean_lambda2\tmin_lambda2\tmax_lambda2\tmean_normalized_gap\tmin_normalized_gap\tmax_normalized_gap").map_err(|e| e.to_string())?;
    for ((group, degree, genus_value, kind), items) in genus_groups {
        let (_, mean_lambda2, min_lambda2, max_lambda2) = finite_stats(items.iter().map(|x| x.eigen.lambda2));
        let (_, mean_gap, min_gap, max_gap) = finite_stats(items.iter().map(|x| 1.0 - x.eigen.lambda2 / 6.0));
        writeln!(genus, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", group, degree, genus_value, kind, items.len(), items.iter().map(|x| x.vertices).sum::<usize>(), items.iter().filter(|x| !x.eigen.values.is_empty()).count(), items.iter().filter(|x| x.eigen.values.is_empty()).count(), format_float(mean_lambda2), format_float(min_lambda2), format_float(max_lambda2), format_float(mean_gap), format_float(min_gap), format_float(max_gap)).map_err(|e| e.to_string())?;
    }
    let mut summary = BufWriter::new(File::create(out.join("SUMMARY.txt")).map_err(|e| e.to_string())?);
    writeln!(summary, "Hurwitz quotient spectral census").map_err(|e| e.to_string())?;
    writeln!(summary, "input component rows = {input_rows}").map_err(|e| e.to_string())?;
    writeln!(summary, "absolute/component cases = {}", computed.len()).map_err(|e| e.to_string())?;
    writeln!(summary, "vertices across cases = {total_states}").map_err(|e| e.to_string())?;
    writeln!(summary, "exact Jacobi spectra = {full_cases}").map_err(|e| e.to_string())?;
    writeln!(summary, "Lanczos extremal spectra = {iterative_cases}").map_err(|e| e.to_string())?;
    writeln!(summary, "full-spectrum vertex threshold = {full_max}").map_err(|e| e.to_string())?;
    writeln!(summary, "degree-5 full-spectrum vertex threshold = {degree5_full_max}").map_err(|e| e.to_string())?;
    writeln!(summary, "Lanczos steps = {lanczos_steps}").map_err(|e| e.to_string())?;
    writeln!(summary, "random baseline replicates per exact component = {}", options.baseline_reps).map_err(|e| e.to_string())?;
    let mut baseline_kinds = options.baseline_kinds.iter().cloned().collect::<Vec<_>>(); baseline_kinds.sort();
    writeln!(summary, "random baseline kinds = {}", baseline_kinds.join(",")).map_err(|e| e.to_string())?;
    writeln!(summary, "trace moments = orders 1..{} with {} probes for iterative spectra", options.trace_order, options.trace_probes).map_err(|e| e.to_string())?;
    writeln!(summary, "max-cut local-search restarts = {}", options.maxcut_restarts).map_err(|e| e.to_string())?;
    writeln!(summary, "exact labeled canonical-code vertex threshold = {}", options.canonical_max).map_err(|e| e.to_string())?;
    writeln!(summary, "regular-null edge-switch attempts per edge = {}", options.switch_attempts_per_edge).map_err(|e| e.to_string())?;
    writeln!(summary, "Lanczos convergence steps = {}", options.lanczos_checks.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")).map_err(|e| e.to_string())?;
    writeln!(summary, "raw graph tables written = {}", !summary_only).map_err(|e| e.to_string())?;
    writeln!(summary, "GENUS_SUMMARY.tsv blocks genus comparisons by group and degree.").map_err(|e| e.to_string())?;
    writeln!(summary, "MATCHED_BASELINES.tsv compares each component with same-size random graphs; Wasserstein distances omit the trivial first eigenvalue.").map_err(|e| e.to_string())?;
    writeln!(summary, "Simple 6-regular baselines use a randomly relabeled circulant start followed by the configured double-edge-switch chain.").map_err(|e| e.to_string())?;
    writeln!(summary, "Cycle-matched Schreier baselines independently conjugate each positive generator, preserving its full cycle type while intentionally breaking braid relations.").map_err(|e| e.to_string())?;
    writeln!(summary, "All graph edges are the six labeled B4 moves (three generators and inverses).").map_err(|e| e.to_string())?;
    writeln!(summary, "The quotient is simultaneous S_d conjugation. Spectral gaps are componentwise; disconnected baselines are retained and reported with their component count.").map_err(|e| e.to_string())?;
    Ok(())
}

fn spectrum(args: &[String]) -> Result<(), String> {
    let components_arg = optional_arg(args, "--components");
    let certificates_arg = optional_arg(args, "--certificates");
    if components_arg.is_none() && certificates_arg.is_none() {
        return Err("spectrum needs --components FILE[,FILE...] and/or --certificates FILE_OR_DIR".into());
    }
    let out = PathBuf::from(get_arg(args, "--out", Some("runs/spectra/spectrum"))?);
    let full_max: usize = parse_num(get_arg(args, "--full-max", Some("256"))?, "--full-max")?;
    let degree5_full_max: usize = parse_num(get_arg(args, "--degree5-full-max", Some(&full_max.to_string()))?, "--degree5-full-max")?;
    let lanczos_steps: usize = parse_num(get_arg(args, "--lanczos-steps", Some("128"))?, "--lanczos-steps")?;
    let baseline_reps: usize = parse_num(get_arg(args, "--baseline-reps", Some("8"))?, "--baseline-reps")?;
    let baseline_kinds = get_arg(args, "--baseline-kinds", Some("random_6_regular,random_schreier"))?.split(',').filter(|x| !x.is_empty()).map(str::to_string).collect::<HashSet<_>>();
    let allowed_baselines = ["random_6_regular", "random_schreier", "cycle_matched_schreier", "loop_matched_schreier"].into_iter().collect::<HashSet<_>>();
    if let Some(kind) = baseline_kinds.iter().find(|kind| !allowed_baselines.contains(kind.as_str())) { return Err(format!("unknown baseline kind {kind}")); }
    let trace_order: usize = parse_num(get_arg(args, "--trace-order", Some("0"))?, "--trace-order")?;
    let trace_probes: usize = parse_num(get_arg(args, "--trace-probes", Some("16"))?, "--trace-probes")?;
    let maxcut_restarts: usize = parse_num(get_arg(args, "--maxcut-restarts", Some("0"))?, "--maxcut-restarts")?;
    let canonical_max: usize = parse_num(get_arg(args, "--canonical-max", Some("0"))?, "--canonical-max")?;
    let switch_attempts_per_edge: usize = parse_num(get_arg(args, "--switch-attempts-per-edge", Some("20"))?, "--switch-attempts-per-edge")?;
    let lanczos_checks = optional_arg(args, "--lanczos-checks").unwrap_or_default().split(',').filter(|x| !x.is_empty()).map(|x| parse_num(x.to_string(), "--lanczos-checks")).collect::<Result<Vec<usize>, _>>()?;
    let requested_threads: usize = parse_num(get_arg(args, "--threads", Some("1"))?, "--threads")?;
    let summary_only = args.iter().any(|x| x == "--summary-only");
    if full_max == 0 || degree5_full_max == 0 || lanczos_steps == 0 || requested_threads == 0 || trace_order > 16 || (trace_order > 0 && trace_probes == 0) {
        return Err("full-max, degree5-full-max, lanczos-steps, and threads must be positive".into());
    }
    let options = SpectrumOptions { baseline_reps, baseline_kinds, trace_order, trace_probes, maxcut_restarts, canonical_max, switch_attempts_per_edge, lanczos_checks };
    let mut cases = Vec::new();
    let mut roots = HashMap::new();
    let mut input_rows = 0;
    if let Some(paths) = components_arg {
        for raw in paths.split(',') {
            let parsed = parse_components_file(Path::new(raw))?;
            input_rows += parsed.len();
            for case in parsed { add_spectrum_case(&mut cases, &mut roots, case); }
        }
    }
    if let Some(raw) = certificates_arg {
        let mut paths = Vec::new();
        for item in raw.split(',').filter(|item| !item.is_empty()) {
            collect_certificate_paths(Path::new(item), &mut paths)?;
        }
        paths.sort();
        let mut names = HashSet::new();
        for path in paths {
            let name = path.file_name().and_then(|x| x.to_str()).unwrap_or("").to_string();
            if !names.insert(name) { continue; }
            input_rows += 1;
            add_spectrum_case(&mut cases, &mut roots, parse_certificate_case(&path)?);
        }
    }
    if let Some(path) = optional_arg(args, "--select-case-ids") {
        let selected = fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?.lines().map(str::trim).filter(|x| !x.is_empty() && !x.starts_with('#')).map(str::to_string).collect::<HashSet<_>>();
        cases.retain(|case| selected.contains(&case.id));
    }
    if cases.is_empty() { return Err("spectrum inputs contained no cases".into()); }
    let threads = requested_threads.min(cases.len());
    eprintln!("spectrum cases={} input_rows={} threads={} full_max={} degree5_full_max={} lanczos_steps={} baselines={} trace_order={} maxcut_restarts={}", cases.len(), input_rows, threads, full_max, degree5_full_max, lanczos_steps, baseline_reps, trace_order, maxcut_restarts);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let finished = std::sync::atomic::AtomicUsize::new(0);
    let completed_work = AtomicU64::new(0);
    let total_work = cases.iter().map(|case| case.input_states.unwrap_or(1) as u64).sum::<u64>();
    let started = Instant::now();
    let results = Mutex::new(Vec::<ComputedSpectrum>::with_capacity(cases.len()));
    let errors = Mutex::new(Vec::<String>::new());
    let cases_ref = &cases;
    thread::scope(|scope| {
        let options_ref = &options;
        for _ in 0..threads {
            let next = &next;
            let results = &results;
            let errors = &errors;
            let cases = cases_ref;
            let finished = &finished;
            let completed_work = &completed_work;
            let options = options_ref;
            scope.spawn(move || loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= cases.len() { break; }
                let work = cases[i].input_states.unwrap_or(1) as u64;
                let case_full_max = if cases[i].d == 5 { degree5_full_max } else { full_max };
                match compute_spectrum_case(cases[i].clone(), case_full_max, lanczos_steps, options, summary_only) {
                    Ok(mut result) => {
                        if summary_only { result.graph.states.clear(); result.graph.transitions.clear(); }
                        results.lock().unwrap().push(result);
                        let done_work = completed_work.fetch_add(work, std::sync::atomic::Ordering::Relaxed) + work;
                        let done_cases = finished.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                        if done_cases % 50 == 0 || done_cases == cases.len() {
                            let elapsed = started.elapsed().as_secs_f64();
                            let fraction = done_work as f64 / total_work.max(1) as f64;
                            let eta = if fraction > 0.0 { elapsed * (1.0 / fraction - 1.0) } else { f64::INFINITY };
                            eprintln!("spectrum progress cases={done_cases}/{} work={:.1}% elapsed={elapsed:.0}s eta={eta:.0}s", cases.len(), 100.0 * fraction.min(1.0));
                        }
                    },
                    Err(error) => errors.lock().unwrap().push(format!("{}: {error}", cases[i].id)),
                }
            });
        }
    });
    let errors = errors.into_inner().unwrap();
    if !errors.is_empty() { return Err(errors.join("; ")); }
    let mut computed = results.into_inner().unwrap();
    computed.sort_by(|a, b| a.case.id.cmp(&b.case.id));
    write_spectrum_outputs(&out, &computed, input_rows, &options, full_max, degree5_full_max, lanczos_steps, summary_only)?;
    eprintln!("spectrum complete: {}", out.display());
    Ok(())
}

fn low_index_seed(d: usize, ambient: Ambient, weak_only:bool, rng: &mut SplitMix64) -> Option<[Perm;4]> {
    low_index_seed_filtered(d,ambient,weak_only,false,rng)
}

fn low_index_seed_filtered(d: usize, ambient: Ambient, weak_only:bool, zero_only:bool, rng: &mut SplitMix64) -> Option<[Perm;4]> {
    // Research heuristic, not a completeness claim. Known A8/A9 failures have
    // small total permutation index. Mix involutions and disjoint 3-cycles,
    // and retain total index <= 2d+2 (genus at most 2 for connected covers).
    for _ in 0..200 {
        let mut t=[Perm::identity(d);4];
        for g in &mut t[..3] {
            let mut points = std::array::from_fn::<_, MAX_D, _>(|i| i);
            rng.shuffle(&mut points[..d]);
            match rng.below(10) {
                0..=5 => {
                    let mut pairs=1+rng.below(d/2);
                    if ambient==Ambient::Alternating {pairs-=pairs%2;}
                    for j in 0..pairs {g.img.swap(points[2*j],points[2*j+1]);}
                }
                6..=8 => {
                    let cycles=1+rng.below(d/3);
                    for j in 0..cycles {
                        let a=points[3*j]; let b=points[3*j+1]; let c=points[3*j+2];
                        g.img[a]=b as u8; g.img[b]=c as u8; g.img[c]=a as u8;
                    }
                }
                _ => *g=random_perm(d,ambient,rng),
            }
        }
        t[3]=product_one_fourth(d,t[0],t[1],t[2]);
        let total:usize=t.iter().map(|&p|index(p,d)).sum();
        if total<2*d-2 || total>2*d+2 {continue;}
        if weak_only && has_transitive_pair(d,&t) {continue;}
        if zero_only && (total==2*d-2 || t.iter().any(|p|p.is_identity(d)) || defect_of(&t,d).0!=0) {continue;}
        if tuple_generates_dynamic(d,ambient,&t) {return Some(t);}
    }
    None
}

type Passport = [Vec<usize>; 4];

fn cycle_type(p: Perm, d: usize) -> Vec<usize> {
    let mut seen = [false; MAX_D];
    let mut lengths = Vec::new();
    for x in 0..d {
        if seen[x] { continue; }
        let (mut y, mut n) = (x, 0);
        while !seen[y] { seen[y] = true; n += 1; y = p.img[y] as usize; }
        lengths.push(n);
    }
    lengths.sort_unstable();
    lengths
}

fn parse_passport(text: &str, d: usize, ambient: Ambient) -> Result<Passport,String> {
    let parts: Vec<_> = text.split('/').collect();
    if parts.len() != 4 { return Err("passport needs four slash-separated cycle types, e.g. 3.3/3.3/3.3/3.3; fixed points are implicit".into()); }
    let mut out: Passport = std::array::from_fn(|_| Vec::new());
    for (i, part) in parts.iter().enumerate() {
        for field in part.split('.') {
            let n: usize = field.parse().map_err(|_| "invalid cycle length")?;
            if n == 0 || n > d { return Err("cycle lengths must be in 1..d".into()); }
            out[i].push(n);
        }
        let total: usize = out[i].iter().sum();
        if total > d { return Err("cycle type exceeds degree".into()); }
        out[i].extend(std::iter::repeat(1).take(d-total));
        out[i].sort_unstable();
    }
    let indices: Vec<_> = out.iter().map(|x| d-x.len()).collect();
    let total: usize = indices.iter().sum();
    if total % 2 != 0 || total < 2*d-2 { return Err("passport cannot be a transitive product-one tuple (parity/index obstruction)".into()); }
    if ambient == Ambient::Alternating && indices.iter().any(|n| n%2 != 0) {
        return Err("odd conjugacy class in A_d passport".into());
    }
    if ambient == Ambient::Symmetric && indices.iter().all(|n| n%2 == 0) {
        return Err("all-even passport cannot generate S_d".into());
    }
    Ok(out)
}

fn perm_of_type(lengths: &[usize], d: usize, rng: Option<&mut SplitMix64>) -> Perm {
    let mut points = std::array::from_fn::<_, MAX_D, _>(|i| i);
    if let Some(rng) = rng { rng.shuffle(&mut points[..d]); }
    let mut g = Perm::identity(d);
    let mut offset = 0;
    for &n in lengths {
        for j in 0..n { g.img[points[offset+j]] = points[offset+(j+1)%n] as u8; }
        offset += n;
    }
    g
}

fn passport_seed(d: usize, ambient: Ambient, passport: &Passport, weak_only:bool, rng: &mut SplitMix64) -> Option<[Perm;4]> {
    passport_seed_filtered(d,ambient,passport,weak_only,false,rng)
}

fn passport_seed_filtered(d: usize, ambient: Ambient, passport: &Passport, weak_only:bool, zero_only:bool, rng: &mut SplitMix64) -> Option<[Perm;4]> {
    // Fixing the first conjugacy-class representative loses no S_d-quotient
    // components. This is still sampling g2/g3, NOT an exhaustive passport census.
    let a = perm_of_type(&passport[0], d, None);
    for _ in 0..200 {
        let b = perm_of_type(&passport[1], d, Some(rng));
        let c = perm_of_type(&passport[2], d, Some(rng));
        let t = [a,b,c,product_one_fourth(d,a,b,c)];
        if weak_only && has_transitive_pair(d,&t) {continue;}
        if zero_only && (genus_of(&t,d)<=0 || t.iter().any(|p|p.is_identity(d)) || defect_of(&t,d).0!=0) {continue;}
        if cycle_type(t[3],d) == passport[3] && tuple_generates_dynamic(d,ambient,&t) { return Some(t); }
    }
    None
}

// Construct sigma_1-fixed quotient states, rather than hoping to hit them.
// For involutions a,b, find h with (a,b)^h=(aba,a). Then c in C(h)
// gives (a,b,c,(abc)^-1)^h = sigma_1(a,b,c,(abc)^-1).
// Components of the two-coloured (a,b) graph make the conjugator inexpensive.
// Structural obstruction: h centralizes ab, c, and hence d=(abc)^-1.
// Every permutation commuting with h preserves the set of points in h-cycles
// of each length. If h has different cycle lengths, <c,d> is intransitive.
// In A25, h has type 2^5 3^5: this explains the observed 10+15 partition.
// This proves a bad pair, NOT a bad entire braid orbit; certify closure below.
// If sigma_2 is fixed too, its relabeling k sends c->b and fixes a,d,
// whereas h sends b->a and fixes c,d. Hence ab~ac~bc and cd~bd~ad;
// ab=(cd)^-1 makes all six product cycle types equal. This explains the
// fingerprint in the rigid example without assuming it for general failures.
fn pair_braid_conjugator(d: usize, a: Perm, b: Perm, rng: &mut SplitMix64) -> Option<Perm> {
    let source = [a,b];
    let target = [b.conjugate_by(a,d),a];
    let mut h = Perm::identity(d);
    let mut assigned = [false; MAX_D];
    let mut used = [false; MAX_D];
    for root in 0..d {
        if assigned[root] { continue; }
        let mut choices = Vec::new();
        for image_root in 0..d {
            if used[image_root] { continue; }
            let mut images = [u8::MAX; MAX_D];
            let mut inverse = [u8::MAX; MAX_D];
            let mut queue = vec![root];
            images[root] = image_root as u8;
            inverse[image_root] = root as u8;
            let (mut head, mut valid) = (0, true);
            while head < queue.len() && valid {
                let x = queue[head]; head += 1;
                for j in 0..2 {
                    let xx = source[j].img[x] as usize;
                    let yy = target[j].img[images[x] as usize] as usize;
                    if used[yy] || (images[xx] != u8::MAX && images[xx] as usize != yy)
                        || (inverse[yy] != u8::MAX && inverse[yy] as usize != xx) {
                        valid = false; break;
                    }
                    if images[xx] == u8::MAX {
                        images[xx] = yy as u8; inverse[yy] = xx as u8; queue.push(xx);
                    }
                }
            }
            if valid { choices.push((queue,images)); }
        }
        if choices.is_empty() { return None; }
        let (points,images) = &choices[rng.below(choices.len())];
        for &x in points { h.img[x] = images[x]; assigned[x] = true; used[images[x] as usize] = true; }
    }
    debug_assert_eq!(source.map(|g|g.conjugate_by(h,d)),target);
    Some(h)
}

// An involution commuting with h fixes an h-cycle (identity or a half-turn),
// or exchanges two equal-length cycles with one of length many phases.
// Count completions before sampling, so impossible transposition counts are
// rejected immediately. Counts fit u128 through degree 34. The overall
// construction is deliberately NOT a uniform tuple draw.
fn commuting_counts(lengths: &[usize], m: usize, memo: &mut HashMap<(Vec<usize>,usize),u128>) -> u128 {
    if lengths.is_empty() { return u128::from(m==0); }
    if 2*m > lengths.iter().sum::<usize>() { return 0; }
    let key = (lengths.to_vec(),m);
    if let Some(&n) = memo.get(&key) { return n; }
    let l = lengths[0];
    let mut n = commuting_counts(&lengths[1..],m,memo);
    if l%2==0 && m>=l/2 { n += commuting_counts(&lengths[1..],m-l/2,memo); }
    let partners = lengths[1..].iter().take_while(|&&v|v==l).count();
    if partners>0 && m>=l { n += (partners*l) as u128*commuting_counts(&lengths[2..],m-l,memo); }
    memo.insert(key,n);
    n
}

fn below_u128(rng: &mut SplitMix64, n: u128) -> u128 {
    assert!(n>0);
    if n<=usize::MAX as u128 {return rng.below(n as usize) as u128;}
    let zone=u128::MAX-u128::MAX%n;
    loop {
        let value=((rng.next_u64() as u128)<<64) | rng.next_u64() as u128;
        if value<zone {return value%n;}
    }
}

fn permutation_cycles(d: usize, h: Perm) -> Vec<Vec<usize>> {
    let mut seen = [false; MAX_D];
    let mut cycles = Vec::new();
    for root in 0..d {
        if seen[root] { continue; }
        let mut cycle = Vec::new();
        let mut x = root;
        while !seen[x] { seen[x] = true; cycle.push(x); x = h.img[x] as usize; }
        cycles.push(cycle);
    }
    cycles.sort_by_key(Vec::len);
    cycles
}

fn commuting_involution(d: usize, h: Perm, mut m: usize, rng: &mut SplitMix64) -> Option<Perm> {
    let mut cycles = permutation_cycles(d,h);
    let mut memo = HashMap::new();
    let mut c = Perm::identity(d);
    while !cycles.is_empty() {
        let lengths: Vec<_> = cycles.iter().map(Vec::len).collect();
        let total = commuting_counts(&lengths,m,&mut memo);
        if total==0 { return None; }
        let l = lengths[0];
        let fixed = commuting_counts(&lengths[1..],m,&mut memo);
        let half = if l%2==0 && m>=l/2 {commuting_counts(&lengths[1..],m-l/2,&mut memo)} else {0};
        let mut choice = below_u128(rng,total);
        let first = cycles.remove(0);
        if choice < fixed { continue; }
        choice -= fixed;
        if choice < half {
            for i in 0..l { c.img[first[i]] = first[(i+l/2)%l] as u8; }
            m -= l/2;
        } else {
            choice -= half;
            let completions = commuting_counts(&lengths[2..],m-l,&mut memo);
            let option = (choice/completions) as usize;
            let second = cycles.remove(option/l);
            let phase = option%l;
            for i in 0..l {
                c.img[first[i]] = second[(i+phase)%l] as u8;
                c.img[second[(i+phase)%l]] = first[i] as u8;
            }
            m -= l;
        }
    }
    debug_assert_eq!(m,0);
    debug_assert_eq!(c.compose(h,d),h.compose(c,d));
    Some(c)
}

fn enumerate_commuting_involutions(cycles: &[Vec<usize>], m: usize, c: Perm,
    visit: &mut impl FnMut(Perm)->Result<bool,String>) -> Result<bool,String> {
    if cycles.is_empty() { return if m==0 {visit(c)} else {Ok(true)}; }
    if 2*m>cycles.iter().map(Vec::len).sum::<usize>() { return Ok(true); }
    let first=&cycles[0]; let l=first.len();
    if !enumerate_commuting_involutions(&cycles[1..],m,c,visit)? { return Ok(false); }
    if l%2==0 && m>=l/2 {
        let mut half=c;
        for i in 0..l {half.img[first[i]]=first[(i+l/2)%l] as u8;}
        if !enumerate_commuting_involutions(&cycles[1..],m-l/2,half,visit)? {return Ok(false);}
    }
    if m>=l {
        for j in 1..cycles.len() {
            if cycles[j].len()!=l {break;}
            let remaining:Vec<_>=cycles[1..].iter().enumerate().filter(|(k,_)|k+1!=j).map(|(_,v)|v.clone()).collect();
            for phase in 0..l {
                let mut paired=c;
                for i in 0..l {
                    paired.img[first[i]]=cycles[j][(i+phase)%l] as u8;
                    paired.img[cycles[j][(i+phase)%l]]=first[i] as u8;
                }
                if !enumerate_commuting_involutions(&remaining,m-l,paired,visit)? {return Ok(false);}
            }
        }
    }
    Ok(true)
}

// C(h) is a product of wreath products C_l wr S_r. A t-cycle of h-cycles
// has order t*l/gcd(l, sum(phases)). Anchor the least unused cycle: every
// commuting permutation of order dividing q is emitted exactly once.
fn gcd_usize(mut a:usize,mut b:usize)->usize {while b!=0 {(a,b)=(b,a%b);} a}

fn commuting_power_count(lengths:&[usize],q:usize,memo:&mut HashMap<Vec<usize>,u128>)->u128 {
    if lengths.is_empty() {return 1;}
    if let Some(&n)=memo.get(lengths) {return n;}
    let l=lengths[0];let r=lengths.iter().take_while(|&&n|n==l).count();
    let mut count=0;let mut arrangements=1u128;
    for t in 1..=q.min(r) {
        if t>1 {arrangements*=((r-t+1)*l) as u128;}
        if q%t!=0 {continue;}
        let phases=(0..l).filter(|&s|q%(t*l/gcd_usize(l,s))==0).count() as u128;
        count+=arrangements*phases*commuting_power_count(&lengths[t..],q,memo);
    }
    memo.insert(lengths.to_vec(),count);count
}

fn enumerate_commuting_power(cycles: &[Vec<usize>], used: u64, q: usize, c: Perm,
    deadline: Instant, skip:&mut u128,memo:&mut HashMap<Vec<usize>,u128>,
    visit: &mut impl FnMut(Perm)->Result<bool,String>) -> Result<bool,String> {
    if Instant::now()>=deadline {return Ok(false);}
    if *skip>0 {
        let lengths:Vec<_>=cycles.iter().enumerate().filter(|(i,_)|used&(1<<i)==0).map(|(_,v)|v.len()).collect();
        let count=commuting_power_count(&lengths,q,memo);
        if *skip>=count {*skip-=count;return Ok(true);}
    }
    let Some(first)=(0..cycles.len()).find(|&i|used&(1<<i)==0) else {return visit(c);};
    let l=cycles[first].len();
    fn phases(cycles:&[Vec<usize>], selected:&[usize], pos:usize, sum:usize, q:usize,
        used:u64, c:Perm, deadline:Instant,skip:&mut u128,memo:&mut HashMap<Vec<usize>,u128>,
        visit:&mut impl FnMut(Perm)->Result<bool,String>) -> Result<bool,String> {
        let l=cycles[selected[0]].len(); let t=selected.len();
        for phase in 0..l {
            if pos+1==t && q%(t*l/gcd_usize(l,(sum+phase)%l))!=0 {continue;}
            let mut next=c;
            for i in 0..l {next.img[cycles[selected[pos]][i]]=cycles[selected[(pos+1)%t]][(i+phase)%l] as u8;}
            let complete=if pos+1==t {enumerate_commuting_power(cycles,used,q,next,deadline,skip,memo,visit)?}
                else {phases(cycles,selected,pos+1,sum+phase,q,used,next,deadline,skip,memo,visit)?};
            if !complete {return Ok(false);}
        }
        Ok(true)
    }
    fn select(cycles:&[Vec<usize>], selected:&mut Vec<usize>, t:usize, used:u64,
        q:usize,c:Perm,deadline:Instant,skip:&mut u128,memo:&mut HashMap<Vec<usize>,u128>,
        visit:&mut impl FnMut(Perm)->Result<bool,String>) -> Result<bool,String> {
        if Instant::now()>=deadline {return Ok(false);}
        if selected.len()==t {return phases(cycles,selected,0,0,q,used,c,deadline,skip,memo,visit);}
        let l=cycles[selected[0]].len();
        for i in 0..cycles.len() {
            if used&(1<<i)!=0 || cycles[i].len()!=l {continue;}
            selected.push(i);
            let complete=select(cycles,selected,t,used|(1<<i),q,c,deadline,skip,memo,visit)?;
            selected.pop();
            if !complete {return Ok(false);}
        }
        Ok(true)
    }
    let mut selected=vec![first];
    for t in 1..=q.min(cycles.iter().filter(|v|v.len()==l).count()) {
        if q%t==0 && !select(cycles,&mut selected,t,used|(1<<first),q,c,deadline,skip,memo,visit)? {return Ok(false);}
    }
    Ok(true)
}

fn power_identity(d:usize,p:Perm,q:usize)->bool {
    (0..q).fold(Perm::identity(d),|a,_|a.compose(p,d)).is_identity(d)
}

fn preserves_partition(d:usize,p:Perm,parts:&[Vec<usize>])->bool {
    let mut labels=[0;MAX_D];
    for (i,part) in parts.iter().enumerate() {for &x in part {labels[x]=i;}}
    parts.iter().all(|part|part.iter().all(|&x|labels[p.img[x] as usize]==labels[p.img[part[0]] as usize]))
        && d==parts.iter().map(Vec::len).sum::<usize>()
}

fn known_weak_states(args:&[String],d:usize)->Result<HashSet<Key>,String> {
    let mut states=HashSet::new();
    if let Ok(dir)=get_arg(args,"--known-certificates",None) {
        for entry in fs::read_dir(dir).map_err(|e|e.to_string())? {
            let path=entry.map_err(|e|e.to_string())?.path();
            if !path.file_name().unwrap().to_string_lossy().starts_with("CERTIFIED_") {continue;}
            let text=fs::read_to_string(&path).map_err(|e|e.to_string())?;
            if !text.lines().any(|s|s==format!("d={d}")) {continue;}
            verify(&path)?;
            let case=parse_certificate_case(&path)?;
            if case.d==d {states.extend(case.states.unwrap_or_default());}
        }
    }
    Ok(states)
}

fn braid_restrict(args: &[String]) -> Result<(),String> {
    let input=PathBuf::from(get_arg(args,"--certificate",None)?);
    verify(&input)?;
    let text=fs::read_to_string(&input).map_err(|e|e.to_string())?;
    let value=|prefix|text.lines().find_map(|s|s.strip_prefix(prefix)).ok_or("missing certificate field");
    let d:usize=value("d=")?.parse().map_err(|_|"invalid degree")?;
    let t=parse_tuple(value("seed=")?,d)?;
    if !(0..3).all(|i|t[i].compose(t[i],d).is_identity(d) && cycle_type(t[i],d)==cycle_type(t[0],d)) {
        return Err("requires first three entries to be equal-type involutions".into());
    }
    let seconds:u64=parse_num(get_arg(args,"--seconds",Some("300"))?,"--seconds")?;
    let cap:usize=parse_num(get_arg(args,"--cap",Some("200000"))?,"--cap")?;
    let limit:u64=parse_num(get_arg(args,"--max-involutions",Some("1000000"))?,"--max-involutions")?;
    if seconds==0 || cap==0 || limit==0 {return Err("positive limits required".into());}
    let h=tuple_conjugator(d,t,hurwitz_move_dynamic(d,t,0))?;
    let blocks=point_orbits(d,&[t[0],t[1],h]);
    if blocks.len()>16 {return Err("too many invariant subsets for this finite experiment".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let create=|name|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|format!("fresh --out required: {e}"));
    let mut rows=create("RESTRICTIONS.tsv")?;
    writeln!(rows,"mask\tdegree\tgroup\ttranspositions\tcompatible_thirds\tvisited\tgenerating\talso_sigma2_fixed\tweak_success\tweak_failure_draws\tunresolved\tcomplete").map_err(|e|e.to_string())?;
    let mut results=create("RESULTS.tsv")?;
    writeln!(results,"status\tquotient_states\tweak\tseed").map_err(|e|e.to_string())?;
    let mut rigid=create("RIGID.tsv")?;
    writeln!(rigid,"degree\tgroup\tgenus\tstatus\tstates\tseed").map_err(|e|e.to_string())?;
    let mut summary=create("SUMMARY.txt")?;
    writeln!(summary,"experiment = braid symmetry restriction\nsource = {}\nconjugator = {}\nconjugator cycle type = {:?}\ninvariant blocks (one-based) = {:?}\nseconds = {seconds}\ncap = {cap}\nmax involutions per subset = {limit}",input.display(),h.cycle_string(d),cycle_type(h,d),blocks.iter().map(|b|b.iter().map(|x|x+1).collect::<Vec<_>>()).collect::<Vec<_>>()).map_err(|e|e.to_string())?;
    let deadline=Instant::now()+Duration::from_secs(seconds);
    let mut all_complete=true;
    let mut certificates=0;
    let mut closed=HashSet::new();
    for mask in 1..(1usize<<blocks.len()) {
        let mut points:Vec<_>=blocks.iter().enumerate().filter(|(i,_)|mask&(1<<i)!=0).flat_map(|(_,v)|v.iter().copied()).collect();
        points.sort_unstable();
        let n=points.len();
        if n<4 {continue;}
        let mut labels=[u8::MAX;MAX_D];
        for (i,&x) in points.iter().enumerate() {labels[x]=i as u8;}
        let restrict=|p:Perm| {
            let mut q=Perm::identity(n);
            for (i,&x) in points.iter().enumerate() {q.img[i]=labels[p.img[x] as usize];}
            q
        };
        let (a,b,hh)=(restrict(t[0]),restrict(t[1]),restrict(h));
        let m=index(a,n);
        let group=if m%2==0 {Ambient::Alternating} else {Ambient::Symmetric};
        let expected=commuting_counts(&cycle_type(hh,n),m,&mut HashMap::new());
        let (mut visited,mut generating,mut fixed,mut success,mut failure,mut unresolved)=(0u64,0u64,0u64,0u64,0u64,0u64);
        let complete=if expected>limit as u128 || Instant::now()>=deadline {false} else {
            enumerate_commuting_involutions(&permutation_cycles(n,hh),m,Perm::identity(n),&mut |c| {
                if Instant::now()>=deadline {return Ok(false);}
                visited+=1;
                let seed=[a,b,c,product_one_fourth(n,a,b,c)];
                if !tuple_generates_dynamic(n,group,&seed) {return Ok(true);}
                generating+=1;
                let key=canonical(&seed,n);
                let second_fixed=key==canonical(&hurwitz_move_dynamic(n,seed,1),n);
                fixed+=u64::from(second_fixed);
                let result=explore_goal(n,seed,cap,deadline.min(Instant::now()+Duration::from_secs(5)),true);
                if result.weak_witness.is_some() {success+=1;}
                else if result.status.starts_with("certified_") {
                    failure+=1;
                    let representative=*result.states.iter().min().unwrap();
                    if closed.insert((n,representative)) {
                        let path=out.join(format!("CERTIFIED_{}{}_{}.txt",group.short(),n,certificates));
                        write_certificate(&path,n,group,&seed,&result).map_err(|e|e.to_string())?;
                        verify(&path)?;
                        certificates+=1;
                    }
                } else {unresolved+=1;}
                writeln!(results,"{}\t{}\t{}\t{}",result.status,result.states.len(),result.weak_witness.is_some(),tuple_string_dynamic(n,&seed)).map_err(|e|e.to_string())?;
                if second_fixed {
                    writeln!(rigid,"{n}\t{}\t{}\t{}\t{}\t{}",group.short(),genus_of(&seed,n),result.status,result.states.len(),tuple_string_dynamic(n,&seed)).map_err(|e|e.to_string())?;
                }
                Ok(true)
            })?
        };
        all_complete &= complete && visited as u128==expected && unresolved==0;
        writeln!(rows,"{mask}\t{n}\t{}\t{m}\t{expected}\t{visited}\t{generating}\t{fixed}\t{success}\t{failure}\t{unresolved}\t{}",group.short(),complete && visited as u128==expected && unresolved==0).map_err(|e|e.to_string())?;
        eprintln!("restriction {}{n}: visited={visited}/{expected} generating={generating} sigma2_fixed={fixed} weak_failures={failure} unresolved={unresolved}",group.short());
    }
    writeln!(summary,"all subsets resolved = {all_complete}\ncertified weak failure components = {certificates}\nOnly invariant restrictions of this seed and this conjugator are covered; not a degree/passport census.").map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())?;
    Ok(())
}

// Every pair of fixed-point-free involutions with these two-colour even-cycle
// sizes is conjugate to this pair. Its centralizer gives all sigma_1 relabelings.
fn matching_pair(sizes:&[usize]) -> (Perm,Perm) {
    let d=sizes.iter().sum();
    let (mut a,mut b)=(Perm::identity(d),Perm::identity(d));
    let mut o=0;
    for &size in sizes {
        assert!(size>=2 && size%2==0);
        for x in (0..size).step_by(2) {a.img[o+x]=(o+x+1) as u8;a.img[o+x+1]=(o+x) as u8;}
        for x in (1..size).step_by(2) {let y=(x+1)%size;b.img[o+x]=(o+y) as u8;b.img[o+y]=(o+x) as u8;}
        o+=size;
    }
    (a,b)
}

fn matching_pair_with_paths(sizes: &[usize], paths: &[usize]) -> (Perm,Perm) {
    let (mut a,mut b)=matching_pair(sizes);
    let mut offset:usize=sizes.iter().sum();
    for &size in paths {
        assert!(size>0 && size%2==1);
        for x in offset..offset+size {a.img[x]=x as u8;b.img[x]=x as u8;}
        for x in 0..size-1 {
            let p=if x%2==0 {&mut a} else {&mut b};
            p.img[offset+x]=(offset+x+1) as u8;
            p.img[offset+x+1]=(offset+x) as u8;
        }
        offset+=size;
    }
    (a,b)
}

fn matching_pair_with_even_paths(sizes:&[usize], paths:&[usize], even:&[usize]) -> (Perm,Perm) {
    let (mut a,mut b)=matching_pair_with_paths(sizes,paths);
    let mut offset:usize=sizes.iter().chain(paths).sum();
    // Each label adds TWO paths of the same even size, with opposite colors.
    // A Hurwitz move exchanges their colored isomorphism types. Matching the
    // multiplicities at EACH size is necessary for a pair-fixed relabeling.
    for &size in even {
        assert!(size>=2 && size%2==0);
        for orientation in 0..2 {
            for x in offset..offset+size {a.img[x]=x as u8;b.img[x]=x as u8;}
            for x in 0..size-1 {
                let p=if (x+orientation)%2==0 {&mut a} else {&mut b};
                p.img[offset+x]=(offset+x+1) as u8;p.img[offset+x+1]=(offset+x) as u8;
            }
            offset+=size;
        }
    }
    (a,b)
}

fn enumerate_pair_centralizer(blocks:&[Vec<usize>],a:Perm,b:Perm,block:usize,used:u64,g:Perm,
    visit:&mut impl FnMut(Perm)->Result<bool,String>) -> Result<bool,String> {
    if block==blocks.len() {return visit(g);}
    let root=blocks[block][0];let size=blocks[block].len();
    for target in 0..blocks.len() {
        if used&(1<<target)!=0 || blocks[target].len()!=size {continue;}
        for &y in &blocks[target] {
            let mut map=[u8::MAX;MAX_D];let mut queue=vec![root];map[root]=y as u8;
            let mut head=0;let mut valid=true;
            while head<queue.len() && valid {
                let x=queue[head];head+=1;
                for p in [a,b] {
                    let xx=p.img[x] as usize;let yy=p.img[map[x] as usize];
                    if map[xx]==u8::MAX {map[xx]=yy;queue.push(xx);}
                    else if map[xx]!=yy {valid=false;break;}
                }
            }
            if !valid || queue.len()!=size {continue;}
            let mut next=g;
            for &x in &blocks[block] {next.img[x]=map[x];}
            if !enumerate_pair_centralizer(blocks,a,b,block+1,used|(1<<target),next,visit)? {return Ok(false);}
        }
    }
    Ok(true)
}

fn hexagon_search(args:&[String])->Result<(),String> {
    let profile=has_flag(args,"--pair-blocks") || has_flag(args,"--pair-paths") || has_flag(args,"--even-paths");
    let k=if profile {0} else {parse_num(get_arg(args,"--k",None)?,"--k")?};
    if !profile && !(2..=5).contains(&k) {return Err("hexagon-search requires k=2..5 (degree <=30)".into());}
    let sizes:Vec<usize>=if profile {
        get_arg(args,"--pair-blocks",Some(""))?.split(',').filter(|s|!s.is_empty()).map(|s|s.parse::<usize>().map_err(|_|"invalid pair block size")).collect::<Result<_,_>>()?
    } else {vec![6;k]};
    let paths:Vec<usize>=get_arg(args,"--pair-paths",Some(""))?.split(',').filter(|s|!s.is_empty())
        .map(|s|s.parse().map_err(|_|"invalid pair path size")).collect::<Result<_,_>>()?;
    let even:Vec<usize>=get_arg(args,"--even-paths",Some(""))?.split(',').filter(|s|!s.is_empty())
        .map(|s|s.parse().map_err(|_|"invalid even path size")).collect::<Result<_,_>>()?;
    let d=sizes.iter().chain(&paths).chain(&even).chain(&even)
        .try_fold(0usize,|sum,&n|sum.checked_add(n)).ok_or("pair degree overflow")?;
    if !(4..=MAX_D).contains(&d) || sizes.iter().any(|&s|s<2 || s%2!=0)
        || paths.iter().any(|&s|s==0 || s%2==0) || even.iter().any(|&s|s<2 || s%2!=0) || (profile && has_flag(args,"--k")) {
        return Err("pair blocks must be positive even cycles, pair paths positive odd paths, total degree 4..34; do not combine with --k".into());
    }
    let m=(d-paths.len()-2*even.len())/2;
    let third_order:usize=parse_num(get_arg(args,"--third-order",Some("2"))?,"--third-order")?;
    let target_group=get_arg(args,"--group",Some("auto"))?;
    if !(2..=6).contains(&third_order) || !["A","S","auto"].contains(&target_group.as_str())
        || (third_order!=2 && (has_flag(args,"--third-pairs") || !profile)) {
        return Err("third-order must be 2..6; higher orders require a pair profile, no third-pairs; group is A/S/auto".into());
    }
    let third_pairs:usize=parse_num(get_arg(args,"--third-pairs",Some(&m.to_string()))?,"--third-pairs")?;
    if third_pairs>d/2 || (has_flag(args,"--third-pairs") && !profile)
        || (target_group=="A" && (m%2!=0 || (third_order==2 && third_pairs%2!=0))) {
        return Err("third-pairs requires a pair profile, 0..d/2; an alternating tuple requires even third-pairs".into());
    }
    let all_h=profile || has_flag(args,"--all-h");
    let quotient_h=has_flag(args,"--quotient-h");
    if quotient_h && !profile {return Err("quotient-h requires a fixed pair profile".into());}
    let all_types=profile || has_flag(args,"--all-types");
    let seconds:u64=parse_num(get_arg(args,"--seconds",Some("600"))?,"--seconds")?;
    let cap:usize=parse_num(get_arg(args,"--cap",Some("200000"))?,"--cap")?;
    if seconds==0 || cap==0 {return Err("positive limits required".into());}
    let start_h:u64=parse_num(get_arg(args,"--start-h",Some("0"))?,"--start-h")?;
    let start_third:u64=parse_num(get_arg(args,"--start-third",Some("0"))?,"--start-third")?;
    let max_thirds:u64=parse_num(get_arg(args,"--max-thirds",Some("18446744073709551615"))?,"--max-thirds")?;
    if max_thirds==0 || ((start_h!=0 || start_third!=0) && !profile) {return Err("positive max-thirds and explicit pair profile required for resumption".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let create=|name|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|format!("fresh --out required: {e}"));
    let mut summary=create("SUMMARY.txt")?;
    let mut results=create("RESULTS.tsv")?;
    writeln!(results,"status\tstates\tseed").map_err(|e|e.to_string())?;
    let mut survivors=create("survivors.txt")?;
    let mut hard=create("hard_candidates.txt")?;
    let (a,b)=matching_pair_with_even_paths(&sizes,&paths,&even);
    let blocks=point_orbits(d,&[a,b]);
    if blocks.len()>16 {return Err("too many pair components".into());}
    let mut rng=SplitMix64::new(20260926);
    let h0=pair_braid_conjugator(d,a,b,&mut rng).ok_or("hexagon pair has no braid relabeling")?;
    let target=[b.conjugate_by(a,d),a];
    let ambient=if m%2==0 && third_pairs%2==0 {Ambient::Alternating} else {Ambient::Symmetric};
    let mut exceptional=vec![1;2*k];exceptional.extend(vec![4;k]);exceptional.sort();
    if !all_h && !all_types && k%3!=0 {
        writeln!(summary,"construction = k identical hexagon components of <a,b>\ngroup = {}_{}\npassport = (2^{m})^3 (4^{k} 1^{})\nsearch complete = true\nstructural exclusion = sigma1 relabeling has only 3-cycles, so the commuting fourth entry's 2k fixed points form a union of 3-cycles; 3 must divide k\nrelabelings visited = 0\ncompatible third involutions visited = 0\nfull-group seeds = 0\ncertified weak-failure components = 0\nThis excludes only the order-three relabeling slice.",ambient.short(),d,2*k).map_err(|e|e.to_string())?;
        summary.sync_all().map_err(|e|e.to_string())?;
        return Ok(());
    }
    let start=Instant::now();let deadline=start+Duration::from_secs(seconds);
    let (mut relabelings,mut eligible_h,mut thirds,mut matching,mut generating,mut good,mut bad,mut unresolved)=(0u64,0u64,0u64,0u64,0u64,0u64,0u64,0u64);
    let mut seen=HashSet::new();let mut closed=known_weak_states(args,d)?;
    let mut direct_good=0u64;
    let (mut index_pruned,mut block_pruned,mut known_hits)=(0u64,0u64,0u64);
    let mut completed_h=HashSet::new();let mut equivalent_h=0u64;
    let mut genera=BTreeMap::<i64,u64>::new();
    let mut passports=BTreeMap::<String,u64>::new();
    let mut group_cache=HashMap::new();
    let (mut ordinal,mut next_h,mut next_third)=(0u64,start_h,start_third);
    let complete=enumerate_pair_centralizer(&blocks,a,b,0,0,Perm::identity(d),&mut |g| {
        let current=ordinal;ordinal+=1;
        if current<start_h {return Ok(true);}
        if Instant::now()>=deadline {return Ok(false);}
        relabelings+=1;
        let h=h0.compose(g,d);
        debug_assert_eq!([a,b].map(|p|p.conjugate_by(h,d)),target);
        let h_key=if quotient_h {Some(colored_representation_key(d,[a,b,h]))} else {None};
        if h_key.as_ref().is_some_and(|key|completed_h.contains(key)) {
            equivalent_h+=1;next_h=current+1;next_third=0;return Ok(true);
        }
        let h_type=cycle_type(h,d);
        if !all_h && h_type!=vec![3;2*k] {next_h=current+1;next_third=0;return Ok(true);}
        // A commuting fourth entry has an h-invariant fixed-point set.
        if !all_types {
            let mut possible=1u64;
            for &length in &h_type {possible|=possible<<length;}
            if possible&(1u64<<(2*k))==0 {next_h=current+1;next_third=0;return Ok(true);}
        }
        eligible_h+=1;
        let cycles=permutation_cycles(d,h);
        // This partition is cheap and exact. Reject only a nontrivial equal-
        // sized block system preserved by all entries, never mere disconnection.
        let parts=point_orbits(d,&[a.compose(b,d),h]);
        let block_filter=parts.len()>1 && parts[0].len()>1
            && parts.iter().all(|p|p.len()==parts[0].len())
            && preserves_partition(d,a,&parts) && preserves_partition(d,b,&parts);
        let mut third_ordinal=if third_order>2 && current==start_h {start_third} else {0u64};
        let mut visit=|c:Perm| {
            let c_index=third_ordinal;third_ordinal+=1;
            if current==start_h && c_index<start_third {return Ok(true);}
            if Instant::now()>=deadline || thirds>=max_thirds {return Ok(false);}
            next_h=current;next_third=c_index+1;
            thirds+=1;
            if third_order>2 && (1..third_order).any(|q|third_order%q==0 && power_identity(d,c,q)) {return Ok(true);}
            let ambient=if a.is_even(d) && c.is_even(d) {Ambient::Alternating} else {Ambient::Symmetric};
            if target_group!="auto" && target_group!=ambient.short() {return Ok(true);}
            let t=[a,b,c,product_one_fourth(d,a,b,c)];
            if !all_types && cycle_type(t[3],d)!=exceptional {return Ok(true);}
            let fourth_index=index(t[3],d);
            if fourth_index>d-3 || 2*m+index(c,d)+fourth_index<2*d-2 {index_pruned+=1;return Ok(true);}
            matching+=1;
            if has_transitive_pair(d,&t) {direct_good+=1;return Ok(true);}
            if !generators_transitive(d,&t[..3]) {return Ok(true);}
            if block_filter && preserves_partition(d,c,&parts) {block_pruned+=1;return Ok(true);}
            let key=canonical(&t,d);
            if closed.contains(&key) {known_hits+=1;return Ok(true);}
            let generates=*group_cache.entry(key).or_insert_with(||tuple_generates_dynamic(d,ambient,&t));
            // A cache miss still uses the exact subgroup-order computation.
            // ponytail: bounded cache; clearing affects speed, never coverage.
            if group_cache.len()>500_000 {group_cache.clear();}
            if !generates {return Ok(true);}
            generating+=1;
            *genera.entry(genus_of(&t,d)).or_default()+=1;
            *passports.entry(passport_string(&t,d)).or_default()+=1;
            if !seen.insert(key) {return Ok(true);}
            // A positive walk witness is conclusive. Only witness-free starts
            // need the expensive exact quotient closure.
            if probe(d,t,&mut rng,64,true).1 {
                good+=1;
                writeln!(results,"weak_found\t0\t{}",tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
                return Ok(true);
            }
            writeln!(survivors,"group={}_{} seed={}",ambient.short(),d,tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
            survivors.sync_data().map_err(|e|e.to_string())?;
            let result=explore_goal(d,t,cap,deadline.min(Instant::now()+Duration::from_secs(10)),true);
            writeln!(results,"{}\t{}\t{}",result.status,result.states.len(),tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
            if result.weak_witness.is_some() {
                good+=1;seen.extend(result.states.iter().copied());
                if seen.len()>500_000 {seen.clear();}
            }
            else if result.status.starts_with("certified_") {
                let path=out.join(format!("CERTIFIED_{}{}_{}.txt",ambient.short(),d,bad));
                write_certificate(&path,d,ambient,&t,&result).map_err(|e|e.to_string())?;
                verify(&path)?;
                bad+=1;closed.extend(result.states);
                eprintln!("HEXAGON CERTIFICATE {}",path.display());
            } else {
                unresolved+=1;
                writeln!(hard,"group={}_{} status={} seed={}",ambient.short(),d,result.status,tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
            }
            Ok(true)
        };
        let mut skip=if current==start_h {start_third as u128} else {0};
        let exhausted=if third_order==2 {enumerate_commuting_involutions(&cycles,third_pairs,Perm::identity(d),&mut visit)?}
            else {enumerate_commuting_power(&cycles,0,third_order,Perm::identity(d),deadline,&mut skip,&mut HashMap::new(),&mut visit)?};
        if third_order>2 && exhausted && skip>0 {return Err("start-third exceeds profile's enumeration".into());}
        if current==start_h && start_third>third_ordinal {return Err("start-third exceeds profile's enumeration".into());}
        if exhausted {
            next_h=current+1;next_third=0;
            if let Some(key)=h_key {completed_h.insert(key);}
        }
        Ok(exhausted)
    })?;
    if start_h>ordinal {return Err("start-h exceeds profile's enumeration".into());}
    writeln!(summary,"quotient conjugate relabelings = {quotient_h}\nconjugate completed relabelings skipped = {equivalent_h}\nDeduplication covers conjugate (a,b,h) triples only after their entire third-entry enumeration completes. Raw third-visit counts can depend on chunk boundaries; coverage remains exact.").map_err(|e|e.to_string())?;
    writeln!(summary,"third entry order filter = {third_order}\ntarget ambient group = {target_group}\nindex obstructions pruned = {index_pruned}\npreserved block systems pruned = {block_pruned}\nknown failure states skipped = {known_hits}\naccepted generating seeds = {generating}\nsolved seeds = {good}\nhard unsolved candidates = {unresolved}").map_err(|e|e.to_string())?;
    writeln!(summary,"first involution transpositions = {m}").map_err(|e|e.to_string())?;
    if third_order==2 {writeln!(summary,"third involution transpositions = {third_pairs}").map_err(|e|e.to_string())?;}
    writeln!(summary,"direct-pair completions proven weak-good (generation not tested) = {direct_good}\ngenerating genus counts after pair pruning = {genera:?}\ngenerating passport counts after pair pruning = {passports:?}").map_err(|e|e.to_string())?;
    let group_label=if target_group!="auto" {target_group.as_str()} else if third_order==2 {ambient.short()}
        else if third_order%2==1 {if m%2==0 {"A"} else {"S"}} else {"auto"};
    let third_label=if third_order==2 {"compatible third involutions visited"} else {"compatible commuting third entries visited"};
    writeln!(summary,"construction = fixed involution pair graph\npair block sizes = {sizes:?}\npair odd path sizes = {paths:?}\npair balanced even path sizes = {even:?}\nfixed points per involution = {}\ngroup = {group_label}_{d}\nfourth type = {}\nall sigma1 relabelings = {all_h}\nrelabelings visited = {relabelings}\nrelabelings passing structural filters = {eligible_h}\n{third_label} = {thirds}\nfourth-type matches = {matching}\nfull-group seeds = {generating}\nweak-success starting states = {good}\ncertified weak-failure components = {bad}\nunresolved candidate starting states = {unresolved}\nsearch complete = {}\ntraversal exhausted = {complete}\nscan start h = {start_h}\nscan start third = {start_third}\nscan next h = {next_h}\nscan next third = {next_third}\nscan range resolved = {}\nelapsed seconds = {:.3}\nThis is a fixed pair-graph symmetry slice, not a degree census. Nonzero starts require gap-free verified prefix coverage.",paths.len()+2*even.len(),if all_types {"all".to_string()} else {format!("4^{k} 1^{}",2*k)},complete && unresolved==0 && start_h==0 && start_third==0,unresolved==0,start.elapsed().as_secs_f64()).map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())?;
    Ok(())
}

fn braid_fixed_seed(d: usize, ambient: Ambient, m: usize, genus: usize, exceptional: Option<&[usize]>,
                    rng: &mut SplitMix64, counts: &[AtomicU64;6]) -> Option<[Perm;4]> {
    let mut lengths = vec![1;d-2*m]; lengths.extend(vec![2;m]);
    let a = perm_of_type(&lengths,d,None);
    for _ in 0..200 {
        counts[0].fetch_add(1,Ordering::Relaxed);
        let b = perm_of_type(&lengths,d,Some(rng));
        let Some(h) = pair_braid_conjugator(d,a,b,rng) else { continue; };
        counts[1].fetch_add(1,Ordering::Relaxed);
        let Some(c) = commuting_involution(d,h,m,rng) else { continue; };
        counts[2].fetch_add(1,Ordering::Relaxed);
        let t = [a,b,c,product_one_fourth(d,a,b,c)];
        if 3*m+index(t[3],d) != 2*d-2+2*genus { continue; }
        if exceptional.is_some_and(|p|cycle_type(t[3],d)!=p) {continue;}
        counts[3].fetch_add(1,Ordering::Relaxed);
        // Avoid canonicalization on disconnected tuples; full group order is
        // needed anyway. No pair-product fingerprint is assumed necessary.
        if !tuple_generates_dynamic(d,ambient,&t) { continue; }
        counts[4].fetch_add(1,Ordering::Relaxed);
        if canonical(&t,d)==canonical(&hurwitz_move_dynamic(d,t,1),d) {
            counts[5].fetch_add(1,Ordering::Relaxed);
        }
        debug_assert_eq!(t.map(|g|g.conjugate_by(h,d)),hurwitz_move_dynamic(d,t,0));
        return Some(t);
    }
    None
}

// Two braid symmetries, not merely one. Let hkh=khk, r=h^-1 k^-1,
// [a,k]=1 and (ar)^3=1. Put b=h^-1 a h, c=k^-1 b k.
// Then abc=(hk)^3, which centralizes h,k. The braid relation gives [c,h]=1;
// hence (a,b,c,e)^h=sigma1(t), (a,b,c,e)^k=sigma2(t), e=(abc)^-1.
// We solve these equations by propagation, NOT by enumerating S_d or tuples.
// An additional a^order=1 restricts the desired entry types; it is NOT
// required by the derivation, so orders 3,4,... explore beyond involutions.
fn impose_commuting(k: Perm, order: usize, a: &mut [u8;MAX_D], inverse: &mut [u8;MAX_D], mut x: usize, mut y: usize) -> Option<bool> {
    let (root_x,root_y)=(x,y);
    let mut changed=false;
    loop {
        for (from,to) in [(x,y),(y,x)].into_iter().take(if order==2 {2} else {1}) {
            if inverse[to]!=u8::MAX && inverse[to]!=from as u8 {return None;}
            if a[from]==u8::MAX {a[from]=to as u8; inverse[to]=from as u8; changed=true;}
            else if a[from]!=to as u8 {return None;}
        }
        x=k.img[x] as usize; y=k.img[y] as usize;
        if x==root_x {break;}
    }
    if y!=root_y {return None;}
    Some(changed)
}

fn propagate_rigid(d: usize, k: Perm, r: Perm, order: usize, a: &mut [u8;MAX_D], inverse: &mut [u8;MAX_D]) -> bool {
    loop {
        let mut changed=false;
        for x in 0..d {
            if order>2 {
                let mut end=x;
                let mut ready=true;
                for _ in 1..order {
                    if a[end]==u8::MAX {ready=false; break;}
                    end=a[end] as usize;
                }
                if ready {
                    match impose_commuting(k,order,a,inverse,end,x) {
                        None=>return false,
                        Some(c)=>changed|=c,
                    }
                }
            }
            let y=a[r.img[x] as usize];
            if y==u8::MAX {continue;}
            let z=a[r.img[y as usize] as usize];
            if z==u8::MAX {continue;}
            // q(x)=y, q(y)=z forces q(z)=x for q=ar of order dividing 3.
            match impose_commuting(k,order,a,inverse,r.img[z as usize] as usize,x) {
                None=>return false,
                Some(c)=>changed|=c,
            }
        }
        if !changed {return true;}
    }
}

fn solve_rigid(d: usize, k: Perm, r: Perm, order: usize, lengths: &[usize;MAX_D], mut a: [u8;MAX_D],
    nodes: &mut u64, max_nodes: u64, deadline: Instant,
    visit: &mut impl FnMut(Perm)->Result<bool,String>) -> Result<bool,String> {
    if *nodes>=max_nodes || Instant::now()>=deadline {return Ok(false);}
    *nodes+=1;
    let mut inverse=[u8::MAX;MAX_D];
    for x in 0..d {if a[x]!=u8::MAX {inverse[a[x] as usize]=x as u8;}}
    if !propagate_rigid(d,k,r,order,&mut a,&mut inverse) {return Ok(true);}
    let root=(0..d).filter(|&x|a[x]==u8::MAX).min_by_key(|&x|
        (0..d).filter(|&y|inverse[y]==u8::MAX && lengths[y]==lengths[x]).count());
    let Some(x)=root else {
        let mut p=Perm::identity(d); p.img[..d].copy_from_slice(&a[..d]);
        return visit(p);
    };
    for y in 0..d {
        if inverse[y]!=u8::MAX || lengths[y]!=lengths[x] {continue;}
        let mut next=a;
        let mut next_inverse=inverse;
        if impose_commuting(k,order,&mut next,&mut next_inverse,x,y).is_some()
            && !solve_rigid(d,k,r,order,lengths,next,nodes,max_nodes,deadline,visit)? {return Ok(false);}
    }
    Ok(true)
}

fn rigid_branches(d:usize,k:Perm,r:Perm,order:usize,lengths:&[usize;MAX_D],start:u64,
    nodes:&mut u64,max_nodes:u64,deadline:Instant,
    visit:&mut impl FnMut(Perm)->Result<bool,String>)->Result<(bool,u64,u64),String> {
    let mut a=[u8::MAX;MAX_D];let mut inverse=[u8::MAX;MAX_D];
    if !propagate_rigid(d,k,r,order,&mut a,&mut inverse) {return Ok((true,0,0));}
    let root=(0..d).filter(|&x|a[x]==u8::MAX).min_by_key(|&x|
        (0..d).filter(|&y|inverse[y]==u8::MAX && lengths[y]==lengths[x]).count());
    let Some(x)=root else {
        if start>1 {return Err("rigid branch cursor out of range".into());}
        if start==1 {return Ok((true,1,1));}
        let mut p=Perm::identity(d);p.img[..d].copy_from_slice(&a[..d]);
        let complete=visit(p)?;return Ok((complete,u64::from(complete),1));
    };
    let targets:Vec<_>=(0..d).filter(|&y|inverse[y]==u8::MAX && lengths[y]==lengths[x]).collect();
    let total=targets.len() as u64;
    if start>total {return Err("rigid branch cursor out of range".into());}
    let mut next=start;
    for &y in targets.iter().skip(start as usize) {
        if Instant::now()>=deadline || *nodes>=max_nodes {return Ok((false,next,total));}
        let mut partial=a;let mut inv=inverse;
        if impose_commuting(k,order,&mut partial,&mut inv,x,y).is_some()
            && !solve_rigid(d,k,r,order,lengths,partial,nodes,max_nodes,deadline,visit)? {return Ok((false,next,total));}
        next+=1;
    }
    Ok((true,next,total))
}

fn rigid_tuple(d: usize, h: Perm, k: Perm, a: Perm) -> [Perm;4] {
    let b=a.conjugate_by(h.inverse(d),d);
    let c=b.conjugate_by(k.inverse(d),d);
    [a,b,c,product_one_fourth(d,a,b,c)]
}

fn commuting_entry_indices(d: usize, k: Perm, order: usize) -> u64 {
    // A t-cycle on k-cycles of length l, with cumulative phase of gcd g,
    // contributes index t*l-g and has order t*l/g. Pack these components by
    // dynamic programming; no permutation enumeration is needed for the bound.
    let mut counts=BTreeMap::new();
    for l in cycle_type(k,d) {*counts.entry(l).or_insert(0usize)+=1;}
    let mut indices=1u64;
    for (l,count) in counts {
        let mut dp=vec![0u64;count+1]; dp[0]=1;
        for used in 0..count {
            for t in 1..=count-used {
                if order!=0 && order%t!=0 {continue;}
                for g in 1..=l {
                    if l%g==0 && (order==0 || order%(t*l/g)==0) {dp[used+t]|=dp[used]<<(t*l-g);}
                }
            }
        }
        let mut combined=0;
        for i in 0..=d {if indices&(1<<i)!=0 {combined|=dp[count]<<i;}}
        indices=combined;
    }
    indices
}

fn braid_block_library(source_d: usize, h: Perm, k: Perm) -> BTreeMap<usize,Vec<(Perm,Perm)>> {
    let mut library:BTreeMap<usize,Vec<(Perm,Perm)>>=BTreeMap::new();
    let mut insert=|degree:usize,pair:(Perm,Perm)| {
        let encode=|(a,b):(Perm,Perm)|canonical(&[a,b,Perm::identity(degree),product_one_fourth(degree,a,b,Perm::identity(degree))],degree);
        let key=encode(pair);
        let entries=library.entry(degree).or_default();
        if !entries.iter().any(|&p|encode(p)==key) {entries.push(pair);}
    };
    for block in point_orbits(source_d,&[h,k]) {
        let n=block.len();
        let mut labels=[u8::MAX;MAX_D];
        for (i,&x) in block.iter().enumerate() {labels[x]=i as u8;}
        let restrict=|p:Perm| {
            let mut q=Perm::identity(n);
            for (i,&x) in block.iter().enumerate() {q.img[i]=labels[p.img[x] as usize];}
            q
        };
        let (hh,kk)=(restrict(h),restrict(k));
        let hk=hh.compose(kk,n);
        let center=hk.compose(hk,n).compose(hk,n);
        // A central permutation has equal cycle lengths on a transitive block.
        let order=cycle_type(center,n)[0];
        for step in 1..=order {
            if order%step!=0 {continue;}
            let power=(0..step).fold(Perm::identity(n),|p,_|p.compose(center,n));
            let fibers=permutation_cycles(n,power);
            let degree=fibers.len();
            let mut projection=[0;MAX_D];
            for (i,fiber) in fibers.iter().enumerate() {for &x in fiber {projection[x]=i as u8;}}
            let quotient=|p:Perm| {
                let mut q=Perm::identity(degree);
                for (i,fiber) in fibers.iter().enumerate() {
                    q.img[i]=projection[p.img[fiber[0]] as usize];
                    assert!(fiber.iter().all(|&x|projection[p.img[x] as usize]==q.img[i]));
                }
                q
            };
            insert(degree,(quotient(hh),quotient(kk)));
        }
    }
    // Independent, parameterized B3 representations on three cyclic fibers.
    // h swaps fibers 0,1 and advances fiber 2; k swaps 1,2 and advances 0.
    // Both hkh and khk swap 0,2 and advance EVERY fiber once. Thus braid
    // holds for every modulus m, and (hk)^3 advances every fiber twice.
    // This varies the symmetry group itself, not only representations of the
    // known A25 subgroup. It gives blocks of degrees 3m up to the exact limit.
    for m in 1..=MAX_D/3 {
        let n=3*m; let (mut h,mut k)=(Perm::identity(n),Perm::identity(n));
        for y in 0..m {
            h.img[y]=(m+y) as u8; h.img[m+y]=y as u8; h.img[2*m+y]=(2*m+(y+1)%m) as u8;
            k.img[y]=((y+1)%m) as u8; k.img[m+y]=(2*m+y) as u8; k.img[2*m+y]=(m+y) as u8;
        }
        insert(n,(h,k));
    }
    library
}

fn tetrahedral_braid_block(m: usize, phase: usize) -> (Perm,Perm) {
    // Central lift of x=(01)(23), y=(012) on four cyclic fibers.
    // x^2=y^3=z advances each fiber by 3; y advances its fixed fiber by 1.
    // h=y^-1 x, k=x^-1 y^2 then satisfy hk=y and hkh=khk=x.
    let n=4*m; assert!(m>0 && n<=MAX_D);
    let (mut x,mut y)=(Perm::identity(n),Perm::identity(n));
    for j in 0..m {
        x.img[j]=(m+(j+phase)%m) as u8; x.img[m+j]=((j+3+m-phase)%m) as u8;
        x.img[2*m+j]=(3*m+j) as u8; x.img[3*m+j]=(2*m+(j+3)%m) as u8;
        y.img[j]=(m+j) as u8; y.img[m+j]=(2*m+j) as u8;
        y.img[2*m+j]=((j+3)%m) as u8; y.img[3*m+j]=(3*m+(j+1)%m) as u8;
    }
    (y.inverse(n).compose(x,n),x.inverse(n).compose(y.compose(y,n),n))
}

fn octahedral_braid_block(m: usize) -> (Perm,Perm) {
    // S4 acts on its six two-subsets. Lift x=(01), y=(123) so that
    // x^2=y^3 advances each cyclic fiber by 3. Odd m coprime to 3 gives
    // a transitive lift. h=y^-1 x, k=x^-1 y^2 satisfy hkh=khk=x.
    assert!(m>0 && m%2==1 && m%3!=0 && 6*m<=MAX_D);
    let pairs:Vec<_>=(0..4).flat_map(|i|(i+1..4).map(move |j|[i,j])).collect();
    let induced=|p: [usize;4]| {
        let mut q=Perm::identity(6);
        for (i,pair) in pairs.iter().enumerate() {
            let mut target=[p[pair[0]],p[pair[1]]];target.sort();
            q.img[i]=pairs.iter().position(|&v|v==target).unwrap() as u8;
        }
        q
    };
    let lift=|base:Perm,power:usize| {
        let mut q=Perm::identity(6*m);
        for cycle in permutation_cycles(6,base) {
            let shift=if cycle.len()==1 {(0..m).find(|s|power*s%m==3%m).unwrap()} else {3%m};
            for (i,&x) in cycle.iter().enumerate() {
                for j in 0..m {
                    q.img[x*m+j]=(cycle[(i+1)%cycle.len()]*m+(j+if i+1==cycle.len() {shift} else {0})%m) as u8;
                }
            }
        }
        q
    };
    let x=lift(induced([1,0,2,3]),2);let y=lift(induced([0,2,3,1]),3);
    (y.inverse(6*m).compose(x,6*m),x.inverse(6*m).compose(y.compose(y,6*m),6*m))
}

// Independent modular representations, not a new order cutoff inside the
// exhausted A25 representation. Choose transitive x0^2=y0^3=1 on a base,
// then solve fiber-voltage equations x^2=y^3=z (one common cyclic shift).
// h=y^-1*x, k=x^-1*y^2 satisfy braid for ANY such x,y.
fn modular_braid_block(n:usize,m:usize,seed:u64)->Result<(Perm,Perm),String> {
    if m<2 || n%m!=0 || n>MAX_D || n/m<3 {return Err("modular label requires degree divisible by modulus >=2, base >=3".into());}
    let base=n/m; let mut rng=SplitMix64::new(seed);
    let mut possibilities=Vec::new();
    for pairs in 1..=base/2 {for triples in 1..=base/3 {
        if (base==2*pairs || m%2==1) && (base==3*triples || m%3!=0) {
            possibilities.push((pairs,triples));
        }
    }}
    if possibilities.is_empty() {return Err("no compatible fixed-point phases for this modulus/base".into());}
    for _ in 0..4096 {
        let (pairs,triples)=possibilities[rng.below(possibilities.len())];
        let mut xt=vec![1;base-2*pairs];xt.extend(vec![2;pairs]);
        let mut yt=vec![1;base-3*triples];yt.extend(vec![3;triples]);
        let x0=perm_of_type(&xt,base,None);let y0=perm_of_type(&yt,base,Some(&mut rng));
        if !generators_transitive(base,&[x0,y0]) {continue;}
        let mut lift=|p:Perm,power:usize| {
            let mut q=Perm::identity(n);
            for cycle in permutation_cycles(base,p) {
                let last_phase=(0..m).find(|s|(power/cycle.len())*s%m==1).unwrap();
                let mut phases=vec![0;cycle.len()];let mut sum=0;
                for phase in phases.iter_mut().take(cycle.len()-1) {*phase=rng.below(m);sum+=*phase;}
                phases[cycle.len()-1]=(last_phase+m-sum%m)%m;
                for (i,&x) in cycle.iter().enumerate() {for j in 0..m {
                    q.img[x*m+j]=(cycle[(i+1)%cycle.len()]*m+(j+phases[i])%m) as u8;
                }}
            }
            q
        };
        let x=lift(x0,2);let y=lift(y0,3);
        let h=y.inverse(n).compose(x,n);let k=x.inverse(n).compose(y.compose(y,n),n);
        if !generators_transitive(n,&[h,k]) {continue;}
        return Ok((h,k));
    }
    Err("no transitive modular lift found in bounded construction attempts".into())
}

fn glue_braid_blocks(source_d: usize, h: Perm, k: Perm, blocks: &[String]) -> Result<(usize,Perm,Perm),String> {
    let sizes=blocks.iter().map(|s| {
        let n=s.strip_prefix(['t','o','b']).unwrap_or(s).split(':').next().unwrap();
        n.parse::<usize>().map_err(|_|"invalid block degree")
    }).collect::<Result<Vec<_>,_>>()?;
    let d=sizes.iter().try_fold(0usize,|s,&n|s.checked_add(n)).ok_or("degree overflow")?;
    if !(4..=MAX_D).contains(&d) {return Err(format!("glued degree must be 4..{MAX_D}"));}
    let library=braid_block_library(source_d,h,k);
    let (mut hh,mut kk,mut offset)=(Perm::identity(d),Perm::identity(d),0);
    for (label,&size) in blocks.iter().zip(&sizes) {
        let (bh,bk)=if let Some(s)=label.strip_prefix('b') {
            let fields=s.split(':').collect::<Vec<_>>();
            if fields.len()!=3 {return Err("modular label requires bDEGREE:MODULUS:SEED".into());}
            let modulus=fields[1].parse().map_err(|_|"invalid modular modulus")?;
            let seed=fields[2].parse().map_err(|_|"invalid modular seed")?;
            modular_braid_block(size,modulus,seed)?
        } else if let Some(s)=label.strip_prefix('t') {
            let (_,phase)=s.split_once(':').ok_or("tetrahedral label requires tDEGREE:PHASE")?;
            let phase:usize=phase.parse().map_err(|_|"invalid tetrahedral phase")?;
            if size<4 || size%4!=0 || phase>=size/4 {return Err("invalid tetrahedral degree/phase".into());}
            tetrahedral_braid_block(size/4,phase)
        } else if label.starts_with('o') {
            if ![6,30].contains(&size) || label.contains(':') {return Err("octahedral label must be o6 or o30".into());}
            octahedral_braid_block(size/6)
        } else {
            if label.contains(':') {return Err("invalid standard block label".into());}
            let matching=library.get(&size).ok_or_else(||format!("no symmetry block of size {size}"))?;
            if matching.len()!=1 {return Err(format!("ambiguous symmetry block of size {size}"));}
            matching[0]
        };
        for i in 0..size {
            hh.img[offset+i]=(offset+bh.img[i] as usize) as u8;
            kk.img[offset+i]=(offset+bk.img[i] as usize) as u8;
        }
        offset+=size;
    }
    Ok((d,hh,kk))
}

fn colored_representation_key(d:usize,generators:[Perm;3])->String {
    let mut keys=Vec::new();
    for points in point_orbits(d,&generators) {
        let n=points.len();let mut labels=[0;MAX_D];
        for (i,&x) in points.iter().enumerate() {labels[x]=i as u8;}
        let restrict=|p:Perm| {
            let mut q=Perm::identity(n);
            for (i,&x) in points.iter().enumerate() {q.img[i]=labels[p.img[x] as usize];}q
        };
        let [a,b,c]=generators.map(restrict);
        let key=canonical(&[a,b,c,product_one_fourth(n,a,b,c)],n);
        keys.push(format!("{n}:{:032x}/{:032x}/{:032x}",key[0],key[1],key[2]));
    }
    keys.sort();keys.join(";")
}

fn braid_representation_key(d:usize,h:Perm,k:Perm)->String {
    colored_representation_key(d,[h,k,Perm::identity(d)])
}

fn braid_library(args:&[String])->Result<(),String> {
    let input=PathBuf::from(get_arg(args,"--certificate",None)?);verify(&input)?;
    let case=parse_certificate_case(&input)?;let d=case.d;
    let h=tuple_conjugator(d,case.seed,hurwitz_move_dynamic(d,case.seed,0))?;
    let k=tuple_conjugator(d,case.seed,hurwitz_move_dynamic(d,case.seed,1))?;
    let profiles=fs::read_to_string(get_arg(args,"--profiles",None)?).map_err(|e|e.to_string())?;
    let out=PathBuf::from(get_arg(args,"--out",None)?);fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let mut rows=OpenOptions::new().create_new(true).write(true).open(out.join("REPRESENTATIONS.tsv")).map_err(|e|e.to_string())?;
    writeln!(rows,"profile\tdegree\trepresentation_key").map_err(|e|e.to_string())?;
    for line in profiles.lines().filter(|s|!s.is_empty()) {
        let blocks=line.split(',').map(str::to_string).collect::<Vec<_>>();
        let (n,hh,kk)=glue_braid_blocks(d,h,k,&blocks)?;
        if hh.compose(kk,n).compose(hh,n)!=kk.compose(hh,n).compose(kk,n) {return Err("invalid constructed braid representation".into());}
        writeln!(rows,"{line}\t{n}\t{}",braid_representation_key(n,hh,kk)).map_err(|e|e.to_string())?;
    }
    rows.sync_all().map_err(|e|e.to_string())?;Ok(())
}

fn rigid_glue(args: &[String]) -> Result<(),String> {
    let external=optional_arg(args,"--braid-action");
    let coupled=has_flag(args,"--representation-blocks");
    if (external.is_some() || coupled) && (has_flag(args,"--certificate") || has_flag(args,"--blocks"))
        || (external.is_some() && coupled) {return Err("choose exactly one source: certificate/blocks, braid-action, or representation-blocks".into());}
    let (text,source_d,d,h,k,sizes)=if coupled {
        let (d,h,k,labels)=braid_frames::coupled(args)?;(String::new(),d,d,h,k,labels)
    } else if let Some(path)=&external {
        let (d,h,k)=braid_frames::read(Path::new(path))?;
        (String::new(),d,d,h,k,vec![format!("explicit:{path}")])
    } else {
        let input=PathBuf::from(get_arg(args,"--certificate",None)?);verify(&input)?;
        let text=fs::read_to_string(&input).map_err(|e|e.to_string())?;
        let value=|prefix|text.lines().find_map(|s|s.strip_prefix(prefix)).ok_or("missing certificate field");
        let source_d:usize=value("d=")?.parse().map_err(|_|"invalid source degree")?;
        let source=parse_tuple(value("seed=")?,source_d)?;
        let h=tuple_conjugator(source_d,source,hurwitz_move_dynamic(source_d,source,0))?;
        let k=tuple_conjugator(source_d,source,hurwitz_move_dynamic(source_d,source,1))?;
        let sizes:Vec<String>=get_arg(args,"--blocks",Some("9,16"))?.split(',').map(str::to_string).collect();
        let (d,h,k)=glue_braid_blocks(source_d,h,k,&sizes)?;(text,source_d,d,h,k,sizes)
    };
    if h.compose(k,d).compose(h,d)!=k.compose(h,d).compose(k,d) {return Err("source symmetries do not satisfy the braid relation".into());}
    let seconds:u64=parse_num(get_arg(args,"--seconds",Some("300"))?,"--seconds")?;
    let max_nodes:u64=parse_num(get_arg(args,"--max-nodes",Some("1000000000"))?,"--max-nodes")?;
    let cap:usize=parse_num(get_arg(args,"--cap",Some("64"))?,"--cap")?;
    let entry_order:usize=parse_num(get_arg(args,"--entry-order",Some("2"))?,"--entry-order")?;
    let diagnose=has_flag(args,"--diagnose");
    if seconds==0 || max_nodes==0 || cap==0 {return Err("positive limits required".into());}
    if entry_order!=0 && !(2..=6).contains(&entry_order) {return Err("entry-order must be 0 (unrestricted) or 2..6".into());}
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let create=|name|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|format!("fresh --out required: {e}"));
    let mut results=create("RESULTS.tsv")?;
    writeln!(results,"status\tquotient_states\tweak\tgenus\tseed").map_err(|e|e.to_string())?;
    let mut survivors=create("survivors.txt")?;
    let mut hard=create("hard_candidates.txt")?;
    let mut summary=create("SUMMARY.txt")?;
    writeln!(summary,"explicit braid action = {}\nrepresentation key = {}",external.as_deref().unwrap_or("disabled"),braid_representation_key(d,h,k)).map_err(|e|e.to_string())?;
    writeln!(summary,"entry order unrestricted = {}",entry_order==0).map_err(|e|e.to_string())?;
    if let Some(path)=optional_arg(args,"--exclude-completed-representations") {
        let text=fs::read_to_string(path).map_err(|e|e.to_string())?;
        if text.lines().next()!=Some("representation_key\tentry_order") {return Err("invalid completed-representation ledger".into());}
        let key=braid_representation_key(d,h,k);
        for line in text.lines().skip(1) {
            let (old,order)=line.split_once('\t').ok_or("malformed representation ledger")?;
            let order:usize=order.parse().map_err(|_|"invalid completed entry order")?;
            if old==key && (order==0 || (entry_order>0 && order%entry_order==0)) {
                writeln!(summary,"search mode = rigid-glue\nentry order divides = {entry_order}\nsearch complete = true\nrigid traversal exhausted = true\ncompleted representation skipped = true\naccepted generating seeds = 0\ncertified weak failures = 0\nhard unsolved candidates = 0").map_err(|e|e.to_string())?;
                return Ok(());
            }
        }
    }
    let mut lengths=[0;MAX_D];
    for cycle in permutation_cycles(d,k) {for &x in &cycle {lengths[x]=cycle.len();}}
    let hk=h.compose(k,d);
    let center=hk.compose(hk,d).compose(hk,d);
    let ambient=if center.is_even(d) {Ambient::Alternating} else {Ambient::Symmetric};
    if coupled && center.is_identity(d) {
        writeln!(summary,"search mode = rigid-glue\ngroup = {}_{}\nentry order divides = {entry_order}\nrigid traversal exhausted = true\nsearch complete = true\nstructural exclusion = identity exceptional entry; separate triple stratum\naccepted generating seeds = 0\ncertified weak failures = 0\nhard unsolved candidates = 0",ambient.short(),d).map_err(|e|e.to_string())?;
        return Ok(());
    }
    let indices=commuting_entry_indices(d,k,entry_order);
    let max_index=(0..d).filter(|i|i%2==usize::from(ambient==Ambient::Symmetric) && indices&(1<<i)!=0).max().unwrap_or(0);
    if 3*max_index+index(center,d)<2*d-2 {
        writeln!(summary,"search mode = rigid-glue\ngroup = {}_{}\nblocks = {sizes:?}\nentry order divides = {entry_order}\nsearch complete = true\nstructural exclusion = Riemann-Hurwitz index upper bound\nmax entry index = {max_index}\nexceptional index = {}\nconstraint nodes = 0\nalgebraic solutions = not enumerated\naccepted generating seeds = 0\ncertified weak failures = 0\nnew weak failure components relative to source = 0\nhard unsolved candidates = 0\nNo transitive tuple can satisfy this representation's constraints.",ambient.short(),d,index(center,d)).map_err(|e|e.to_string())?;
        summary.sync_all().map_err(|e|e.to_string())?;
        eprintln!("rigid index exclusion {}{d} blocks={sizes:?} max_index={max_index}",ambient.short());
        return Ok(());
    }
    let known:HashSet<Key>=text.lines().filter_map(|s|s.strip_prefix("STATE\t")).map(|s|parse_hex_key(&s.split('\t').collect::<Vec<_>>())).collect::<Result<_,_>>()?;
    let mut closed=known_weak_states(args,d)?;
    if coupled || external.is_some() || sizes.iter().any(|s|s.starts_with('b')) {
        let mut file=create("SYMMETRY.json")?;
        writeln!(file,"{{\"degree\":{d},\"h\":{:?},\"k\":{:?},\"center\":{:?}}}",h.img[..d].to_vec(),k.img[..d].to_vec(),center.img[..d].to_vec()).map_err(|e|e.to_string())?;
    }
    let (mut nodes,mut leaves,mut transitive,mut generating,mut weak,mut failures,mut new_failures,mut unresolved)=(0u64,0u64,0u64,0u64,0u64,0u64,0u64,0u64);
    let mut direct_pair_pruned=0u64;
    let mut genera=BTreeMap::new();
    let mut passports=BTreeMap::new();
    let mut partitions=BTreeMap::new();
    let mut proper_groups=BTreeMap::new();
    let mut obstruction=if diagnose {Some(create("OBSTRUCTIONS.tsv")?)} else {None};
    if let Some(file)=&mut obstruction {writeln!(file,"order\tblock_sizes\tblocks_zero_based\tseed").map_err(|e|e.to_string())?;}
    let start=Instant::now(); let deadline=start+Duration::from_secs(seconds);
    let mut last_progress=start;
    let branch_start:u64=parse_num(get_arg(args,"--start-branch",Some("0"))?,"--start-branch")?;
    eprintln!("rigid glue {}{d} blocks={sizes:?} entry_order={entry_order} seconds={seconds} max_nodes={max_nodes}",ambient.short());
    let mut visit=|a:Perm| {
        leaves+=1;
        let t=rigid_tuple(d,h,k,a);
        // Fail loudly if propagation or the algebraic construction is wrong.
        if !(0..entry_order).fold(Perm::identity(d),|p,_|p.compose(a,d)).is_identity(d) || a.compose(k,d)!=k.compose(a,d)
            || t[3]!=center.inverse(d) || t.map(|p|p.conjugate_by(h,d))!=hurwitz_move_dynamic(d,t,0)
            || t.map(|p|p.conjugate_by(k,d))!=hurwitz_move_dynamic(d,t,1) {
            return Err("rigid constraint solver produced an invalid tuple".into());
        }
        if last_progress.elapsed()>=Duration::from_secs(5) {
            eprintln!("rigid progress: algebraic={leaves} transitive={transitive} generating={generating} new_failures={new_failures} elapsed={:.0}s",start.elapsed().as_secs_f64());
            last_progress=Instant::now();
        }
        let orbits=point_orbits(d,&t[..3]);
        if diagnose {
            let mut sizes:Vec<_>=orbits.iter().map(Vec::len).collect();sizes.sort();
            *partitions.entry(sizes).or_insert(0u64)+=1;
        }
        if orbits.len()!=1 {return Ok(true);}
        transitive+=1;
        if has_transitive_pair(d,&t) {
            direct_pair_pruned+=1;
            writeln!(results,"weak_found\t1\ttrue\t{}\t{}",genus_of(&t,d),tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
            return Ok(true);
        }
        let order=subgroup_order(d,&t[..3]);
        if order!=ambient.target_order(d) {
            if diagnose {
                let blocks=invariant_blocks(d,&t[..3]);
                let sizes=blocks.as_ref().map(|b|b.iter().map(Vec::len).collect::<Vec<_>>());
                let count=proper_groups.entry((order,sizes.clone())).or_insert(0u64);
                if *count==0 {writeln!(obstruction.as_mut().unwrap(),"{order}\t{sizes:?}\t{blocks:?}\t{}",tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;}
                *count+=1;
            }
            return Ok(true);
        }
        generating+=1;
        let genus=genus_of(&t,d); *genera.entry(genus).or_insert(0u64)+=1;
        *passports.entry(passport_string(&t,d)).or_insert(0u64)+=1;
        let key=canonical(&t,d);
        if closed.contains(&key) {return Ok(true);}
        writeln!(survivors,"group={}_{} seed={}",ambient.short(),d,tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
        survivors.sync_data().map_err(|e|e.to_string())?;
        let result=explore_goal(d,t,cap,deadline.min(Instant::now()+Duration::from_secs(10)),true);
        writeln!(results,"{}\t{}\t{}\t{genus}\t{}",result.status,result.states.len(),result.weak_witness.is_some(),tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
        if result.weak_witness.is_some() {weak+=1;}
        else if result.status.starts_with("certified_") {
            let is_known=d==source_d && result.states.iter().all(|key|known.contains(key));
            let path=out.join(format!("CERTIFIED_{}{}_{}.txt",ambient.short(),d,failures));
            write_certificate(&path,d,ambient,&t,&result).map_err(|e|e.to_string())?;
            verify(&path)?;
            failures+=1; new_failures+=u64::from(!is_known);
            closed.extend(result.states);
            eprintln!("RIGID CERTIFICATE {} known_source_orbit={is_known}",path.display());
        } else {
            unresolved+=1;
            writeln!(hard,"group={}_{} status={} seed={}",ambient.short(),d,result.status,tuple_string_dynamic(d,&t)).map_err(|e|e.to_string())?;
        }
        Ok(Instant::now()<deadline)
    };
    let r=h.inverse(d).compose(k.inverse(d),d);
    let complete=if has_flag(args,"--start-branch") {
        let (complete,next,total)=rigid_branches(d,k,r,entry_order,&lengths,branch_start,&mut nodes,max_nodes,deadline,&mut visit)?;
        writeln!(summary,"rigid start branch = {branch_start}\nrigid next branch = {next}\nrigid total branches = {total}\nrigid traversal exhausted = {complete}").map_err(|e|e.to_string())?;
        complete
    } else {solve_rigid(d,k,r,entry_order,&lengths,[u8::MAX;MAX_D],&mut nodes,max_nodes,deadline,&mut visit)?};
    writeln!(summary,"search mode = rigid-glue\ngroup = {}_{}\nblocks = {:?}\nentry order divides = {entry_order}\nsearch complete = {}\nconstraint nodes = {nodes}\nalgebraic solutions = {leaves}\ntransitive solutions = {transitive}\naccepted generating seeds = {generating}\nweak witnesses (excluding cached failures) = {weak}\ncertified weak failures = {failures}\nnew weak failure components relative to source = {new_failures}\nhard unsolved candidates = {unresolved}\ngenerating genus counts = {genera:?}\nelapsed seconds = {:.3}\nCoverage is only the specified glued symmetry representation; not a degree census.",ambient.short(),d,sizes,complete && unresolved==0,start.elapsed().as_secs_f64()).map_err(|e|e.to_string())?;
    writeln!(summary,"generating passport counts = {passports:?}").map_err(|e|e.to_string())?;
    writeln!(summary,"direct-pair completions proven weak-good (generation not tested) = {direct_pair_pruned}").map_err(|e|e.to_string())?;
    if diagnose {writeln!(summary,"algebraic point-orbit partition counts = {partitions:?}\nproper transitive group order and block-size counts = {proper_groups:?}\nThese are per-solution obstructions in the enumerated representation, not a uniform theorem for untested representations.").map_err(|e|e.to_string())?;}
    summary.sync_all().map_err(|e|e.to_string())?;
    eprintln!("rigid complete={} nodes={nodes} algebraic={leaves} generating={generating} new_weak_components={new_failures}",complete && unresolved==0);
    Ok(())
}

fn conjugacy_class(d: usize, lengths: &[usize]) -> Vec<Perm> {
    // The least unused point anchors its cycle. Choosing its length from a
    // multiset and ordering the remaining points generates each permutation once.
    fn cycles(d:usize, unused:u64, counts:&mut [usize;MAX_D+1], p:Perm, out:&mut Vec<Perm>) {
        if unused==0 {out.push(p);return;}
        let root=unused.trailing_zeros() as usize;
        for length in 1..=d {
            if counts[length]==0 {continue;}
            counts[length]-=1;
            extend(d,unused&!(1<<root),counts,p,root,root,length-1,out);
            counts[length]+=1;
        }
    }
    fn extend(d:usize, unused:u64, counts:&mut [usize;MAX_D+1], mut p:Perm,
              root:usize, last:usize, remaining:usize, out:&mut Vec<Perm>) {
        if remaining==0 {p.img[last]=root as u8;cycles(d,unused,counts,p,out);return;}
        let mut choices=unused;
        while choices!=0 {
            let x=choices.trailing_zeros() as usize;choices&=choices-1;
            p.img[last]=x as u8;
            extend(d,unused&!(1<<x),counts,p,root,x,remaining-1,out);
        }
    }
    let mut counts=[0usize;MAX_D+1];
    for &length in lengths {counts[length]+=1;}
    let mut out=Vec::new();
    cycles(d,(1u64<<d)-1,&mut counts,Perm::identity(d),&mut out);
    // Preserve historical lexicographic scan indices for existing campaigns.
    out.sort_unstable_by_key(|p|p.img);
    out
}

fn small_perm_key(d:usize,p:Perm)->u64 {
    (0..d).fold(0,|key,i|key|((p.img[i] as u64)<<(4*i)))
}

fn centralizer_rows(d:usize,a:Perm,rows:Vec<Perm>)->Vec<Perm> {
    assert!(d<=11);
    let cycles=permutation_cycles(d,a);
    let mut generators=Vec::new();
    let mut previous=HashMap::<usize,Vec<usize>>::new();
    for cycle in cycles {
        if cycle.len()>1 {generators.push(cycle_perm(d,&cycle));}
        if let Some(old)=previous.insert(cycle.len(),cycle.clone()) {
            let mut swap=Perm::identity(d);
            for (&x,&y) in old.iter().zip(&cycle) {swap.img[x]=y as u8;swap.img[y]=x as u8;}
            generators.push(swap);
        }
    }
    let generators:Vec<_>=generators.into_iter().map(|g|(g,g.inverse(d))).collect();
    let mut seen=HashSet::with_capacity(rows.len());
    let mut representatives=Vec::new();
    let mut orbit=Vec::new();
    for b in rows {
        if !seen.insert(small_perm_key(d,b)) {continue;}
        representatives.push(b);orbit.clear();orbit.push(b);
        let mut head=0;
        while head<orbit.len() {
            let x=orbit[head];head+=1;
            for &(g,inverse) in &generators {
                let y=g.compose(x.compose(inverse,d),d);
                if seen.insert(small_perm_key(d,y)) {orbit.push(y);}
            }
        }
    }
    representatives
}

fn small_cycle_signature(d:usize,p:Perm)->u64 {
    let mut unused=(1u64<<d)-1;
    let mut signature=0;
    while unused!=0 {
        let mut x=unused.trailing_zeros() as usize;
        let mut length=0;
        while unused&(1<<x)!=0 {unused&=!(1<<x);length+=1;x=p.img[x] as usize;}
        signature+=1u64<<(4*length);
    }
    signature
}

fn read_parents(path: &Path, d: usize) -> Result<Vec<(usize,[Perm;4])>,String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut parents = Vec::new();
    let mut seen = HashSet::new();
    for line in text.lines().filter(|s| !s.trim().is_empty()) {
        let (header, tuple) = line.split_once(" seed=").ok_or("parent needs group=A_9 seed=[...] header")?;
        let (group, degree) = header.strip_prefix("group=").ok_or("invalid parent group")?.split_once('_').ok_or("invalid parent group/degree")?;
        let n: usize = degree.parse().map_err(|_| "invalid parent degree")?;
        if !(4..=MAX_D).contains(&n) { return Err("invalid parent degree".into()); }
        let t = parse_tuple(tuple,n)?;
        validate_tuple(n,Ambient::parse(group)?,&t)?;
        if n <= d && d <= n+2 && seen.insert((n,canonical(&t,n))) { parents.push((n,t)); }
    }
    if parents.is_empty() { return Err("no parents within zero, one, or two points of target degree".into()); }
    Ok(parents)
}

fn deform_seed(d: usize, ambient: Ambient, parents: &[(usize,[Perm;4])], rng: &mut SplitMix64) -> Option<[Perm;4]> {
    for _ in 0..200 {
        let (n, mut t) = parents[rng.below(parents.len())];
        // Explore different positions in the parent orbit before changing it.
        for _ in 0..rng.below(24) { t = hurwitz_move_dynamic(n,t,rng.below(6)); }
        for g in &mut t { for x in n..d { g.img[x] = x as u8; } }
        for x in n..d {
            let a = rng.below(n);
            let b = (a+1+rng.below(n-1))%n;
            let j = rng.below(3);
            t[j] = cycle_perm(d,&[a,b,x]).compose(t[j],d);
        }
        // Same-degree neighbours need a genuine change; extended seeds also
        // sometimes get one. This is a heuristic, not a braid equivalence claim.
        if n == d || rng.below(2) == 0 {
            let mut points = std::array::from_fn::<_, MAX_D, _>(|i| i);
            rng.shuffle(&mut points[..d]);
            let j = rng.below(3);
            t[j] = cycle_perm(d,&points[..3]).compose(t[j],d);
        }
        if ambient == Ambient::Symmetric && t[..3].iter().all(|g| g.is_even(d)) {
            let a = rng.below(d);
            let b = (a+1+rng.below(d-1))%d;
            let j = rng.below(3);
            t[j] = cycle_perm(d,&[a,b]).compose(t[j],d);
        }
        if ambient == Ambient::Alternating && t[..3].iter().any(|g| !g.is_even(d)) { continue; }
        t[3] = product_one_fourth(d,t[0],t[1],t[2]);
        if tuple_generates_dynamic(d,ambient,&t) { return Some(t); }
    }
    None
}

// Change one entry within its conjugacy class, keeping the other two fixed.
// This leaves the parent passport in three positions, but is not a braid move.
// The initial pair must be intransitive in a weak failure. An optional
// all-six-pair filter is a heuristic and does NOT establish exhaustive coverage.
fn local_switch_seed(d: usize, ambient: Ambient, states: &[[Perm;4]], switches: usize,
                     max_genus: i64, all_pairs: bool, preserve:bool, rng: &mut SplitMix64, counts: &[AtomicU64;6]) -> Option<[Perm;4]> {
    counts[0].fetch_add(1,Ordering::Relaxed);
    let mut t=states[rng.below(states.len())];
    let original_fourth=if preserve {Some(cycle_type(t[3],d))} else {None};
    let j=rng.below(3);
    let mut r=Perm::identity(d);
    for _ in 0..switches {
        let x=rng.below(d);
        let y=(x+1+rng.below(d-1))%d;
        let mut transposition=Perm::identity(d);
        transposition.img[x]=y as u8;
        transposition.img[y]=x as u8;
        r=transposition.compose(r,d);
    }
    let changed=t[j].conjugate_by(r,d);
    if changed==t[j] {return None;}
    t[j]=changed;
    t[3]=product_one_fourth(d,t[0],t[1],t[2]);
    if original_fourth.is_some_and(|p|p!=cycle_type(t[3],d)) {return None;}
    counts[1].fetch_add(1,Ordering::Relaxed);
    if genus_of(&t,d)>max_genus {return None;}
    counts[2].fetch_add(1,Ordering::Relaxed);
    if pair_transitive(d,t[0],t[1]) || (all_pairs && (0..4).any(|a|(a+1..4).any(|b|pair_transitive(d,t[a],t[b])))) {return None;}
    counts[3].fetch_add(1,Ordering::Relaxed);
    if point_orbits(d,&t[..3]).len()!=1 {return None;}
    counts[4].fetch_add(1,Ordering::Relaxed);
    if subgroup_order(d,&t[..3])!=ambient.target_order(d) {return None;}
    counts[5].fetch_add(1,Ordering::Relaxed);
    Some(t)
}

fn probe(d: usize, t: [Perm;4], rng: &mut SplitMix64, steps: usize, weak_only: bool) -> (bool,bool,usize) {
    probe_with_word(d,t,rng,steps,weak_only,None)
}

fn probe_with_word(d: usize, mut t: [Perm;4], rng: &mut SplitMix64, steps: usize,
                   weak_only: bool, mut word: Option<&mut Vec<u8>>) -> (bool,bool,usize) {
    if let Some(w)=word.as_deref_mut() {w.clear();}
    let mut weak=false;
    let mut pair_changed=true;
    for step in 0..=steps {
        if pair_changed && pair_transitive(d,t[0],t[1]) {
            weak=true;
            if weak_only || subgroup_order(d,&t[..2])>=factorial_u128(d)/2 {return (true,true,step);}
        }
        if step==steps {break;}
        let mv=rng.below(6);
        if let Some(w)=word.as_deref_mut() {w.push(mv as u8);}
        t=hurwitz_move_dynamic(d,t,mv);
        // Moves on entries (1,2) preserve their subgroup; moves on (3,4)
        // leave it unchanged. Only the middle braid generator changes it.
        pair_changed=mv%3==1;
    }
    (false,weak,steps+1)
}

// A partial search proves connectivity when it hits a previously visited
// canonical state. Disjoint partial searches do NOT prove distinct orbits.
struct OrbitMapResult {
    status: &'static str,
    states: Vec<Key>,
    weak: bool,
    overlap: Option<usize>,
}

fn map_orbit(d: usize, seed: [Perm; 4], cap: usize, deadline: Instant, known: &HashMap<Key, usize>) -> OrbitMapResult {
    let root = canonical(&seed, d);
    let mut out = OrbitMapResult { status: "capped", states: vec![root], weak: false, overlap: None };
    let mut seen = HashSet::from([root]);
    let mut head = 0;
    while head < out.states.len() {
        if Instant::now() >= deadline { out.status = "timeout"; return out; }
        let t = unpack(out.states[head], d);
        head += 1;
        out.weak |= pair_transitive(d, t[0], t[1]);
        for mv in 0..3 {
            let next = canonical(&hurwitz_move_dynamic(d, t, mv), d);
            if seen.contains(&next) { continue; }
            if let Some(&id) = known.get(&next) {
                out.status = "overlap";
                out.overlap = Some(id);
                return out;
            }
            if seen.len() >= cap { return out; }
            seen.insert(next);
            out.states.push(next);
        }
    }
    out.status = if out.weak { "closed_weak_success" } else { "closed_weak_failure" };
    out
}

fn orbit_map(args: &[String], d: usize, ambient: Ambient, cap: usize, seconds: u64) -> Result<(), String> {
    let passport_text = get_arg(args, "--passport", None)?;
    let passport = parse_passport(&passport_text, d, ambient)?;
    let master: u64 = parse_num(get_arg(args, "--seed", Some("7640891576956012809"))?, "--seed")?;
    let trials: u64 = parse_num(get_arg(args, "--seeds", Some("1000000000"))?, "--seeds")?;
    let total_cap: usize = parse_num(get_arg(args, "--max-states", Some("250000"))?, "--max-states")?;
    if total_cap == 0 { return Err("max-states must be positive".into()); }
    let mut replay = Vec::new();
    if let Ok(path) = get_arg(args, "--replay", None) {
        for (line_no, line) in fs::read_to_string(path).map_err(|e| e.to_string())?.lines().enumerate() {
            if line.trim().is_empty() { continue; }
            let text = line.split_once("seed=").map(|(_, s)| s).unwrap_or(line);
            let seed = parse_tuple(text, d).map_err(|e| format!("replay line {}: {e}", line_no + 1))?;
            validate_tuple(d, ambient, &seed).map_err(|e| format!("replay line {}: {e}", line_no + 1))?;
            if seed.iter().enumerate().any(|(i, &p)| cycle_type(p, d) != passport[i]) {
                return Err(format!("replay line {} has the wrong passport", line_no + 1));
            }
            replay.push(seed);
        }
    }
    let out_dir = PathBuf::from(get_arg(args, "--out", None)?);
    fs::create_dir_all(&out_dir).map_err(|e| e.to_string())?;
    let mut rows = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out_dir.join("ORBITS.tsv")).map_err(|e| format!("use a fresh --out: {e}"))?);
    writeln!(rows, "sample\tstatus\tquotient_states_seen\tweak_seen\tconnected_to_sample\tseed").map_err(|e| e.to_string())?;
    let start = Instant::now();
    let deadline = start + Duration::from_secs(seconds);
    let mut rng = SplitMix64::new(master);
    let mut known: HashMap<Key, usize> = HashMap::new();
    let mut attempted = 0u64;
    let mut accepted = 0usize;
    let mut overlaps = 0usize;
    let mut closed_good = 0usize;
    let mut closed_bad = 0usize;
    let mut capped = 0usize;
    let mut weak_seen = 0usize;
    let limit = trials.saturating_add(replay.len() as u64);
    while attempted < limit && known.len() < total_cap && Instant::now() < deadline {
        let trial = attempted;
        attempted += 1;
        let seed = if trial < replay.len() as u64 { Some(replay[trial as usize]) }
                   else { passport_seed(d, ambient, &passport, false, &mut rng) };
        let Some(seed) = seed else { continue; };
        accepted += 1;
        let root = canonical(&seed, d);
        let result = if let Some(&id) = known.get(&root) {
            OrbitMapResult { status: "overlap", states: Vec::new(), weak: false, overlap: Some(id) }
        } else {
            map_orbit(d, seed, cap.min(total_cap - known.len()), deadline, &known)
        };
        let id = result.overlap.unwrap_or(accepted);
        for &key in &result.states { known.insert(key, id); }
        if result.weak { weak_seen += 1; }
        match result.status {
            "overlap" => overlaps += 1,
            "closed_weak_success" => closed_good += 1,
            "closed_weak_failure" => {
                closed_bad += 1;
                let mut orders = BTreeMap::new();
                for &key in &result.states {
                    let t = unpack(key, d);
                    *orders.entry(subgroup_order(d, &t[..2])).or_insert(0) += 1;
                }
                let proof = Exploration { status: "certified_weak_and_strong_failure", weak_witness: None, weak_word: None, states: result.states.clone(), orders };
                let representative = *result.states.iter().min().unwrap();
                let path = out_dir.join(format!("CERTIFIED_{}{}_{:032x}_{:032x}_{:032x}.txt", ambient.short(), d, representative[0], representative[1], representative[2]));
                write_certificate(&path, d, ambient, &seed, &proof).map_err(|e| e.to_string())?;
                eprintln!("CLOSED WEAK FAILURE: {}", path.display());
            }
            _ => capped += 1,
        }
        writeln!(rows, "{accepted}\t{}\t{}\t{}\t{}\t{}", result.status, result.states.len(), result.weak, result.overlap.map_or(String::new(), |x| x.to_string()), tuple_string_dynamic(d, &seed)).map_err(|e| e.to_string())?;
        rows.flush().map_err(|e| e.to_string())?;
        if accepted % 100 == 0 {
            eprintln!("orbit map: accepted={accepted} states={} proven_overlaps={overlaps} closed_good={closed_good} closed_bad={closed_bad} elapsed={:.0}s", known.len(), start.elapsed().as_secs_f64());
        }
    }
    let mut summary = OpenOptions::new().write(true).create_new(true).open(out_dir.join("SUMMARY.txt")).map_err(|e| e.to_string())?;
    writeln!(summary, "search mode = orbit-map\npassport = {passport_text}\ngroup = {}_{}\nfamily = passport\nweak only = true\nmaster seed = {master}\nper-seed quotient-state cap = {cap}\nglobal quotient-state cap = {total_cap}\nattempted seeds = {attempted}\naccepted generating seeds = {accepted}\nsolved seeds = {weak_seen}\nproven overlapping samples = {overlaps}\nclosed weak-success components = {closed_good}\ncertified weak failures = {closed_bad}\nhard unsolved candidates = {capped}\nglobal stored quotient states = {}\nelapsed seconds = {:.3}\nOnly overlap and closed-component claims are proven. Disjoint partial searches are not distinct orbits. A weak_seen count excludes samples skipped by proven overlap.", ambient.short(), d, known.len(), start.elapsed().as_secs_f64()).map_err(|e| e.to_string())?;
    summary.sync_all().map_err(|e| e.to_string())?;
    Ok(())
}

fn smart_hunt(args: &[String], d: usize, ambient: Ambient, cap: usize, seconds: u64) -> Result<(),String> {
    if has_flag(args,"--joint-benchmark") {return weak_discovery::benchmark_joint(args,d,ambient);}
    if has_flag(args,"--audit-construction") {return weak_discovery::audit(args,d,ambient);}
    let threads:usize=parse_num(get_arg(args,"--threads",Some("8"))?,"--threads")?;
    let seeds:u64=parse_num(get_arg(args,"--seeds",Some("1000000000"))?,"--seeds")?;
    let seed_seconds:u64=parse_num(get_arg(args,"--candidate-seconds",Some("10"))?,"--candidate-seconds")?;
    let master:u64=parse_num(get_arg(args,"--seed",Some("7640891576956012809"))?,"--seed")?;
    let steps:usize=parse_num(get_arg(args,"--steps",Some("64"))?,"--steps")?;
    let exact=args[1]=="passport-exact";
    let family=if exact {"passport".to_string()} else {get_arg(args,"--family",Some("low-index"))?};
    let weak_only=has_flag(args,"--weak-only");
    let finite=matches!(family.as_str(),"s3-exact"|"pair-fiber-exact");
    if finite && (!weak_only || has_flag(args,"--no-braid-fixed")) {return Err("finite construction requires weak-only; no-braid-fixed is not an enumeration filter".into());}
    let mut symmetry=if matches!(family.as_str(),"square-symmetry"|"symmetry-centralizer"|"triangle-symmetry"|"pure-symmetry"|"pure-centralizer"|"marked-symmetry"|"dessin-glue"|"pair-glue"|"resistant-glue"|"cut-join-glue"|"pair-fiber"|"independent-symmetry"|"independent-glue"|"double-coset-glue"|"coordinated-glue"|"coordinated-fiber"|"joint-pure") {
        if !weak_only || (!matches!(family.as_str(),"dessin-glue"|"pair-glue"|"resistant-glue"|"cut-join-glue"|"pair-fiber"|"independent-glue"|"double-coset-glue"|"coordinated-glue"|"coordinated-fiber") && has_flag(args,"--no-braid-fixed")) {return Err("symmetry families require weak-only and retain fixed starts".into());}
        Some(weak_discovery::Factory::new(args,d,ambient,&family)?)
    } else {None};
    if threads==0 || threads>128 || seed_seconds==0 || steps>10000 {return Err("require 1..128 threads, positive candidate-seconds, steps <= 10000".into());}
    let passport = if family=="passport" {Some(parse_passport(&get_arg(args,"--passport",None)?,d,ambient)?)} else {None};
    let structural = if family=="braid-fixed" {
        let m:usize=parse_num(get_arg(args,"--involution-pairs",None)?,"--involution-pairs")?;
        let g:usize=parse_num(get_arg(args,"--genus",Some("0"))?,"--genus")?;
        if m==0 || m>d/2 || m%2!=usize::from(ambient==Ambient::Symmetric) || g>d
            || 3*m>2*d-2+2*g || 2*d-2+2*g-3*m>d-1 {
            return Err("infeasible equal-involution m/genus/group combination".into());
        }
        Some((m,g))
    } else {None};
    let exceptional=if has_flag(args,"--exceptional-type") {
        let Some((m,g))=structural else {return Err("exceptional-type requires braid-fixed".into());};
        let mut p=get_arg(args,"--exceptional-type",None)?.split('.').map(|s|s.parse::<usize>().map_err(|_|"invalid exceptional cycle length")).collect::<Result<Vec<_>,_>>()?;
        if p.iter().any(|&l|l==0 || l>d) || p.iter().sum::<usize>()>d {return Err("invalid exceptional cycle type".into());}
        p.extend(vec![1;d-p.iter().sum::<usize>()]);p.sort();
        if p.iter().map(|l|l-1).sum::<usize>()!=2*d-2+2*g-3*m {return Err("exceptional type disagrees with genus".into());}
        Some(p)
    } else {None};
    let structural_counts: [AtomicU64;6] = std::array::from_fn(|_|AtomicU64::new(0));
    let local_counts: [AtomicU64;6] = std::array::from_fn(|_|AtomicU64::new(0));
    let all_pairs=has_flag(args,"--all-pairs-intransitive");
    let no_fixed=has_flag(args,"--no-braid-fixed");
    let preserve=has_flag(args,"--preserve-passport");
    if (no_fixed || preserve) && (!weak_only || exact) {return Err("diversity filters require non-exact weak-only search".into());}
    if preserve && family!="local-switch" {return Err("preserve-passport requires local-switch".into());}
    if exact && (d>11 || has_flag(args,"--replay") || has_flag(args,"--replay-only")) {
        return Err("passport-exact supports d<=11 and does not accept replay inputs".into());
    }
    let mut pair_pruned=0u64;
    let mut conjugate_pruned=0u64;
    let quotient_rows=has_flag(args,"--centralizer-rows");
    if quotient_rows && !exact {return Err("centralizer-rows requires passport-exact".into());}
    let classes=if exact {
        let p=passport.as_ref().unwrap();
        let mut b=conjugacy_class(d,&p[1]);
        let mut c=if p[1]==p[2] {b.clone()} else {conjugacy_class(d,&p[2])};
        let a=perm_of_type(&p[0],d,None);
        if weak_only {
            let raw=(b.len() as u64)*(c.len() as u64);
            b.retain(|&x|!pair_transitive(d,a,x));
            c.retain(|&x|!pair_transitive(d,a,x));
            pair_pruned=raw-(b.len() as u64)*(c.len() as u64);
        }
        if quotient_rows {
            let before=b.len();
            b=centralizer_rows(d,a,b);
            conjugate_pruned=((before-b.len()) as u64)*(c.len() as u64);
        }
        Some((a,b,c))
    } else {None};
    // A transitive first pair is already a weak witness for every generating
    // completion, so those full rows are proven good without scanning them.
    let scan_total=classes.as_ref().map(|(_,b,c)| b.len() as u64*c.len() as u64).unwrap_or(0);
    let fourth_signature=if exact {passport.as_ref().unwrap()[3].iter().fold(0u64,|s,&l|s+(1u64<<(4*l)))} else {0};
    let scan_start:u64=parse_num(get_arg(args,"--start",Some("0"))?,"--start")?;
    let scan_count:u64=parse_num(get_arg(args,"--count",Some("18446744073709551615"))?,"--count")?;
    if exact && (scan_count==0 || (scan_total==0 && scan_start>0) || (scan_total>0 && scan_start>=scan_total)) { return Err("empty or out-of-range passport scan".into()); }
    let parents = if family=="deform" {read_parents(Path::new(&get_arg(args,"--parents",None)?),d)?} else {Vec::new()};
    let local = if family=="local-switch" {
        let path=PathBuf::from(get_arg(args,"--certificate",None)?);
        verify(&path)?;
        let source=parse_certificate_case(&path)?;
        if source.d!=d || source.ambient!=ambient || !weak_only || source.states.is_none() {
            return Err("local-switch requires a matching weak-only certificate".into());
        }
        let switches:usize=parse_num(get_arg(args,"--switches",Some("1"))?,"--switches")?;
        let max_genus:i64=parse_num(get_arg(args,"--max-genus",Some("100"))?,"--max-genus")?;
        if switches==0 || switches>64 || max_genus<0 {return Err("local-switch requires 1..64 switches and nonnegative max-genus".into());}
        let states=source.states.unwrap().into_iter().map(|key|unpack(key,d)).collect::<Vec<_>>();
        Some((path,states,switches,max_genus))
    } else {None};
    let family_parsed=if finite || symmetry.is_some() || matches!(family.as_str(),"low-index"|"passport"|"deform"|"braid-fixed"|"local-switch") {None} else {Some(SeedFamily::parse(&family)?)};
    if let Some(f)=family_parsed {
        let valid=match f {
            SeedFamily::Imprimitive => !nontrivial_divisors(d).is_empty(),
            SeedFamily::Affine => is_prime_usize(d) && d>=5,
            SeedFamily::Cyclic => ambient==Ambient::Symmetric || d%2==1,
            _ => true,
        };
        if !valid {return Err("family is unavailable at this degree/group".into());}
    }
    let mut replay=Vec::new();
    let mut replay_keys=HashSet::new();
    if let Ok(path)=get_arg(args,"--replay",None) {
        for (line_no,line) in fs::read_to_string(path).map_err(|e|e.to_string())?.lines().enumerate() {
            if line.trim().is_empty() {continue;}
            let text=line.split_once("seed=").map(|(_,s)|s).unwrap_or(line);
            let t=parse_tuple(text,d).map_err(|e|format!("replay line {}: {e}",line_no+1))?;
            validate_tuple(d,ambient,&t).map_err(|e|format!("replay line {}: {e}",line_no+1))?;
            if replay_keys.insert(canonical(&t,d)) {replay.push(t);}
        }
    }
    let out_dir=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out_dir).map_err(|e|e.to_string())?;
    // A fresh directory/attempt is required. Resume is at completed-job level
    // in the overnight runner; no ledger or partial result is truncated.
    let create=|name:&str| OpenOptions::new().write(true).create_new(true).open(out_dir.join(name));
    if finite {for t in if family=="s3-exact" {weak_discovery::s3_population(args,d,ambient,&out_dir)?}
        else {weak_discovery::fiber_population(args,d,ambient,&out_dir)?} {
        if replay_keys.insert(canonical(&t,d)) {replay.push(t);}
    }}
    let mut results=create("RESULTS.tsv").map_err(|e|format!("use a fresh --out: {e}"))?;
    writeln!(results,"status\tquotient_states\tweak\tseed").map_err(|e|e.to_string())?;
    let results=Mutex::new(results);
    let candidates=Mutex::new(create("hard_candidates.txt").map_err(|e|e.to_string())?);
    let pending=Mutex::new(create("survivors.txt").map_err(|e|e.to_string())?);
    let mut hints=create("PROBE_RECORDS.tsv").map_err(|e|e.to_string())?;
    writeln!(hints,"worker\tgoal\tfirst_witness_walk_step\tseed").map_err(|e|e.to_string())?;
    let hints=Mutex::new(hints);
    let histogram=Mutex::new(vec![0u64;steps+2]);
    let done=std::sync::atomic::AtomicBool::new(false);
    let claimed=AtomicU64::new(0);
    let attempted=AtomicU64::new(0);
    let generated=AtomicU64::new(0);
    let solved=AtomicU64::new(0);
    let weak=AtomicU64::new(0);
    let hard=AtomicU64::new(0);
    let certified=AtomicU64::new(0);
    let weak_fail=AtomicU64::new(0);
    let direct_pair_pruned=AtomicU64::new(0);
    let short_words_enabled=exact && weak_only && has_flag(args,"--short-weak-witnesses");
    let short_word_pruned=AtomicU64::new(0);
    let repeats=AtomicU64::new(0);
    let promoted_repeats=AtomicU64::new(0);
    let known=known_weak_states(args,d)?;
    if let Some(factory)=&mut symmetry {factory.set_known(&known);}
    let closed_nonempty=std::sync::atomic::AtomicBool::new(!known.is_empty());
    let closed=Mutex::new(known);
    let promoted=Mutex::new(HashSet::<Key>::new());
    let positive=if weak_only && has_flag(args,"--cache-positive-witnesses") {Some(std::sync::Arc::new(positive_witness::Cache::new(args,d,&out_dir)?))} else {None};
    if let Some(factory)=&mut symmetry {factory.set_positive(positive.clone())?;}
    let probe_witness_limit:u64=parse_num(get_arg(args,"--probe-witness-limit",Some("0"))?,"--probe-witness-limit")?;
    if probe_witness_limit>1024 || (probe_witness_limit>0 && (positive.is_none() || steps>4096)) {
        return Err("probe-witness-limit requires weak positive caching, steps<=4096 and limit<=1024".into());
    }
    let probe_witness_exports=AtomicU64::new(0);
    let start=Instant::now();
    let deadline=start+Duration::from_secs(seconds);
    let replay_only=finite || has_flag(args,"--replay-only");
    let work_count=if exact {scan_count.min(scan_total-scan_start)} else if replay_only {replay.len() as u64}
        else {(replay.len() as u64).saturating_add(if symmetry.as_ref().is_some_and(|f|f.is_empty()) {0} else {seeds})};
    eprintln!("smart hunt {}{} family={family} weak_only={weak_only} replay={} threads={threads} seconds={seconds} cap={cap}",ambient.short(),d,replay.len());
    thread::scope(|scope| -> Result<(),String> {
        let monitor=scope.spawn(|| {
            while !done.load(Ordering::Relaxed) {
                thread::park_timeout(Duration::from_secs(5));
                if done.load(Ordering::Relaxed) {break;}
                eprintln!("smart progress: generating={} solved={} hard={} certified={} repeat_bad_seeds={} elapsed={:.0}s",generated.load(Ordering::Relaxed),solved.load(Ordering::Relaxed),hard.load(Ordering::Relaxed),certified.load(Ordering::Relaxed),repeats.load(Ordering::Relaxed),start.elapsed().as_secs_f64());
            }
        });
        let mut workers=Vec::new();
        for tid in 0..threads {
            let replay=&replay; let closed=&closed; let out_dir=&out_dir;
            let promoted=&promoted;
            let passport=&passport; let parents=&parents; let hints=&hints; let histogram=&histogram;
            let classes=&classes;
            let structural_counts=&structural_counts;
            let local_counts=&local_counts;
            let local=&local;
            let symmetry=&symmetry;
            let exceptional=&exceptional;
            let claimed=&claimed; let attempted=&attempted; let generated=&generated; let solved=&solved;
            let weak=&weak; let hard=&hard; let certified=&certified; let weak_fail=&weak_fail; let repeats=&repeats;
            let promoted_repeats=&promoted_repeats;
            let direct_pair_pruned=&direct_pair_pruned; let closed_nonempty=&closed_nonempty;
            let short_word_pruned=&short_word_pruned;
            let results=&results; let candidates=&candidates; let pending=&pending;
            let positive=&positive;
            let probe_witness_exports=&probe_witness_exports;
            workers.push(scope.spawn(move || -> io::Result<()> {
                let mut rng=SplitMix64::new(master ^ (tid as u64+1).wrapping_mul(0xD1B54A32D192ED03));
                let mut walk_counts=vec![0u64;steps+2];
                let mut record=0;
                let mut probe_word=Vec::new();
                let (mut next_job,mut block_end)=(0u64,0u64);
                loop {
                    if next_job==block_end {
                        if Instant::now()>=deadline {break;}
                        next_job=claimed.fetch_add(if exact {128} else {1},Ordering::Relaxed);
                        block_end=next_job.saturating_add(if exact {128} else {1}).min(work_count);
                        if next_job>=work_count {break;}
                    }
                    // Finish every claimed block, including at the deadline:
                    // scan next remains a gap-free prefix across all workers.
                    let job=next_job;next_job+=1;
                    if job>=work_count {break;}
                    attempted.fetch_add(1,Ordering::Relaxed);
                    let replayed=job<replay.len() as u64;
                    let seed=if let Some((a,bs,cs))=classes {
                        let index=scan_start+job;
                        let b=bs[(index/cs.len() as u64) as usize];
                        let c=cs[(index%cs.len() as u64) as usize];
                        let product=a.compose(b,d).compose(c,d);
                        // Inversion preserves cycle type; construct g4 only on matches.
                        if small_cycle_signature(d,product)!=fourth_signature {continue;}
                        let t=[*a,b,c,product.inverse(d)];
                        // Any chosen pair can be moved to the first two slots
                        // by Hurwitz moves, up to simultaneous conjugacy.
                        if weak_only && has_transitive_pair(d,&t) {
                            direct_pair_pruned.fetch_add(1,Ordering::Relaxed);
                            continue;
                        }
                        if short_words_enabled && positive_witness::short_witness(d,t).is_some() {
                            short_word_pruned.fetch_add(1,Ordering::Relaxed);
                            continue;
                        }
                        if !tuple_generates_dynamic(d,ambient,&t) {continue;}
                        t
                    } else if replayed {replay[job as usize]} else {
                        let t=match family_parsed {
                            Some(f)=>make_seed(d,ambient,f,&mut rng),
                            None=>if let Some(factory)=symmetry {factory.generate(d,ambient,&mut rng)}
                                  else if let Some((m,g))=structural {braid_fixed_seed(d,ambient,m,g,exceptional.as_deref(),&mut rng,structural_counts)}
                                  else if let Some(p)=passport {passport_seed(d,ambient,p,weak_only,&mut rng)}
                                  else if let Some((_,states,switches,max_genus))=local {local_switch_seed(d,ambient,states,*switches,*max_genus,all_pairs,preserve,&mut rng,local_counts)}
                                  else if !parents.is_empty() {deform_seed(d,ambient,parents,&mut rng)}
                                  else {low_index_seed(d,ambient,weak_only,&mut rng)},
                        };
                        match t {Some(t)=>t,None=>continue}
                    };
                    if !replayed && weak_only && all_pairs && has_transitive_pair(d,&seed) {continue;}
                    if !replayed && no_fixed {
                        let key=canonical(&seed,d);
                        if (0..3).any(|mv|canonical(&hurwitz_move_dynamic(d,seed,mv),d)==key) {continue;}
                    }
                    generated.fetch_add(1,Ordering::Relaxed);
                    // A census revisits a small quotient many times under the
                    // centralizer of g1. Skip certified components before walks.
                    let exact_key=if (exact || symmetry.is_some()) && closed_nonempty.load(Ordering::Acquire) {Some(canonical(&seed,d))} else {None};
                    if closed_nonempty.load(Ordering::Acquire)
                            && exact_key.map(|key|closed.lock().unwrap().contains(&key)).unwrap_or(false) {
                        repeats.fetch_add(1,Ordering::Relaxed);
                        continue;
                    }
                    let trace=!replayed && probe_witness_exports.load(Ordering::Relaxed)<probe_witness_limit;
                    let (success,weak_seen,walk_step)=if replayed {(false,false,steps+1)} else {
                        probe_with_word(d,seed,&mut rng,steps,weak_only,if trace {Some(&mut probe_word)} else {None})
                    };
                    if success && trace && probe_witness_exports.fetch_add(1,Ordering::Relaxed)<probe_witness_limit {
                        // Cache insertion literally replays the word. Export is bounded;
                        // tracing uses the same RNG sequence and changes no rejection.
                        if let Some(cache)=positive {cache.insert(d,seed,&probe_word)?;}
                    }
                    if !replayed {
                        walk_counts[walk_step]+=1;
                        if success && walk_step>record {
                            record=walk_step;
                            writeln!(hints.lock().unwrap(),"{tid}\t{}\t{walk_step}\t{}",if weak_only {"weak"} else {"strong"},tuple_string_dynamic(d,&seed))?;
                        }
                    }
                    if success {solved.fetch_add(1,Ordering::Relaxed); weak.fetch_add(1,Ordering::Relaxed); continue;}
                    let key=exact_key.unwrap_or_else(||canonical(&seed,d));
                    if positive.as_ref().is_some_and(|cache|cache.contains(key)) {
                        solved.fetch_add(1,Ordering::Relaxed); weak.fetch_add(1,Ordering::Relaxed); continue;
                    }
                    if closed_nonempty.load(Ordering::Acquire) && closed.lock().unwrap().contains(&key) {
                        repeats.fetch_add(1,Ordering::Relaxed); continue;
                    }
                    // Deduplicate promoted quotient starts, including capped ones.
                    // They are durably replayed, not mistaken for proven outcomes.
                    if !promoted.lock().unwrap().insert(key) {promoted_repeats.fetch_add(1,Ordering::Relaxed);continue;}
                    // Persist every promoted seed before starting expensive
                    // exploration, so an interrupted candidate can be replayed.
                    let text=tuple_string_dynamic(d,&seed);
                    writeln!(pending.lock().unwrap(),"group={}_{} seed={text}",ambient.short(),d)?;
                    let result=explore_goal(d,seed,cap,deadline.min(Instant::now()+Duration::from_secs(seed_seconds)),weak_only);
                    if let (Some(cache),Some(word))=(positive,&result.weak_word) {cache.insert(d,seed,word)?;}
                    let seen_weak=weak_seen || result.weak_witness.is_some();
                    if seen_weak {weak.fetch_add(1,Ordering::Relaxed);}
                    writeln!(results.lock().unwrap(),"{}\t{}\t{}\t{}",result.status,result.states.len(),seen_weak,text)?;
                    if result.status.starts_with("certified_") {
                        let mut cache=closed.lock().unwrap();
                        let representative=*result.states.iter().min().unwrap();
                        if cache.contains(&representative) {repeats.fetch_add(1,Ordering::Relaxed); continue;}
                        let path=out_dir.join(format!("CERTIFIED_{}{}_{}_{:032x}_{:032x}.txt",ambient.short(),d,format!("{:032x}",representative[0]),representative[1],representative[2]));
                        write_certificate(&path,d,ambient,&seed,&result)?;
                        cache.extend(result.states.iter().copied());
                        closed_nonempty.store(true,Ordering::Release);
                        certified.fetch_add(1,Ordering::Relaxed);
                        if result.weak_witness.is_none() {weak_fail.fetch_add(1,Ordering::Relaxed);}
                        eprintln!("CERTIFIED {}{} status={} quotient_states={} file={}",ambient.short(),d,result.status,result.states.len(),path.display());
                    } else if result.status.ends_with("_found") {solved.fetch_add(1,Ordering::Relaxed);} else {
                        hard.fetch_add(1,Ordering::Relaxed);
                        writeln!(candidates.lock().unwrap(),"group={}_{} weak_seen={} status={} states={} seed={text}",ambient.short(),d,seen_weak,result.status,result.states.len())?;
                    }
                }
                for (total,n) in histogram.lock().unwrap().iter_mut().zip(walk_counts) { *total+=n; }
                Ok(())
            }));
        }
        let mut error=None;
        for worker in workers {
            match worker.join() {Ok(Ok(()))=>{},Ok(Err(e))=>error=Some(e.to_string()),Err(_)=>error=Some("worker panicked".into())}
        }
        done.store(true,Ordering::Relaxed);
        monitor.thread().unpark();
        if let Some(e)=error {Err(e)} else {Ok(())}
    })?;
    // A time-limited replay can end before workers claim every input. Retain
    // those leads too; otherwise a completed timed batch silently loses them.
    for seed in replay.iter().skip((claimed.load(Ordering::Relaxed) as usize).min(replay.len())) {
        writeln!(candidates.lock().unwrap(),"group={}_{} status=unprocessed seed={}",ambient.short(),d,tuple_string_dynamic(d,seed)).map_err(|e|e.to_string())?;
        hard.fetch_add(1,Ordering::Relaxed);
    }
    let mut walk_file=create("PROBE_HISTOGRAM.tsv").map_err(|e|e.to_string())?;
    writeln!(walk_file,"first_witness_walk_step\tsamples").map_err(|e|e.to_string())?;
    for (step,n) in histogram.lock().unwrap().iter().enumerate() {
        if *n>0 { writeln!(walk_file,"{}\t{n}",if step>steps {"not_found".to_string()} else {step.to_string()}).map_err(|e|e.to_string())?; }
    }
    let mut summary=create("SUMMARY.txt").map_err(|e|e.to_string())?;
    if let Some(cache)=&positive {cache.summary(&mut summary).map_err(|e|e.to_string())?;}
    writeln!(summary,"short braid witnesses pruned before generation = {}\nFull-group generation is not tested for these literal positive witnesses.",short_word_pruned.load(Ordering::Relaxed)).map_err(|e|e.to_string())?;
    writeln!(summary,"probe witness export limit = {probe_witness_limit}\nprobe witness export attempts = {}\nProbe words are selected successful walks, not shortest words or random state frequencies.",probe_witness_exports.load(Ordering::Relaxed).min(probe_witness_limit)).map_err(|e|e.to_string())?;
    if finite {
        let construction=fs::read_to_string(out_dir.join(if family=="s3-exact" {"S3_SUMMARY.txt"} else {"FIBER_SUMMARY.txt"})).map_err(|e|e.to_string())?;
        write!(summary,"{construction}").map_err(|e|e.to_string())?;
        writeln!(summary,"finite population resolved = {}",construction.lines().any(|s|s=="finite construction complete = true") && attempted.load(Ordering::Relaxed)==work_count && hard.load(Ordering::Relaxed)==0).map_err(|e|e.to_string())?;
    }
    writeln!(summary,"search mode = {}",if exact {"passport-exact"} else {"smart-hunt"}).map_err(|e|e.to_string())?;
    if let Some((m,g))=structural {
        if let Some(p)=exceptional {writeln!(summary,"exceptional cycle type = {p:?}").map_err(|e|e.to_string())?;}
        writeln!(summary,"involution transpositions = {m}\ntarget genus = {g}\nConstruction is biased sampling of sigma_1-fixed states, NOT a census.").map_err(|e|e.to_string())?;
        for (label,count) in ["raw involution pairs","pairs admitting braid conjugator","commuting third involutions","target genus draws","generating sigma1 fixed draws","also sigma2 fixed draws"].iter().zip(&structural_counts) {
            writeln!(summary,"{label} = {}",count.load(Ordering::Relaxed)).map_err(|e|e.to_string())?;
        }
    }
    if let Some((path,states,switches,max_genus))=&local {
        writeln!(summary,"source certificate = {}\nsource quotient states = {}\nconjugator transpositions = {switches}\nmaximum genus = {max_genus}\nall-six-pairs heuristic = {all_pairs}\nLocal switch is biased sampling, NOT a census.",path.display(),states.len()).map_err(|e|e.to_string())?;
        for (label,count) in ["local proposals","changed entry","genus survivors","pair-filter survivors","point-transitive survivors","full-group survivors"].iter().zip(&local_counts) {
            writeln!(summary,"{label} = {}",count.load(Ordering::Relaxed)).map_err(|e|e.to_string())?;
        }
    }
    if let Some(factory)=&symmetry {factory.save_outputs(&out_dir,d)?;factory.summary(&mut summary).map_err(|e|e.to_string())?;}
    if passport.is_some() && !exact { writeln!(summary,"passport = {}",get_arg(args,"--passport",None)?).map_err(|e|e.to_string())?; }
    if exact {
        let next=scan_start+attempted.load(Ordering::Relaxed);
        writeln!(summary,"scan raw class pairs = {}\nscan total candidate pairs = {scan_total}\nscan start = {scan_start}\nscan end = {}\nscan next = {next}\nscan range resolved = {}\nfixed-entry transitive pairs pruned = {pair_pruned}\ncentralizer-equivalent pairs removed = {conjugate_pruned}\nscan enumeration = {}\ndirect-pair witnesses pruned = {}\npassport = {}",scan_total+pair_pruned+conjugate_pruned,scan_start+work_count,next==scan_start+work_count && hard.load(Ordering::Relaxed)==0,if quotient_rows {"centralizer-v1"} else {"lex-v1"},direct_pair_pruned.load(Ordering::Relaxed),get_arg(args,"--passport",None)?).map_err(|e|e.to_string())?;
    }
    writeln!(summary,"group = {}_{}\nfamily = {family}\nweak only = {weak_only}\nthreads = {threads}\nmaster seed = {master}\nwalk steps = {steps}\ncandidate seconds = {seed_seconds}\ncap = {cap}\nreplay inputs = {}\nattempted seeds = {}\naccepted generating seeds = {}\nsolved seeds = {}\nweak success seen (excluding cached repeats) = {}\ncertified bad components (mod S_d conjugation) = {}\ncertified weak failures = {}\nrepeat bad seeds = {}\nhard unsolved candidates = {}\nelapsed seconds = {:.3}\nSamples are not distinct tuples or orbits. Capped/time-limited cases are NOT counterexamples.",ambient.short(),d,replay.len(),attempted.load(Ordering::Relaxed),generated.load(Ordering::Relaxed),solved.load(Ordering::Relaxed),weak.load(Ordering::Relaxed),certified.load(Ordering::Relaxed),weak_fail.load(Ordering::Relaxed),repeats.load(Ordering::Relaxed),hard.load(Ordering::Relaxed),start.elapsed().as_secs_f64()).map_err(|e|e.to_string())?;
    writeln!(summary,"duplicate promoted starts = {}",promoted_repeats.load(Ordering::Relaxed)).map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())?;
    eprintln!("smart hunt complete: {}",out_dir.display());
    Ok(())
}

// Find h with target_i = h source_i h^-1 for all four entries.  A generating
// tuple is transitive, so the image of point 0 determines h completely.
fn tuple_conjugator(d: usize, source: [Perm; 4], target: [Perm; 4]) -> Result<Perm, String> {
    let mut answer = None;
    for root_image in 0..d {
        let mut image = [u8::MAX; MAX_D];
        let mut queue = Vec::with_capacity(d);
        image[0] = root_image as u8;
        queue.push(0usize);
        let (mut head, mut consistent) = (0, true);
        while head < queue.len() && consistent {
            let x = queue[head];
            head += 1;
            for j in 0..4 {
                let x_next = source[j].img[x] as usize;
                let y_next = target[j].img[image[x] as usize];
                if image[x_next] == u8::MAX {
                    image[x_next] = y_next;
                    queue.push(x_next);
                } else if image[x_next] != y_next {
                    consistent = false;
                    break;
                }
            }
        }
        if !consistent || queue.len() != d { continue; }
        let mut h = Perm::identity(d);
        h.img[..d].copy_from_slice(&image[..d]);
        let mut seen = [false; MAX_D];
        if image[..d].iter().any(|&y| std::mem::replace(&mut seen[y as usize], true)) { continue; }
        if (0..4).all(|j| source[j].conjugate_by(h, d) == target[j]) {
            if answer.is_some_and(|old| old != h) { return Err("nonunique simultaneous conjugator".into()); }
            answer = Some(h);
        }
    }
    answer.ok_or_else(|| "missing simultaneous conjugator for quotient edge".into())
}

struct LiftProfile {
    holonomy_order: u128,
    loop_order: u128,
    generators: usize,
    all_even: bool,
    point_orbits: String,
    edges: Vec<(usize, usize, usize, Perm, Perm)>,
}

fn point_orbits(d: usize, generators: &[Perm]) -> Vec<Vec<usize>> {
    let mut seen = [false; MAX_D];
    let mut orbits = Vec::new();
    for root in 0..d {
        if seen[root] { continue; }
        let mut queue = vec![root];
        seen[root] = true;
        let mut head = 0;
        while head < queue.len() {
            let x = queue[head];
            head += 1;
            for g in generators {
                let y = g.img[x] as usize;
                if !seen[y] { seen[y] = true; queue.push(y); }
            }
        }
        orbits.push(queue);
    }
    orbits
}

fn point_orbit_partition(d: usize, generators: &[Perm]) -> String {
    let mut sizes:Vec<_>=point_orbits(d,generators).iter().map(Vec::len).collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    sizes.iter().map(usize::to_string).collect::<Vec<_>>().join("+")
}

fn invariant_blocks(d:usize, generators:&[Perm]) -> Option<Vec<Vec<usize>>> {
    // For a transitive action, close 0~y under the generators. Every block
    // containing 0 and y must contain this equivalence closure. Trying each y
    // therefore detects imprimitivity exactly (not a randomized test).
    if point_orbits(d,generators).len()!=1 {return None;}
    for y in 1..d {
        let mut labels:Vec<_>=(0..d).collect();labels[y]=0;
        loop {
            let mut changed=false;
            for i in 0..d {for j in i+1..d {
                if labels[i]!=labels[j] {continue;}
                for g in generators {
                    let a=labels[g.img[i] as usize];let b=labels[g.img[j] as usize];
                    if a!=b {
                        for x in &mut labels {if *x==a.max(b) {*x=a.min(b);}}
                        changed=true;
                    }
                }
            }}
            if !changed {break;}
        }
        let mut blocks=BTreeMap::new();
        for (i,label) in labels.into_iter().enumerate() {blocks.entry(label).or_insert_with(Vec::new).push(i);}
        if blocks.len()>1 {return Some(blocks.into_values().collect());}
    }
    None
}

fn braid_symmetry(args:&[String]) -> Result<(),String> {
    let source=PathBuf::from(get_arg(args,"--certificate",None)?);verify(&source)?;
    let case=parse_certificate_case(&source)?;
    if case.states.as_ref().map(Vec::len)!=Some(1) {return Err("braid-symmetry currently requires a one-state certificate".into());}
    let d=case.d;let t=case.seed;let id=Perm::identity(d);
    let mut h=[id;3];
    for i in 0..3 {h[i]=tuple_conjugator(d,t,hurwitz_move_dynamic(d,t,i))?;}
    for i in 0..2 {
        if h[i].compose(h[i+1],d).compose(h[i],d)!=h[i+1].compose(h[i],d).compose(h[i+1],d) {
            return Err("braid relation failed".into());
        }
    }
    if h[0].compose(h[2],d)!=h[2].compose(h[0],d) {return Err("far commutativity failed".into());}
    let delta=h[0].compose(h[1],d).compose(h[2],d);
    if (0..4).fold(id,|p,_|p.compose(delta,d))!=id {return Err("full twist is nontrivial".into());}
    let order=subgroup_order(d,&h);
    let all_odd=h.iter().all(|p|!p.is_even(d));
    let out=PathBuf::from(get_arg(args,"--out",None)?);fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let create=|name|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|e.to_string());
    let mut data=create("BRAID_SYMMETRY.json")?;
    writeln!(data,"{{\"degree\":{d},\"genus\":{},\"group_order\":\"{}\",\"symmetry_order\":\"{order}\",\"tuple\":{:?},\"conjugators\":{:?}}}",genus_of(&t,d),subgroup_order(d,&t),t.map(|p|p.img[..d].to_vec()),h.map(|p|p.img[..d].to_vec())).map_err(|e|e.to_string())?;
    let mut pairs=create("PAIR_OBSTRUCTIONS.tsv")?;
    writeln!(pairs,"first\tsecond\tgroup_order\tpoint_orbits_zero_based").map_err(|e|e.to_string())?;
    for i in 0..4 {for j in i+1..4 {
        writeln!(pairs,"{i}\t{j}\t{}\t{:?}",subgroup_order(d,&[t[i],t[j]]),point_orbits(d,&[t[i],t[j]])).map_err(|e|e.to_string())?;
    }}
    let mut summary=create("SUMMARY.txt")?;
    writeln!(summary,"all braid conjugators odd = {all_odd}").map_err(|e|e.to_string())?;
    if all_odd {
        writeln!(summary,"labeled graph bipartite = true\nlabeled adjacency minimum eigenvalue = -6\nlabeled normalized absolute spectral gap = 0\nProof: uniqueness of tuple relabeling identifies the labeled orbit with the conjugator group. Every generator and its inverse reverses permutation sign, giving an exact -6 eigenfunction. This does not determine the relaxation gap 1-lambda_2/6.").map_err(|e|e.to_string())?;
    }
    writeln!(summary,"experiment = one-state braid symmetry\nsource = {}\ndegree = {d}\ngenus = {}\nsymmetry group order = {order}\nsymmetry point orbits = {:?}\nsymmetry invariant blocks = {:?}\nconjugator types = {:?}\nconjugators all even = {}\nbraid relations = verified\nfar commutativity = verified\nfull twist identity = verified\nsearch complete = true\nThe certificate verifies full tuple generation and one-state six-move closure. Uniqueness of relabelings was checked. Thus every braid image is a simultaneous conjugate of the seed; conjugation preserves its first-pair orbit partition. Block detection is exact for transitive symmetry actions; None for an intransitive action does not mean primitive. Numerical spectra of action quotients do not establish the full labeled spectral gap.",source.display(),genus_of(&t,d),point_orbits(d,&h),invariant_blocks(d,&h),h.map(|p|cycle_type(p,d)),h.iter().all(|p|p.is_even(d))).map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())?;
    Ok(())
}

fn mechanism_audit(args:&[String]) -> Result<(),String> {
    let paths=get_arg(args,"--certificates",None)?;
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let mut summary=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("SUMMARY.tsv")).map_err(|e|e.to_string())?);
    let mut cusps=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("CUSPS.tsv")).map_err(|e|e.to_string())?);
    let mut states=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("STATES.tsv")).map_err(|e|e.to_string())?);
    let mut voltage_rows=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("VOLTAGES.tsv")).map_err(|e|e.to_string())?);
    let mut edge_rows=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("EDGES.tsv")).map_err(|e|e.to_string())?);
    let mut seed_pairs=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("SEED_PAIR_ORBITS.tsv")).map_err(|e|e.to_string())?);
    let mut pair_counts=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("PAIR_PARTITIONS.tsv")).map_err(|e|e.to_string())?);
    let mut incidence=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("INCIDENCE.tsv")).map_err(|e|e.to_string())?);
    let mut decorations=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("DECORATIONS.tsv")).map_err(|e|e.to_string())?);
    let mut middle=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("CUSP_MIDDLE.tsv")).map_err(|e|e.to_string())?);
    let mut symmetry_rows=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("SYMMETRIES.tsv")).map_err(|e|e.to_string())?);
    writeln!(summary,"case\tdegree\tgenus\tstates\tcusps\tall_six_pairs_intransitive_states\tfixed_sigma1\tfixed_sigma2\tfixed_sigma3\tseed_fixed_moves\tseed_conjugator_types\tseed_two_braid_group_order\tseed_two_braid_point_orbits\tposition_fiber\tposition_monodromy_order\tposition_monodromy_blocks\tposition_voltage_types").map_err(|e|e.to_string())?;
    writeln!(cusps,"case\tcusp\tsize\tfirst_pair_point_orbits").map_err(|e|e.to_string())?;
    writeln!(states,"case\tstate\texceptional_position_zero_based\tfirst_pair_point_orbits\tall_six_pair_point_orbits\tfixed_moves").map_err(|e|e.to_string())?;
    writeln!(voltage_rows,"case\tposition_zero_based\tmove_zero_based\timages_zero_based\tcycle_type").map_err(|e|e.to_string())?;
    writeln!(edge_rows,"case\tstate\tsigma1\tsigma2\tsigma3").map_err(|e|e.to_string())?;
    writeln!(seed_pairs,"case\tpair\torbit_sizes\tpoint_orbits_one_based").map_err(|e|e.to_string())?;
    writeln!(pair_counts,"case\tpair\torbit_sizes\tstates").map_err(|e|e.to_string())?;
    writeln!(incidence,"case\tstate\tedge\tleft_vertex\tright_vertex\tlength\tpoints_one_based").map_err(|e|e.to_string())?;
    writeln!(decorations,"case\tstate\tcusp\tordered_passport\tdefect\tleft_genus\tright_genus\tleft_pair_order\tright_pair_order\tleft_orbit_sizes\tright_orbit_sizes\tsigma1_period\tsigma2_period\tsigma3_period").map_err(|e|e.to_string())?;
    writeln!(middle,"case\tfrom_cusp\tto_cusp\tmove\tstate_edges").map_err(|e|e.to_string())?;
    writeln!(symmetry_rows,"case\tdegree\tmove_zero_based\timages_zero_based").map_err(|e|e.to_string())?;
    for raw in paths.split(',') {
        let path=Path::new(raw);
        verify(path)?;
        let case=parse_certificate_case(path)?;
        for a in 0..4 {for b in a+1..4 {
            let mut orbits=point_orbits(case.d,&[case.seed[a],case.seed[b]]);
            for orbit in &mut orbits {orbit.sort_unstable();}
            orbits.sort_by(|x,y|y.len().cmp(&x.len()).then_with(||x.cmp(y)));
            let sets=orbits.iter().map(|orbit|orbit.iter().map(|x|(x+1).to_string()).collect::<Vec<_>>().join(","))
                .collect::<Vec<_>>().join("|");
            writeln!(seed_pairs,"{}\t{}{}\t{}\t{}",path.display(),a+1,b+1,
                orbits.iter().map(|orbit|orbit.len().to_string()).collect::<Vec<_>>().join("+"),sets).map_err(|e|e.to_string())?;
        }}
        let graph=graph_from_states(case.d,case.states.ok_or("certificate has no states")?)?;
        let n=graph.states.len();
        let mut first=Vec::with_capacity(n);
        let mut pair_histogram=BTreeMap::<(usize,String),usize>::new();
        let mut all_intransitive=0usize;
        let mut fixed=[0usize;3];
        for (i,&key) in graph.states.iter().enumerate() {
            let tuple=unpack(key,case.d);
            let mut partitions=Vec::with_capacity(6);
            for a in 0..4 {for b in a+1..4 {
                partitions.push(point_orbit_partition(case.d,&[tuple[a],tuple[b]]));
            }}
            for (j,partition) in partitions.iter().enumerate() {
                *pair_histogram.entry((j,partition.clone())).or_default()+=1;
            }
            if partitions.iter().all(|p|p.contains('+')) {all_intransitive+=1;}
            let types=tuple.map(|p|cycle_type(p,case.d));
            let position=(0..4).find(|&j|types.iter().filter(|t|**t==types[j]).count()==1)
                .map(|j|j.to_string()).unwrap_or_default();
            let fixed_moves=(0..3).filter(|&mv|graph.transitions[i][mv] as usize==i).collect::<Vec<_>>();
            for &mv in &fixed_moves {fixed[mv]+=1;}
            first.push(partitions[0].clone());
            writeln!(states,"{}\t{i}\t{position}\t{}\t{}\t{}",path.display(),partitions[0],partitions.join(";"),
                fixed_moves.iter().map(usize::to_string).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
            writeln!(edge_rows,"{}\t{i}\t{}\t{}\t{}",path.display(),graph.transitions[i][0],graph.transitions[i][1],graph.transitions[i][2]).map_err(|e|e.to_string())?;
        }
        for ((j,partition),count) in pair_histogram {
            writeln!(pair_counts,"{}\t{}\t{}\t{}",path.display(),["12","13","14","23","24","34"][j],partition,count).map_err(|e|e.to_string())?;
        }
        let mut seen=vec![false;n];
        let mut cusp_of=vec![usize::MAX;n];
        let mut cusp_count=0;
        for root in 0..n {
            if seen[root] {continue;}
            let mut queue=vec![root];seen[root]=true;
            let mut head=0;
            while head<queue.len() {
                let i=queue[head];head+=1;
                cusp_of[i]=cusp_count;
                if first[i]!=first[root] {return Err(format!("{}: first-pair partition changed inside a cusp",path.display()));}
                if defect_of(&unpack(graph.states[i],case.d),case.d)!=defect_of(&unpack(graph.states[root],case.d),case.d) {
                    return Err("defect or pair genus changed inside an end-braid cusp".into());
                }
                for mv in [0,2,3,5] {
                    let j=graph.transitions[i][mv] as usize;
                    if !seen[j] {seen[j]=true;queue.push(j);}
                }
            }
            writeln!(cusps,"{}\t{cusp_count}\t{}\t{}",path.display(),queue.len(),first[root]).map_err(|e|e.to_string())?;
            cusp_count+=1;
        }
        let mut transition_counts=BTreeMap::new();
        for (i,&key) in graph.states.iter().enumerate() {
            let t=unpack(key,case.d);let (delta,gl,gr)=defect_of(&t,case.d);
            let left=point_orbits(case.d,&t[..2]);let right=point_orbits(case.d,&t[2..]);
            let mut left_of=[0usize;MAX_D];let mut right_of=[0usize;MAX_D];
            for (j,points) in left.iter().enumerate() {for &x in points {left_of[x]=j;}}
            for (j,points) in right.iter().enumerate() {for &x in points {right_of[x]=j;}}
            let boundaries=permutation_cycles(case.d,t[0].compose(t[1],case.d));
            if boundaries.len() as i64-left.len() as i64-right.len() as i64+1!=delta || delta+gl+gr!=genus_of(&t,case.d) {
                return Err("incidence Euler identity failed".into());
            }
            for (j,points) in boundaries.iter().enumerate() {
                if points.iter().any(|&x|left_of[x]!=left_of[points[0]] || right_of[x]!=right_of[points[0]]) {return Err("boundary crosses pair orbits".into());}
                writeln!(incidence,"{}\t{i}\t{j}\t{}\t{}\t{}\t{}",path.display(),left_of[points[0]],right_of[points[0]],points.len(),points.iter().map(|x|(x+1).to_string()).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
            }
            let periods=std::array::from_fn::<_,3,_>(|mv| {let mut j=graph.transitions[i][mv] as usize;let mut length=1;while j!=i {j=graph.transitions[j][mv] as usize;length+=1;}length});
            writeln!(decorations,"{}\t{i}\t{}\t{}\t{delta}\t{gl}\t{gr}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",path.display(),cusp_of[i],passport_string(&t,case.d),subgroup_order(case.d,&t[..2]),subgroup_order(case.d,&t[2..]),point_orbit_partition(case.d,&t[..2]),point_orbit_partition(case.d,&t[2..]),periods[0],periods[1],periods[2]).map_err(|e|e.to_string())?;
            for mv in [1,4] {*transition_counts.entry((cusp_of[i],cusp_of[graph.transitions[i][mv] as usize],mv)).or_insert(0usize)+=1;}
        }
        for ((from,to,mv),count) in transition_counts {writeln!(middle,"{}\t{from}\t{to}\t{mv}\t{count}",path.display()).map_err(|e|e.to_string())?;}
        let seed_index=graph.states.iter().position(|&key|key==canonical(&case.seed,case.d)).ok_or("missing seed")?;
        let seed_fixed=(0..3).filter(|&mv|graph.transitions[seed_index][mv] as usize==seed_index).collect::<Vec<_>>();
        let conjugators=seed_fixed.iter().map(|&mv|tuple_conjugator(case.d,case.seed,hurwitz_move_dynamic(case.d,case.seed,mv)))
            .collect::<Result<Vec<_>,_>>()?;
        for (&mv,&h) in seed_fixed.iter().zip(&conjugators) {
            writeln!(symmetry_rows,"{}\t{}\t{mv}\t{}",path.display(),case.d,h.img[..case.d].iter().map(u8::to_string).collect::<Vec<_>>().join(",")).map_err(|e|e.to_string())?;
        }
        let two_braid=if seed_fixed.contains(&0) && seed_fixed.contains(&1) {
            let h=conjugators[seed_fixed.iter().position(|&mv|mv==0).unwrap()];
            let k=conjugators[seed_fixed.iter().position(|&mv|mv==1).unwrap()];
            (subgroup_order(case.d,&[h,k]).to_string(),point_orbit_partition(case.d,&[h,k]))
        } else {(String::new(),String::new())};
        let position=if let Ok(fibers)=position_fibers(case.d,&graph) {
            let m=fibers[0].len();
            let mut lifted=vec![vec![usize::MAX;m];4];
            let mut labels=vec![usize::MAX;n];
            for (r,&s) in fibers[0].iter().enumerate() {
                lifted[0][r]=s;labels[s]=r;
                for p in 0..3 {
                    let next=graph.transitions[lifted[p][r]][p] as usize;
                    if labels[next]!=usize::MAX {return Err("position transport is not bijective".into());}
                    lifted[p+1][r]=next;labels[next]=r;
                }
            }
            if labels.contains(&usize::MAX) {return Err("position transport missed a state".into());}
            let mut voltages=Vec::new();
            for p in 0..4 {for mv in 0..3 {
                let mut voltage=Perm::identity(m);
                for r in 0..m {voltage.img[r]=labels[graph.transitions[lifted[p][r]][mv] as usize] as u8;}
                writeln!(voltage_rows,"{}\t{p}\t{mv}\t{}\t{:?}",path.display(),
                    voltage.img[..m].iter().map(u8::to_string).collect::<Vec<_>>().join(","),cycle_type(voltage,m)).map_err(|e|e.to_string())?;
                voltages.push(voltage);
            }}
            let types=voltages.iter().filter(|v|!v.is_identity(m)).map(|&v|format!("{:?}",cycle_type(v,m)))
                .collect::<Vec<_>>().join(";");
            (m.to_string(),subgroup_order(m,&voltages).to_string(),
                format!("{:?}",invariant_blocks(m,&voltages).as_ref().map(|b|b.iter().map(Vec::len).collect::<Vec<_>>())),types)
        } else {(String::new(),String::new(),String::new(),String::new())};
        writeln!(summary,"{}\t{}\t{}\t{n}\t{cusp_count}\t{all_intransitive}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            path.display(),case.d,genus_of(&case.seed,case.d),fixed[0],fixed[1],fixed[2],
            seed_fixed.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
            conjugators.iter().map(|&h|format!("{:?}",cycle_type(h,case.d))).collect::<Vec<_>>().join(";"),
            two_braid.0,two_braid.1,position.0,position.1,position.2,position.3).map_err(|e|e.to_string())?;
    }
    Ok(())
}

fn s18_fiber_census(args:&[String]) -> Result<(),String> {
    let certificate=PathBuf::from(get_arg(args,"--certificate",None)?);
    verify(&certificate)?;
    let case=parse_certificate_case(&certificate)?;
    let graph=graph_from_states(case.d,case.states.ok_or("certificate has no states")?)?;
    let fibers=position_fibers(case.d,&graph)?;
    if case.d!=18 || fibers[0].len()!=10 {return Err("requires the ten-state S18 position fiber".into());}
    let mut fiber_index=vec![usize::MAX;graph.states.len()];
    for (r,&s) in fibers[0].iter().enumerate() {fiber_index[s]=r;}
    let mut generators=[Perm::identity(10);2];
    for (j,mv) in [1,2].into_iter().enumerate() {
        for (r,&s) in fibers[0].iter().enumerate() {
            let target=graph.transitions[s][mv] as usize;
            if fiber_index[target]==usize::MAX {return Err("fiber generator leaves its position".into());}
            generators[j].img[r]=fiber_index[target] as u8;
        }
    }
    let [u,v]=generators;
    if cycle_type(u,10)!=vec![1,2,7] || cycle_type(v,10)!=vec![1,2,7]
        || u.compose(v,10).compose(u,10)!=v.compose(u,10).compose(v,10) {
        return Err("unexpected S18 fiber braid pair".into());
    }
    let mut powers=Vec::with_capacity(14);
    let mut h=Perm::identity(10);
    for _ in 0..14 {powers.push(h);h=h.compose(u,10);}
    if !h.is_identity(10) || powers.iter().copied().collect::<HashSet<_>>().len()!=14 {
        return Err("unexpected fiber centralizer".into());
    }
    let mut class_size=0usize;
    let mut solutions=0usize;
    let mut classes=BTreeMap::<Vec<u8>,(usize,Perm)>::new();
    for fixed in 0..10 {
        for a in 0..10 {for b in a+1..10 {
            if a==fixed || b==fixed {continue;}
            let cycle=(0..10).filter(|&x|x!=fixed && x!=a && x!=b).map(|x|x as u8).collect::<Vec<_>>();
            let mut tail=cycle[1..].to_vec();
            loop {
                let mut w=Perm::identity(10);
                w.img[a]=b as u8;w.img[b]=a as u8;
                let mut previous=cycle[0];
                for &next in &tail {w.img[previous as usize]=next;previous=next;}
                w.img[previous as usize]=cycle[0];
                class_size+=1;
                if u.compose(w,10).compose(u,10)==w.compose(u,10).compose(w,10) {
                    solutions+=1;
                    let representative=powers.iter().map(|&g|g.compose(w,10).compose(g.inverse(10),10).img[..10].to_vec())
                        .min().ok_or("empty centralizer")?;
                    classes.entry(representative).and_modify(|row|row.0+=1).or_insert((1,w));
                }
                if !next_permutation(&mut tail) {break;}
            }
        }}
    }
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let mut rows=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("CLASSES.tsv")).map_err(|e|e.to_string())?);
    writeln!(rows,"centralizer_orbit_size\tgenerated_group_order\trepresentative_images_zero_based").map_err(|e|e.to_string())?;
    let mut full=0;
    for (_key,(size,w)) in &classes {
        let order=subgroup_order(10,&[u,*w]);
        if order==factorial_u128(10) {full+=1;}
        writeln!(rows,"{size}\t{order}\t{}",w.img[..10].iter().map(u8::to_string).collect::<Vec<_>>().join(","))
            .map_err(|e|e.to_string())?;
    }
    if class_size!=259_200 || classes.values().map(|row|row.0).sum::<usize>()!=solutions {
        return Err("fiber braid census count mismatch".into());
    }
    let x=u.compose(v,10).compose(u,10);
    let y=u.compose(v,10);
    if !y.inverse(10).compose(x,10).compose(u.inverse(10),10).is_identity(10) {
        return Err("auxiliary triple is not product-one".into());
    }
    let mut summary=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("SUMMARY.txt")).map_err(|e|e.to_string())?);
    writeln!(summary,"certificate = {}\nfiber size = 10\nu images = {:?}\nv images = {:?}\nu and v cycle type = {:?}\n<u,v> order = {}\nuv cycle type = {:?}\nuvu cycle type = {:?}\n(uv)^3 identity = {}\nauxiliary product-one triple = ((uv)^-1, uvu, u^-1)\nauxiliary three-point genus = {}\nfixed-u candidate class size = {class_size}\nbraid-compatible v = {solutions}\ncentralizer orbits = {}\nfull S10 centralizer orbits = {full}\nsearch complete = true\nAll v of type 7+2+1 are enumerated with one fixed point, one transposition, and one cyclic ordering of the remaining seven points. Uniqueness is only for this cycle type and this ten-point braid representation, not for S18 tuples.",
        certificate.display(),u.img[..10].to_vec(),v.img[..10].to_vec(),cycle_type(u,10),subgroup_order(10,&[u,v]),
        cycle_type(y,10),cycle_type(x,10),y.compose(y,10).compose(y,10).is_identity(10),
        1-10+(index(x,10)+index(y,10)+index(u,10)) as i64/2,classes.len()).map_err(|e|e.to_string())?;
    s18_matching_slice(&case.seed,&out)?;
    Ok(())
}

fn s18_matching_slice(seed:&[Perm;4],out:&Path) -> Result<(),String> {
    let d=18;
    let h=tuple_conjugator(d,*seed,hurwitz_move_dynamic(d,*seed,0))?;
    let cycles=permutation_cycles(d,h);
    if cycles.len()!=6 || cycles.iter().any(|cycle|cycle.len()!=3)
        || seed[..3].iter().any(|&p|cycle_type(p,d)!=vec![2;9]) {
        return Err("matching slice requires the S18 six-triple involution construction".into());
    }
    let nine_orbits=point_orbits(d,&[seed[0].compose(seed[1],d),h]);
    if nine_orbits.len()!=2 || nine_orbits.iter().any(|orbit|orbit.len()!=9)
        || ![seed[0],seed[1]].iter().all(|p|
            nine_orbits[0].iter().all(|&x|nine_orbits[1].contains(&(p.img[x] as usize)))) {
        return Err("expected two nine-point blocks exchanged by the first pair".into());
    }
    let cycle_side:Vec<_>=cycles.iter().map(|cycle|usize::from(nine_orbits[1].contains(&cycle[0]))).collect();
    let ab=seed[0].compose(seed[1],d);
    if ab.compose(h,d)!=h.compose(ab,d) {return Err("ab does not centralize h".into());}
    let mut base_ab=Perm::identity(6);
    for i in 0..6 {
        base_ab.img[i]=cycles.iter().position(|cycle|cycle.contains(&(ab.img[cycles[i][0]] as usize)))
            .ok_or("ab did not map an h-cycle")? as u8;
    }
    if cycle_type(base_ab,6)!=vec![3,3] {return Err("unexpected ab action on six h-cycles".into());}
    let mut rows=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("MATCHING_SLICE.tsv")).map_err(|e|e.to_string())?);
    writeln!(rows,"matching_zero_based_with_phase\tcross_nine_block_edges\tbase_fourth_cycle_type\tfourth_cycle_type\tgenus\ttransitive\tgenerated_group_order\tproper_group_block_sizes\tac_point_orbits\tbc_point_orbits\tcd_point_orbits")
        .map_err(|e|e.to_string())?;
    let mut counts=BTreeMap::<(Vec<usize>,bool,bool),usize>::new();
    let mut cross_counts=[[0usize;4];4];
    let mut visited=0usize;
    let complete=enumerate_commuting_involutions(&cycles,9,Perm::identity(d),&mut |c| {
        visited+=1;
        let t=[seed[0],seed[1],c,product_one_fourth(d,seed[0],seed[1],c)];
        let transitive=point_orbits(d,&t[..3]).len()==1;
        let group_order=subgroup_order(d,&t[..3]);
        let full=group_order==factorial_u128(d);
        let blocks=if transitive && !full {
            invariant_blocks(d,&t[..3]).map(|parts|parts.iter().map(Vec::len).collect::<Vec<_>>())
        } else {None};
        let fourth=cycle_type(t[3],d);
        *counts.entry((fourth.clone(),transitive,full)).or_default()+=1;
        let mut matching=Vec::new();
        let mut cross=0usize;
        let mut base_c=Perm::identity(6);
        for i in 0..6 {
            let image=c.img[cycles[i][0]] as usize;
            let j=cycles.iter().position(|cycle|cycle.contains(&image)).ok_or("missing matching orbit")?;
            base_c.img[i]=j as u8;
            if i<j {
                cross+=usize::from(cycle_side[i]!=cycle_side[j]);
                let phase=cycles[j].iter().position(|&x|x==image).ok_or("missing matching phase")?;
                matching.push(format!("{i}-{j}:{phase}"));
            }
        }
        if matching.len()!=3 || ![1,3].contains(&cross) {return Err("unexpected six-cycle matching".into());}
        cross_counts[cross][0]+=1;
        cross_counts[cross][1]+=usize::from(transitive);
        cross_counts[cross][2]+=usize::from(full);
        cross_counts[cross][3]+=usize::from(full && genus_of(&t,d)==1);
        if cross==3 && full {return Err("full S18 generation preserved a nine-point block system".into());}
        let base_fourth=cycle_type(base_ab.compose(base_c,6),6);
        writeln!(rows,"{}\t{cross}\t{base_fourth:?}\t{:?}\t{}\t{transitive}\t{group_order}\t{blocks:?}\t{}\t{}\t{}",
            matching.join(","),fourth,genus_of(&t,d),
            point_orbit_partition(d,&[t[0],c]),point_orbit_partition(d,&[t[1],c]),
            point_orbit_partition(d,&[c,t[3]])).map_err(|e|e.to_string())?;
        Ok(true)
    })?;
    if !complete || visited!=405 {return Err("S18 matching slice is incomplete".into());}
    let mut skeleton=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("NINE_BLOCK_SUMMARY.tsv")).map_err(|e|e.to_string())?);
    writeln!(skeleton,"cross_edges\tmatching_skeletons\tphase_assignments\ttransitive\tfull_S18\tfull_S18_genus_one").map_err(|e|e.to_string())?;
    for cross in [1,3] {
        let row=cross_counts[cross];
        writeln!(skeleton,"{cross}\t{}\t{}\t{}\t{}\t{}",row[0]/27,row[0],row[1],row[2],row[3]).map_err(|e|e.to_string())?;
    }
    let mut summary=BufWriter::new(OpenOptions::new().write(true).create_new(true)
        .open(out.join("MATCHING_SUMMARY.tsv")).map_err(|e|e.to_string())?);
    writeln!(summary,"fourth_cycle_type\tgenus\ttransitive\tfull_S18\tcount").map_err(|e|e.to_string())?;
    for ((fourth,transitive,full),count) in counts {
        let genus=1-18+(27+18-fourth.len()) as i64/2;
        writeln!(summary,"{fourth:?}\t{genus}\t{transitive}\t{full}\t{count}").map_err(|e|e.to_string())?;
    }
    Ok(())
}

fn lift_profile(d: usize, graph: &SpectrumGraph) -> Result<LiftProfile, String> {
    let n = graph.transitions.len();
    let id = Perm::identity(d);
    let mut voltage = vec![[id; 6]; n];
    for i in 0..n {
        let source = unpack(graph.states[i], d);
        for mv in 0..6 {
            let j = graph.transitions[i][mv] as usize;
            let moved = hurwitz_move_dynamic(d, source, mv);
            voltage[i][mv] = tuple_conjugator(d, unpack(graph.states[j], d), moved)?;
        }
    }
    for i in 0..n {
        for mv in 0..6 {
            let j = graph.transitions[i][mv] as usize;
            if voltage[i][mv].compose(voltage[j][(mv + 3) % 6], d) != id {
                return Err(format!("voltage inverse check failed at vertex {i}, move {mv}"));
            }
        }
    }
    let mut tree = vec![None; n];
    tree[0] = Some(id);
    let mut queue = vec![0usize];
    let mut head = 0;
    while head < queue.len() {
        let i = queue[head];
        head += 1;
        for mv in 0..6 {
            let j = graph.transitions[i][mv] as usize;
            if tree[j].is_none() {
                tree[j] = Some(tree[i].unwrap().compose(voltage[i][mv], d));
                queue.push(j);
            }
        }
    }
    if queue.len() != n { return Err("quotient graph is disconnected".into()); }
    let mut generators = HashSet::new();
    let mut loops = HashSet::new();
    let mut edges = Vec::with_capacity(6 * n);
    for i in 0..n {
        for mv in 0..6 {
            let j = graph.transitions[i][mv] as usize;
            let h = tree[i].unwrap().compose(voltage[i][mv], d).compose(tree[j].unwrap().inverse(d), d);
            edges.push((i, mv, j, voltage[i][mv], h));
            if !h.is_identity(d) { generators.insert(h); }
            if i == j && !h.is_identity(d) { loops.insert(h); }
        }
    }
    let generators: Vec<_> = generators.into_iter().collect();
    let loops: Vec<_> = loops.into_iter().collect();
    Ok(LiftProfile {
        holonomy_order: subgroup_order(d, &generators),
        loop_order: subgroup_order(d, &loops),
        generators: generators.len(),
        all_even: generators.iter().all(|g| g.is_even(d)),
        point_orbits: point_orbit_partition(d, &generators),
        edges,
    })
}

// The regular labeled lift contains every representation of H.  Restricting
// its adjacency to the standard (sum-zero) part of the natural point action
// gives a small symmetric block and a certified upper bound on the full gap.
fn standard_point_lambda_max(d: usize, n: usize, edges: &[(usize, usize, usize, Perm, Perm)]) -> Result<f64, String> {
    let p = d - 1;
    let size = n * p;
    let mut basis = vec![vec![0.0; p]; d];
    for k in 0..p {
        let scale = ((k + 1) * (k + 2)) as f64;
        let scale = scale.sqrt();
        for x in 0..=k { basis[x][k] = 1.0 / scale; }
        basis[k + 1][k] = -((k + 1) as f64) / scale;
    }
    let mut matrix = vec![0.0; size * size];
    for &(i, _, j, voltage, _) in edges {
        let inverse = voltage.inverse(d);
        for a in 0..d {
            let b = inverse.img[a] as usize;
            for r in 0..p {
                let value = basis[a][r];
                if value == 0.0 { continue; }
                for s in 0..p { matrix[(i * p + r) * size + j * p + s] += value * basis[b][s]; }
            }
        }
    }
    for i in 0..size {
        for j in 0..i {
            if (matrix[i * size + j] - matrix[j * size + i]).abs() > 1e-9 {
                return Err(format!("natural representation block is not symmetric at {i},{j}"));
            }
        }
    }
    Ok(jacobi_eigenvalues(matrix, size)[0])
}

fn binomial(n: usize, k: usize) -> usize {
    let k = k.min(n - k);
    (0..k).fold(1usize, |value, i| value * (n - i) / (i + 1))
}

fn subset_action_graph(d: usize, quotient: &SpectrumGraph, edges: &[(usize, usize, usize, Perm, Perm)], k: usize) -> SpectrumGraph {
    let mut positions: Vec<usize> = (0..k).collect();
    let mut subsets = Vec::with_capacity(binomial(d, k));
    loop {
        subsets.push(positions.iter().fold(0u64, |mask, &x| mask | (1u64 << x)));
        let Some(i) = (0..k).rev().find(|&i| positions[i] < d - k + i) else { break; };
        positions[i] += 1;
        for j in (i + 1)..k { positions[j] = positions[j - 1] + 1; }
    }
    let lookup: HashMap<u64, u32> = subsets.iter().enumerate().map(|(i, &mask)| (mask, i as u32)).collect();
    let width = subsets.len();
    let mut transitions = vec![[0u32; 6]; quotient.transitions.len() * width];
    for &(i, mv, j, voltage, _) in edges {
        let inverse = voltage.inverse(d);
        for (subset_index, &mask) in subsets.iter().enumerate() {
            let mut from = mask;
            let mut to = 0u64;
            while from != 0 {
                let x = from.trailing_zeros() as usize;
                from &= from - 1;
                to |= 1u64 << inverse.img[x];
            }
            transitions[i * width + subset_index][mv] = j as u32 * width as u32 + lookup[&to];
        }
    }
    SpectrumGraph { states: Vec::new(), transitions }
}

fn ordered_action_graph(d: usize, quotient: &SpectrumGraph, edges: &[(usize, usize, usize, Perm, Perm)], k: usize) -> SpectrumGraph {
    let bits=if d<=32 {5} else {6};
    let mask=(1u32<<bits)-1;
    let mut tuples = vec![0u32];
    for slot in 0..k {
        let mut next = Vec::with_capacity(tuples.len() * (d - slot));
        for code in tuples {
            for x in 0..d {
                if (0..slot).any(|i| ((code >> (bits * i)) & mask) == x as u32) { continue; }
                next.push(code | ((x as u32) << (bits * slot)));
            }
        }
        tuples = next;
    }
    let lookup: HashMap<u32, u32> = tuples.iter().enumerate().map(|(i, &code)| (code, i as u32)).collect();
    let width = tuples.len();
    let mut transitions = vec![[0u32; 6]; quotient.transitions.len() * width];
    for &(i, mv, j, voltage, _) in edges {
        let inverse = voltage.inverse(d);
        for (tuple_index, &code) in tuples.iter().enumerate() {
            let mut to = 0u32;
            for slot in 0..k {
                let x = ((code >> (bits * slot)) & mask) as usize;
                to |= (inverse.img[x] as u32) << (bits * slot);
            }
            transitions[i * width + tuple_index][mv] = j as u32 * width as u32 + lookup[&to];
        }
    }
    SpectrumGraph { states: Vec::new(), transitions }
}

// An integer, mean-zero test vector gives a rational Rayleigh quotient.
// Its value is a rigorous lower bound on lambda_2 of this action graph and
// hence of the labeled lift; the floating-point iteration only finds a vector.
fn exact_rayleigh(graph: &SpectrumGraph, vector: &[i64]) -> Result<(i128, i128), String> {
    if vector.len() != graph.transitions.len() || vector.len() < 2 { return Err("Rayleigh vector has wrong length".into()); }
    if vector.iter().map(|&x| x as i128).sum::<i128>() != 0 { return Err("Rayleigh vector is not mean-zero".into()); }
    let mut numerator = 0i128;
    let mut denominator = 0i128;
    for (i, row) in graph.transitions.iter().enumerate() {
        let x = vector[i] as i128;
        denominator = denominator.checked_add(x.checked_mul(x).ok_or("Rayleigh square overflow")?).ok_or("Rayleigh denominator overflow")?;
        for (mv, &j) in row.iter().enumerate() {
            let j = j as usize;
            if j >= vector.len() || graph.transitions[j][(mv + 3) % 6] as usize != i {
                return Err(format!("action graph inverse edge failed at {i}, move {mv}"));
            }
            numerator = numerator.checked_add(x.checked_mul(vector[j] as i128).ok_or("Rayleigh product overflow")?).ok_or("Rayleigh numerator overflow")?;
        }
    }
    if denominator == 0 { return Err("Rayleigh vector is zero".into()); }
    Ok((numerator, denominator))
}

fn rayleigh_vector(graph: &SpectrumGraph, iterations: usize, seed: u64) -> Result<Vec<i64>, String> {
    let n = graph.transitions.len();
    if n < 2 || iterations == 0 { return Err("Rayleigh needs at least two vertices and positive iterations".into()); }
    let mut rng = SplitMix64::new(seed);
    let mut v = (0..n).map(|_| (rng.next_u64() as f64 / u64::MAX as f64) - 0.5).collect::<Vec<_>>();
    project_constant(&mut v);
    let mut length = norm(&v);
    for x in &mut v { *x /= length; }
    for step in 0..iterations {
        let mut next = graph_matvec(graph, &v);
        for (y, &x) in next.iter_mut().zip(&v) { *y += 6.0 * x; }
        project_constant(&mut next);
        length = norm(&next);
        if !length.is_finite() || length == 0.0 { return Err("Rayleigh iteration collapsed".into()); }
        for x in &mut next { *x /= length; }
        v = next;
        if (step + 1) % 1000 == 0 { eprintln!("Rayleigh iteration {}/{}", step + 1, iterations); }
    }
    let mut integers = v.iter().map(|&x| (x * 10_000_000_000.0).round() as i64).collect::<Vec<_>>();
    let sum = integers.iter().map(|&x| x as i128).sum::<i128>();
    integers[0] = (integers[0] as i128 - sum).try_into().map_err(|_| "Rayleigh balancing overflow")?;
    Ok(integers)
}

fn rayleigh(args: &[String]) -> Result<(), String> {
    let certificate = get_arg(args, "--certificate", None)?;
    let out = PathBuf::from(get_arg(args, "--out", None)?);
    let k: usize = parse_num(get_arg(args, "--ordered-k", Some("3"))?, "--ordered-k")?;
    let iterations: usize = parse_num(get_arg(args, "--iterations", Some("6000"))?, "--iterations")?;
    let seed: u64 = parse_num(get_arg(args, "--seed", Some("6751043"))?, "--seed")?;
    if !(1..=4).contains(&k) { return Err("ordered-k must be 1..4".into()); }
    let case = parse_certificate_case(Path::new(&certificate))?;
    let quotient = enumerate_spectrum_graph(case.d, case.seed)?;
    let profile = lift_profile(case.d, &quotient)?;
    let graph = ordered_action_graph(case.d, &quotient, &profile.edges, k);
    if graph_components(&graph) != 1 { return Err("ordered action graph is disconnected".into()); }
    let started = Instant::now();
    let vector = if let Some(path) = optional_arg(args, "--vector") {
        fs::read_to_string(path).map_err(|e| e.to_string())?.lines().map(|s| s.parse::<i64>().map_err(|e| e.to_string())).collect::<Result<Vec<_>,_>>()?
    } else { rayleigh_vector(&graph, iterations, seed)? };
    let (num, den) = exact_rayleigh(&graph, &vector)?;
    let gap_num = 6i128.checked_mul(den).and_then(|x| x.checked_sub(num)).ok_or("gap numerator overflow")?;
    let gap_den = 6i128.checked_mul(den).ok_or("gap denominator overflow")?;
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let mut summary = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join("RAYLEIGH_CERT.txt")).map_err(|e| e.to_string())?);
    writeln!(summary, "case={}\ncertificate={}\nordered_k={}\nvertices={}\nquotient_vertices={}\nholonomy_order={}\niterations={}\nseed={}\nvector_sum=0\nrayleigh_numerator={}\nrayleigh_denominator={}\nrayleigh_decimal={:.12}\nfull_relaxation_gap_at_most={}/{}\nfull_relaxation_gap_decimal_at_most={:.12}\nelapsed_seconds={:.3}",
        case.id, certificate, k, graph.transitions.len(), quotient.transitions.len(), profile.holonomy_order, iterations, seed,
        num, den, num as f64 / den as f64, gap_num, gap_den, gap_num as f64 / gap_den as f64, started.elapsed().as_secs_f64()).map_err(|e| e.to_string())?;
    let mut file = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join("RAYLEIGH_VECTOR.txt")).map_err(|e| e.to_string())?);
    for x in vector { writeln!(file, "{x}").map_err(|e| e.to_string())?; }
    eprintln!("Rayleigh certificate complete: {} ({:.9} <= lambda2)", out.display(), num as f64 / den as f64);
    Ok(())
}

fn position_fibers(d:usize,graph:&SpectrumGraph) -> Result<[Vec<usize>;4],String> {
    let seed=unpack(*graph.states.first().ok_or("empty graph")?,d);
    let types=seed.map(|p|cycle_type(p,d));
    let exceptional=types.iter().find(|t|types.iter().filter(|u|u==t).count()==1).ok_or("no exceptional entry type")?;
    if types.iter().filter(|t|*t!=exceptional).any(|t|types.iter().filter(|u|u==&t).count()!=3) {
        return Err("requires a three-plus-one passport".into());
    }
    let mut fibers:[Vec<usize>;4]=Default::default();
    let mut positions=Vec::new();
    for (i,key) in graph.states.iter().enumerate() {
        let t=unpack(*key,d);
        let matches:Vec<_>=(0..4).filter(|&j|cycle_type(t[j],d)==*exceptional).collect();
        if matches.len()!=1 {return Err("exceptional type not unique at a state".into());}
        positions.push(matches[0]);fibers[matches[0]].push(i);
    }
    if fibers.iter().any(|f|f.len()!=fibers[0].len()) {return Err("unequal position fibers".into());}
    for (i,row) in graph.transitions.iter().enumerate() {
        for (mv,&j) in row.iter().enumerate() {
            let p=positions[i];let generator=mv%3;
            let expected=if p==generator {p+1} else if p==generator+1 {p-1} else {p};
            if positions[j as usize]!=expected {return Err("position quotient edge failed".into());}
        }
    }
    Ok(fibers)
}

fn project_fibers(v:&mut [f64],fibers:&[Vec<usize>;4]) {
    for f in fibers {
        let mean=f.iter().map(|&i|v[i]).sum::<f64>()/f.len() as f64;
        for &i in f {v[i]-=mean;}
    }
}

fn position_spectrum(args:&[String]) -> Result<(),String> {
    let certificate=PathBuf::from(get_arg(args,"--certificate",None)?);
    verify(&certificate)?;
    let case=parse_certificate_case(&certificate)?;
    let graph=graph_from_states(case.d,case.states.unwrap())?;
    let n=graph.states.len();
    // ponytail: dense exact-size decomposition capped at 256; use a projected
    // sparse eigensolver if larger certified quotients become relevant.
    if !(8..=256).contains(&n) {return Err("position-spectrum needs 8..256 vertices".into());}
    let fibers=position_fibers(case.d,&graph)?;
    let mut basis=Vec::new();
    for f in &fibers {
        for j in 1..f.len() {
            let mut v=vec![0.0;n];let scale=((j*(j+1)) as f64).sqrt();
            for &i in &f[..j] {v[i]=1.0/scale;}
            v[f[j]]=-(j as f64)/scale;basis.push(v);
        }
    }
    let images:Vec<_>=basis.iter().map(|v|graph_matvec(&graph,v)).collect();
    let size=basis.len();
    let block:Vec<_>=basis.iter().flat_map(|u|images.iter().map(move |v|dot(u,v))).collect();
    let internal=jacobi_eigenvalues(block,size);
    let factor=vec![6.0,2.0+2.0*2f64.sqrt(),2.0,2.0-2.0*2f64.sqrt()];
    let mut joined=factor.clone();joined.extend(&internal);joined.sort_by(|a,b|b.total_cmp(a));
    let mut matrix=vec![0.0;n*n];
    for (i,row) in graph.transitions.iter().enumerate() {for &j in row {matrix[i*n+j as usize]+=1.0;}}
    let full=jacobi_eigenvalues(matrix,n);
    let decomposition_error=full.iter().zip(&joined).map(|(a,b)|(a-b).abs()).fold(0.0,f64::max);
    if decomposition_error>1e-8 {return Err(format!("spectral decomposition mismatch {decomposition_error}"));}
    let mut rng=SplitMix64::new(20260927);
    let mut vector:Vec<_>=(0..n).map(|_|rng.next_u64() as f64/u64::MAX as f64-0.5).collect();
    project_fibers(&mut vector,&fibers);
    for _ in 0..4000 {
        let av=graph_matvec(&graph,&vector);
        for (x,y) in vector.iter_mut().zip(av) {*x=6.0**x-y;}
        project_fibers(&mut vector,&fibers);
        let length=norm(&vector);
        if !length.is_finite() || length==0.0 {return Err("negative mode iteration collapsed".into());}
        for x in &mut vector {*x/=length;}
    }
    let mut integers:Vec<_>=vector.iter().map(|x|(x*1e9).round() as i64).collect();
    for f in &fibers {
        let sum:i64=f.iter().map(|&i|integers[i]).sum();integers[f[0]]-=sum;
    }
    let (num,den)=exact_rayleigh(&graph,&integers)?;
    if num>=0 {return Err("no negative Rayleigh witness found".into());}
    let rayleigh=num as f64/den as f64;
    let av=graph_matvec(&graph,&vector);
    let residual=av.iter().zip(&vector).map(|(a,v)|(a-rayleigh*v).powi(2)).sum::<f64>().sqrt();
    let out=PathBuf::from(get_arg(args,"--out",None)?);
    fs::create_dir_all(&out).map_err(|e|e.to_string())?;
    let create=|name|OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e|e.to_string());
    let mut file=create("SPECTRUM.tsv")?;
    writeln!(file,"block\teigenvalue").map_err(|e|e.to_string())?;
    for (label,values) in [("position",&factor),("within_position",&internal),("full",&full)] {
        for value in values {writeln!(file,"{label}\t{value:.14}").map_err(|e|e.to_string())?;}
    }
    let mut witness=create("NEGATIVE_VECTOR.tsv")?;
    writeln!(witness,"state\tposition\tinteger\tcolor\tcanonical_key").map_err(|e|e.to_string())?;
    for (p,f) in fibers.iter().enumerate() {
        for &i in f {writeln!(witness,"{i}\t{p}\t{}\t{}\t{:?}",integers[i],integers[i]>=0,graph.states[i]).map_err(|e|e.to_string())?;}
    }
    let mut edges=create("EDGES.tsv")?;
    writeln!(edges,"from\tmove\tto\tsame_color").map_err(|e|e.to_string())?;
    let mut conflicts=[0usize;6];
    for (i,row) in graph.transitions.iter().enumerate() {
        for (mv,&j) in row.iter().enumerate() {
            let same=(integers[i]>=0)==(integers[j as usize]>=0);conflicts[mv]+=usize::from(same);
            writeln!(edges,"{i}\t{mv}\t{j}\t{same}").map_err(|e|e.to_string())?;
        }
    }
    let mut summary=create("SUMMARY.txt")?;
    writeln!(summary,"experiment = exceptional-position spectral decomposition\nsource = {}\nvertices = {n}\nfiber sizes = {:?}\nposition eigenvalues = {factor:?}\nwithin-position numerical minimum = {:.14}\nwithin-position numerical maximum = {:.14}\ndecomposition maximum numerical error = {decomposition_error:e}\nnegative numerical residual = {residual:e}\nnegative Rayleigh numerator = {num}\nnegative Rayleigh denominator = {den}\nexact lambda_min upper bound = {num}/{den}\nexact normalized absolute gap upper bound = {}/{}\nnormalized absolute gap upper bound decimal = {:.14}\ninteger fiber sums = [0, 0, 0, 0]\nsame-color directed arcs by move = {conflicts:?}\ntotal directed arcs = {}\nAll edges and inverse edges checked. Rational bound is exact; eigenvalues and residual are floating point. The sign cut is a witness, not a proven maximum cut.\nsearch complete = true",certificate.display(),fibers.iter().map(Vec::len).collect::<Vec<_>>(),internal.last().unwrap(),internal[0],6*den+num,6*den,1.0+rayleigh/6.0,6*n).map_err(|e|e.to_string())?;
    summary.sync_all().map_err(|e|e.to_string())?;
    Ok(())
}

#[derive(Clone, Copy)]
struct NullDatum { replicate: usize, vertices: usize, components: usize, loops: usize, gap: f64, w1: f64 }

#[derive(Default)]
struct NullClass {
    members: Vec<(String, String, i32)>,
    vertices: usize,
    observed: Option<NullDatum>,
    controls: BTreeMap<String, Vec<NullDatum>>,
}

#[derive(Default)]
struct NullSummary { classes: usize, obs_gap: f64, null_gap: f64, delta: f64, below: usize, above: usize, disconnected: usize, replicates: usize, w1: f64, w1_classes: usize }

fn analyze_nulls(args: &[String]) -> Result<(), String> {
    let inputs = get_arg(args, "--inputs", None)?;
    let out = PathBuf::from(get_arg(args, "--out", None)?);
    let mut classes = BTreeMap::<(String, usize), NullClass>::new();
    let mut seen_runs = HashSet::new();
    for raw in inputs.split(',') {
        let dir = Path::new(raw);
        let run = dir.file_name().and_then(|s| s.to_str()).ok_or("input needs a directory name")?.to_string();
        if !seen_runs.insert(run.clone()) { return Err(format!("duplicate input run name: {run}")); }
        let members = fs::read_to_string(dir.join("CLASS_MEMBERS.tsv")).map_err(|e| e.to_string())?;
        for (line_no, line) in members.lines().enumerate().skip(1) {
            let cols = line.split('\t').collect::<Vec<_>>();
            if cols.len() != 7 { return Err(format!("{run} CLASS_MEMBERS line {} has {} fields", line_no + 1, cols.len())); }
            let id = cols[1].parse::<usize>().map_err(|e| format!("{run} class id: {e}"))?;
            let genus = cols[4].parse::<i32>().map_err(|e| format!("{run} genus: {e}"))?;
            let vertices = cols[5].parse::<usize>().map_err(|e| format!("{run} vertices: {e}"))?;
            let class = classes.entry((run.clone(), id)).or_default();
            if class.vertices != 0 && class.vertices != vertices { return Err(format!("{run} class {id} inconsistent size")); }
            class.vertices = vertices;
            class.members.push((cols[0].to_string(), cols[2].to_string(), genus));
        }
        let samples = fs::read_to_string(dir.join("NULL_REPLICATES.tsv")).map_err(|e| e.to_string())?;
        for (line_no, line) in samples.lines().enumerate().skip(1) {
            let cols = line.split('\t').collect::<Vec<_>>();
            if cols.len() != 13 { return Err(format!("{run} NULL_REPLICATES line {} has {} fields", line_no + 1, cols.len())); }
            let id = cols[0].parse::<usize>().map_err(|e| format!("{run} class id: {e}"))?;
            let class = classes.get_mut(&(run.clone(), id)).ok_or_else(|| format!("{run} null class {id} lacks members"))?;
            let parse_usize = |i: usize| cols[i].parse::<usize>().map_err(|e| format!("{run} line {} field {i}: {e}", line_no + 1));
            let parse_float = |i: usize| cols[i].parse::<f64>().map_err(|e| format!("{run} line {} field {i}: {e}", line_no + 1));
            let row = NullDatum { replicate: parse_usize(3)?, vertices: parse_usize(4)?, components: parse_usize(5)?, loops: parse_usize(6)?, gap: parse_float(8)?, w1: parse_float(12)? };
            if row.vertices != class.vertices || !row.gap.is_finite() || row.components == 0 { return Err(format!("{run} class {id} invalid null row {}", line_no + 1)); }
            if cols[2] == "observed" {
                if row.replicate != 0 || row.components != 1 || class.observed.replace(row).is_some() { return Err(format!("{run} class {id} duplicate/invalid observation")); }
            } else { class.controls.entry(cols[2].to_string()).or_default().push(row); }
        }
    }
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let mut detail = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join("CLASS_COMPARISONS.tsv")).map_err(|e| e.to_string())?);
    let mut total = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join("MODEL_SUMMARY.tsv")).map_err(|e| e.to_string())?);
    writeln!(detail, "run\tclass_id\trepresentative\tgroups\tgenera\tmembers\tvertices\tobserved_loops\tmodel\treplicates\tobserved_gap\tmean_null_gap\tgap_minus_null_mean\tnull_gap_q05\tnull_gap_q95\tlower_tail_rank\tupper_tail_rank\tdisconnected_replicates\tmean_wasserstein_if_exact").map_err(|e| e.to_string())?;
    writeln!(total, "model\tclasses\tmean_observed_gap\tclass_mean_null_gap\tmean_gap_minus_null\tbelow_null_5pct\tabove_null_95pct\tdisconnected_replicates\ttotal_replicates\tclasses_with_exact_wasserstein\tclass_mean_wasserstein_if_exact").map_err(|e| e.to_string())?;
    let mut summaries = BTreeMap::<String, NullSummary>::new();
    let mut mixed_genus = 0usize;
    for ((run, id), class) in &mut classes {
        let observed = class.observed.ok_or_else(|| format!("{run} class {id} lacks observation"))?;
        let genera = class.members.iter().map(|x| x.2).collect::<BTreeSet<_>>();
        let groups = class.members.iter().map(|x| x.1.as_str()).collect::<BTreeSet<_>>();
        if genera.len() > 1 { mixed_genus += 1; }
        if class.controls.is_empty() { return Err(format!("{run} class {id} lacks controls")); }
        for (model, rows) in &mut class.controls {
            rows.sort_by_key(|r| r.replicate);
            if rows.iter().enumerate().any(|(i, r)| r.replicate != i) { return Err(format!("{run} class {id} {model} missing/duplicate replicate")); }
            if matches!(model.as_str(), "loop_matched" | "cycle_matched") && rows.iter().any(|r| r.loops != observed.loops) {
                return Err(format!("{run} class {id} {model} changed loop count"));
            }
            let n = rows.len();
            if n == 0 { return Err(format!("{run} class {id} {model} has no replicates")); }
            let mut gaps = rows.iter().map(|r| r.gap).collect::<Vec<_>>();
            gaps.sort_by(|a,b| a.total_cmp(b));
            let mean = gaps.iter().sum::<f64>() / n as f64;
            let lower = (1 + gaps.iter().filter(|&&g| g <= observed.gap + 1e-10).count()) as f64 / (n + 1) as f64;
            let upper = (1 + gaps.iter().filter(|&&g| g >= observed.gap - 1e-10).count()) as f64 / (n + 1) as f64;
            let disconnected = rows.iter().filter(|r| r.components != 1).count();
            let finite_w1 = rows.iter().filter_map(|r| r.w1.is_finite().then_some(r.w1)).collect::<Vec<_>>();
            let mean_w1 = if finite_w1.len() == n { finite_w1.iter().sum::<f64>() / n as f64 } else { f64::NAN };
            writeln!(detail, "{run}\t{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{model}\t{n}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{disconnected}\t{}",
                class.members[0].0, groups.iter().copied().collect::<Vec<_>>().join(","), genera.iter().map(i32::to_string).collect::<Vec<_>>().join(","),
                class.members.len(), class.vertices, observed.loops, format_float(observed.gap), format_float(mean), format_float(observed.gap - mean),
                format_float(gaps[n * 5 / 100]), format_float(gaps[n * 95 / 100]), format_float(lower), format_float(upper), format_float(mean_w1)).map_err(|e| e.to_string())?;
            let summary = summaries.entry(model.clone()).or_default();
            summary.classes += 1; summary.obs_gap += observed.gap; summary.null_gap += mean; summary.delta += observed.gap - mean;
            summary.below += usize::from(lower <= 0.05); summary.above += usize::from(upper <= 0.05);
            summary.disconnected += disconnected; summary.replicates += n;
            if mean_w1.is_finite() { summary.w1 += mean_w1; summary.w1_classes += 1; }
        }
    }
    for (model, s) in &summaries {
        writeln!(total, "{model}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.classes, format_float(s.obs_gap / s.classes as f64), format_float(s.null_gap / s.classes as f64),
            format_float(s.delta / s.classes as f64), s.below, s.above, s.disconnected, s.replicates, s.w1_classes,
            format_float(if s.w1_classes == 0 { f64::NAN } else { s.w1 / s.w1_classes as f64 })).map_err(|e| e.to_string())?;
    }
    let mut methods = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join("METHODS.txt")).map_err(|e| e.to_string())?);
    writeln!(methods, "classes={}\nmixed_genus_classes={}\nunit=exact move-labeled quotient graph class\nprimary_outcome=(6-lambda2)/6, set to zero for disconnected controls\nclass_weight=equal, regardless of case count or replicate count\ntail_rank=(1+number of controls on the indicated side)/(1+replicates), descriptive and unadjusted\nnull_quantiles=sorted samples at floor(0.05*B) and floor(0.95*B)\nexact_wasserstein=reported only when stored full eigenspectra are finite\nwarning=large-graph lambda2 values use finite-step Lanczos and need convergence checks\nwarning=shared graph classes can have multiple genera; these rows do not estimate a genus effect\nwarning=classes cover only the supplied saved catalogs/selection, not a random sample of all degrees or passports\ninputs={inputs}", classes.len(), mixed_genus).map_err(|e| e.to_string())?;
    eprintln!("null analysis complete: {} classes, {} mixed-genus classes, {}", classes.len(), mixed_genus, out.display());
    Ok(())
}

fn labeled_lift_graph(d: usize, seed: [Perm; 4], cap: usize) -> Option<SpectrumGraph> {
    let mut tuples = vec![seed];
    let mut seen = HashMap::from([([seed[0], seed[1], seed[2]], 0u32)]);
    let mut transitions = Vec::new();
    let mut head = 0;
    while head < tuples.len() {
        let mut row = [0u32; 6];
        for (mv, target) in row.iter_mut().enumerate() {
            let t = hurwitz_move_dynamic(d, tuples[head], mv);
            let key = [t[0], t[1], t[2]];
            *target = if let Some(&j) = seen.get(&key) { j } else {
                if tuples.len() >= cap { return None; }
                let j = tuples.len() as u32;
                seen.insert(key, j);
                tuples.push(t);
                j
            };
        }
        transitions.push(row);
        head += 1;
    }
    Some(SpectrumGraph { states: Vec::new(), transitions })
}

fn research(args: &[String]) -> Result<(), String> {
    let components = optional_arg(args, "--components");
    let certificates = optional_arg(args, "--certificates");
    if components.is_none() && certificates.is_none() { return Err("research needs --components and/or --certificates".into()); }
    let out = PathBuf::from(get_arg(args, "--out", None)?);
    let max_lift: usize = parse_num(get_arg(args, "--max-lift", Some("30000"))?, "--max-lift")?;
    let null_max: usize = parse_num(get_arg(args, "--null-max", Some("128"))?, "--null-max")?;
    let null_reps: usize = parse_num(get_arg(args, "--null-reps", Some("64"))?, "--null-reps")?;
    let full_max: usize = parse_num(get_arg(args, "--full-max", Some("32"))?, "--full-max")?;
    let standard_max: usize = parse_num(get_arg(args, "--standard-max", Some("256"))?, "--standard-max")?;
    let subset_max_k: usize = parse_num(get_arg(args, "--subset-max-k", Some("0"))?, "--subset-max-k")?;
    let subset_cap: usize = parse_num(get_arg(args, "--subset-cap", Some("1000000"))?, "--subset-cap")?;
    let ordered_max_k: usize = parse_num(get_arg(args, "--ordered-max-k", Some("0"))?, "--ordered-max-k")?;
    if ordered_max_k > 6 { return Err("ordered-max-k must be at most 6 for packed tuple encoding".into()); }
    let steps: usize = parse_num(get_arg(args, "--lanczos-steps", Some("128"))?, "--lanczos-steps")?;
    let trace_probes: usize = parse_num(get_arg(args, "--trace-probes", Some("8"))?, "--trace-probes")?;
    if max_lift == 0 || full_max == 0 || steps == 0 || (null_reps > 0 && trace_probes == 0) {
        return Err("max-lift, full-max, lanczos-steps, and trace-probes must be positive".into());
    }
    let only = optional_arg(args, "--only").map(|s| s.split(',').map(str::to_string).collect::<HashSet<_>>());
    let mut cases = Vec::new();
    let mut roots = HashMap::new();
    if let Some(paths) = components {
        for raw in paths.split(',') {
            for case in parse_components_file(Path::new(raw))? { add_spectrum_case(&mut cases, &mut roots, case); }
        }
    }
    if let Some(paths) = certificates {
        let mut files = Vec::new();
        for raw in paths.split(',') { collect_certificate_paths(Path::new(raw), &mut files)?; }
        files.sort();
        files.dedup();
        for file in files { add_spectrum_case(&mut cases, &mut roots, parse_certificate_case(&file)?); }
    }
    if let Some(only) = &only { cases.retain(|case| only.contains(&case.id)); }
    if cases.is_empty() { return Err("research inputs selected no cases".into()); }
    fs::create_dir_all(&out).map_err(|e| e.to_string())?;
    let create = |name: &str| -> Result<BufWriter<File>, String> {
        Ok(BufWriter::new(OpenOptions::new().write(true).create_new(true).open(out.join(name)).map_err(|e| format!("{}: {e}", out.join(name).display()))?))
    };
    let mut case_file = create("CASES.tsv")?;
    let mut class_file = create("CLASS_MEMBERS.tsv")?;
    let mut null_file = create("NULL_REPLICATES.tsv")?;
    let mut voltage_file = create("VOLTAGES.tsv")?;
    let mut state_file = create("STATES.tsv")?;
    let mut subset_file = create("SUBSET_SPECTRA.tsv")?;
    writeln!(case_file, "case_id\tgroup\tdegree\tgenus\tpassport\tquotient_vertices\tholonomy_order\tlabeled_vertices_expected\tloop_holonomy_order\tholonomy_generators\tholonomy_even\tholonomy_point_orbits\tlabeled_vertices_computed\tquotient_lambda2\tlabeled_lambda2\tlabeled_lambda2_check\tlabeled_method\tstandard_point_lambda_max").map_err(|e| e.to_string())?;
    writeln!(class_file, "case_id\tclass_id\tgroup\tdegree\tgenus\tvertices\tpassport").map_err(|e| e.to_string())?;
    writeln!(null_file, "class_id\trepresentative\tkind\treplicate\tvertices\tcomponents\tloops\tlambda2\trelaxation_gap\tmoment3\tmoment5\tmoment7\twasserstein_nontrivial").map_err(|e| e.to_string())?;
    writeln!(voltage_file, "case_id\tfrom\tmove\tto\tedge_conjugator\tbased_holonomy\tbased_holonomy_even").map_err(|e| e.to_string())?;
    writeln!(state_file, "case_id\tvertex\tfirst\tsecond\tthird").map_err(|e| e.to_string())?;
    writeln!(subset_file, "case_id\tk\tvertices\tstatus\tcomponents\tlambda2\tlambda2_check\tmethod").map_err(|e| e.to_string())?;
    let mut classes = HashMap::<String, usize>::new();
    let started = Instant::now();
    for (number, case) in cases.iter().enumerate() {
        let graph = match &case.states { Some(states) => graph_from_states(case.d, states.clone())?, None => enumerate_spectrum_graph(case.d, case.seed)? };
        for (i, key) in graph.states.iter().enumerate() {
            writeln!(state_file, "{}\t{}\t{:032x}\t{:032x}\t{:032x}", case.id, i, key[0], key[1], key[2]).map_err(|e| e.to_string())?;
        }
        state_file.flush().map_err(|e| e.to_string())?;
        if let Some(expected) = case.input_states {
            if graph.transitions.len() > expected { return Err(format!("{}: quotient graph exceeds input state count", case.id)); }
        }
        let profile = lift_profile(case.d, &graph)?;
        for &(from, mv, to, edge, based) in &profile.edges {
            writeln!(voltage_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}", case.id, from, mv, to,
                edge.cycle_string(case.d), based.cycle_string(case.d), based.is_even(case.d)).map_err(|e| e.to_string())?;
        }
        voltage_file.flush().map_err(|e| e.to_string())?;
        let n = graph.transitions.len();
        let expected_lift = profile.holonomy_order * n as u128;
        let quotient = eigen_for_graph(&graph, full_max, steps, case_seed(case, 0x1111));
        let standard_lambda = if n <= standard_max / (case.d - 1) {
            standard_point_lambda_max(case.d, n, &profile.edges)?
        } else { f64::NAN };
        let mut lift_n = 0;
        let mut lift_lambda2 = f64::NAN;
        let mut lift_check = f64::NAN;
        let mut lift_method = "skipped-above-cap".to_string();
        if expected_lift <= max_lift as u128 {
            let lift = labeled_lift_graph(case.d, unpack(graph.states[0], case.d), max_lift).ok_or_else(|| format!("{}: lift exceeded exact holonomy size", case.id))?;
            lift_n = lift.transitions.len();
            if lift_n as u128 != expected_lift { return Err(format!("{}: lift size {} disagrees with holonomy prediction {expected_lift}", case.id, lift_n)); }
            let eig = eigen_for_graph(&lift, full_max, steps, case_seed(case, 0x2222));
            lift_lambda2 = eig.lambda2;
            lift_method = eig.method;
            if lift_n > full_max { lift_check = lanczos_extremes(&lift, steps * 2, case_seed(case, 0x2222)).0; }
        }
        writeln!(case_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            case.id, case.ambient.short(), case.d, genus_of(&case.seed, case.d), passport_string(&case.seed, case.d), n,
            profile.holonomy_order, expected_lift, profile.loop_order, profile.generators, profile.all_even, profile.point_orbits,
            lift_n, format_float(quotient.lambda2), format_float(lift_lambda2), format_float(lift_check), lift_method,
            format_float(standard_lambda)).map_err(|e| e.to_string())?;
        case_file.flush().map_err(|e| e.to_string())?;
        for k in 1..=subset_max_k.min(case.d / 2) {
            let vertices = n.checked_mul(binomial(case.d, k)).ok_or("subset action size overflow")?;
            if vertices > subset_cap || vertices > u32::MAX as usize {
                writeln!(subset_file, "{}\t{}\t{}\tcap\t\tnan\tnan\t", case.id, k, vertices).map_err(|e| e.to_string())?;
                continue;
            }
            let subset = subset_action_graph(case.d, &graph, &profile.edges, k);
            let components = graph_components(&subset);
            let eig = eigen_for_graph(&subset, full_max, steps, case_seed(case, 0x5555 + k as u64));
            let check = if vertices > full_max { lanczos_extremes(&subset, steps * 2, case_seed(case, 0x5555 + k as u64)).0 } else { f64::NAN };
            writeln!(subset_file, "{}\t{}\t{}\tcomputed\t{}\t{}\t{}\t{}", case.id, k, vertices, components, format_float(eig.lambda2), format_float(check), eig.method).map_err(|e| e.to_string())?;
            subset_file.flush().map_err(|e| e.to_string())?;
            eprintln!("research subset case={} k={} vertices={} lambda2={} elapsed={:.1}s", case.id, k, vertices, format_float(eig.lambda2), started.elapsed().as_secs_f64());
        }
        for k in 1..=ordered_max_k.min(case.d) {
            let width = (0..k).try_fold(1usize, |value, i| value.checked_mul(case.d - i));
            let vertices = width.and_then(|width| width.checked_mul(n)).unwrap_or(usize::MAX);
            if vertices > subset_cap || vertices > u32::MAX as usize {
                writeln!(subset_file, "{}\tordered-{}\t{}\tcap\t\tnan\tnan\t", case.id, k, vertices).map_err(|e| e.to_string())?;
                continue;
            }
            let ordered = ordered_action_graph(case.d, &graph, &profile.edges, k);
            let components = graph_components(&ordered);
            let eig = eigen_for_graph(&ordered, full_max, steps, case_seed(case, 0x6666 + k as u64));
            let check = if vertices > full_max { lanczos_extremes(&ordered, steps * 2, case_seed(case, 0x6666 + k as u64)).0 } else { f64::NAN };
            writeln!(subset_file, "{}\tordered-{}\t{}\tcomputed\t{}\t{}\t{}\t{}", case.id, k, vertices, components, format_float(eig.lambda2), format_float(check), eig.method).map_err(|e| e.to_string())?;
            subset_file.flush().map_err(|e| e.to_string())?;
            eprintln!("research ordered case={} k={} vertices={} lambda2={} elapsed={:.1}s", case.id, k, vertices, format_float(eig.lambda2), started.elapsed().as_secs_f64());
        }
        if null_reps > 0 && n <= null_max {
            let code = canonical_labeled_code(&graph, null_max);
            let next_id = classes.len() + 1;
            let (class_id, new_class) = match classes.entry(code) {
                std::collections::hash_map::Entry::Occupied(entry) => (*entry.get(), false),
                std::collections::hash_map::Entry::Vacant(entry) => { entry.insert(next_id); (next_id, true) },
            };
            writeln!(class_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}", case.id, class_id, case.ambient.short(), case.d, genus_of(&case.seed, case.d), n, passport_string(&case.seed, case.d)).map_err(|e| e.to_string())?;
            if new_class {
                // Only the first representative of each exact move-labeled class is sampled.
                let mut write_row = |kind: &str, replicate: usize, sample: &SpectrumGraph| -> Result<(), String> {
                    let eig = if kind == "observed" { quotient.clone() } else { eigen_for_graph(sample, full_max, steps, case_seed(case, 0x3333 + replicate as u64)) };
                    let moments = stochastic_trace_moments(sample, &eig, 7, trace_probes, case_seed(case, 0x4444 + replicate as u64));
                    let (loops, _) = graph_shape(sample);
                    let gap = if graph_components(sample) == 1 { (6.0 - eig.lambda2) / 6.0 } else { 0.0 };
                    writeln!(null_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        class_id, case.id, kind, replicate, n, graph_components(sample), loops, format_float(eig.lambda2), format_float(gap),
                        format_float(moments[2].1), format_float(moments[4].1), format_float(moments[6].1),
                        format_float(spectral_wasserstein(&quotient.values, &eig.values))).map_err(|e| e.to_string())
                };
                write_row("observed", 0, &graph)?;
                for replicate in 0..null_reps {
                    for (kind, salt) in [("random_6_regular", 0x5100), ("random_schreier", 0x5200), ("cycle_matched", 0x5300), ("loop_matched", 0x5400)] {
                        let mut rng = SplitMix64::new(case_seed(case, salt + replicate as u64));
                        let sample = match kind {
                            "random_6_regular" => random_regular_graph(n, 6, 20, &mut rng),
                            "random_schreier" => Some(random_schreier_graph(n, &mut rng)),
                            "cycle_matched" => Some(cycle_matched_schreier_graph(&graph, &mut rng)),
                            _ => Some(loop_matched_schreier_graph(&graph, &mut rng)),
                        };
                        if let Some(sample) = sample {
                            if kind == "loop_matched" && graph.transitions.iter().enumerate().any(|(i, row)| (0..6).any(|mv| (row[mv] as usize == i) != (sample.transitions[i][mv] as usize == i))) {
                                return Err("loop-matched null changed a loop position".into());
                            }
                            write_row(kind, replicate, &sample)?;
                        }
                    }
                }
                null_file.flush().map_err(|e| e.to_string())?;
            }
        }
        class_file.flush().map_err(|e| e.to_string())?;
        eprintln!("research progress {}/{} case={} quotient={} lift={} classes={} elapsed={:.1}s", number + 1, cases.len(), case.id, n, expected_lift, classes.len(), started.elapsed().as_secs_f64());
    }
    let mut summary = create("SUMMARY.txt")?;
    writeln!(summary, "cases = {}\nexact move-labeled null classes = {}\nnull replicates per class = {}\nmaximum labeled lift size = {}\nmaximum null quotient size = {}\nfull eigenspectrum threshold = {}\nstandard point block dimension threshold = {}\nmaximum subset size = {}\nmaximum ordered tuple size = {}\naction graph vertex cap = {}\nlanczos steps = {}\ntrace probes = {}\nelapsed seconds = {:.3}", cases.len(), classes.len(), null_reps, max_lift, null_max, full_max, standard_max, subset_max_k, ordered_max_k, subset_cap, steps, trace_probes, started.elapsed().as_secs_f64()).map_err(|e| e.to_string())?;
    eprintln!("research complete: {}", out.display());
    Ok(())
}

// Inspect ALL states, not merely the construction's starting tuple. This is
// a structural fingerprint, not proof that two components share a mechanism.
fn component_profile(d: usize, states: &[Key]) -> ([usize;3], BTreeMap<u8,usize>, BTreeMap<Vec<usize>,usize>) {
    let mut fixed=[0;3];let mut masks=BTreeMap::new();let mut partitions=BTreeMap::new();
    for &key in states {
        let t=unpack(key,d);let mut mask=0u8;
        for i in 0..3 {
            if canonical(&hurwitz_move_dynamic(d,t,i),d)==key {fixed[i]+=1;mask|=1<<i;}
        }
        *masks.entry(mask).or_insert(0)+=1;
        let mut sizes:Vec<_>=point_orbits(d,&t[..2]).iter().map(Vec::len).collect();
        sizes.sort_unstable();*partitions.entry(sizes).or_insert(0)+=1;
    }
    (fixed,masks,partitions)
}

fn certificate_profile(args: &[String]) -> Result<(),String> {
    let path=get_arg(args,"--certificate",None)?;verify(Path::new(&path))?;
    let case=parse_certificate_case(Path::new(&path))?;
    let states=case.states.as_ref().ok_or("missing closed states")?;
    let (fixed,masks,partitions)=component_profile(case.d,states);
    let no_fixed=*masks.get(&0).unwrap_or(&0);
    let mut defects=BTreeMap::new();for &key in states {*defects.entry(defect_of(&unpack(key,case.d),case.d).0).or_insert(0usize)+=1;}
    let out=get_arg(args,"--out",None)?;
    let mut file=OpenOptions::new().write(true).create_new(true).open(&out).map_err(|e|format!("{out}: {e}"))?;
    writeln!(file,"group\tdegree\tgenus\tquotient_states\tpassport\tsigma1_fixed_states\tsigma2_fixed_states\tsigma3_fixed_states\tany_braid_fixed_states\tno_braid_fixed_states\twhole_component_no_braid_fixed\tloop_mask_histogram\tfirst_pair_partition_histogram\tconstant_defect\tdefect_histogram").map_err(|e|e.to_string())?;
    writeln!(file,"{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:?}\t{}\t{:?}",case.ambient.short(),case.d,genus_of(&case.seed,case.d),states.len(),passport_string(&case.seed,case.d),fixed[0],fixed[1],fixed[2],states.len()-no_fixed,no_fixed,no_fixed==states.len(),masks,partitions,defects.len()==1,defects).map_err(|e|e.to_string())?;
    if has_flag(args,"--pure-stability") {
        // Keep the existing TSV schema intact for historical readers.
        let mut pure=OpenOptions::new().write(true).create_new(true).open(format!("{out}.pure.txt")).map_err(|e|e.to_string())?;
        pure_stability::write_profile(&mut pure,case.d,states).map_err(|e|e.to_string())?;
        pure.sync_all().map_err(|e|e.to_string())?;
    }
    file.sync_all().map_err(|e|e.to_string())
}

pub(super) fn run(args: &[String]) -> Result<(), String> {
    if args[1]=="defect-mission" {return defect_mission::run(args);}
    if args[1]=="certificate-profile" {return certificate_profile(args);}
    if args[1]=="defect-census" {return defect_census(args);}
    if args[1]=="hexagon-search" {return hexagon_search(args);}
    if args[1]=="braid-library" {
        if has_flag(args,"--state-actions") && has_flag(args,"--independent-actions") {return Err("choose state-actions OR independent-actions".into());}
        if has_flag(args,"--fiber-actions") {return braid_frames::fiber_actions(args);}
        if has_flag(args,"--state-actions") {return braid_frames::state_actions(args);}
        if has_flag(args,"--independent-actions") {return braid_frames::generate(args);}
        return braid_library(args);
    }
    if args[1]=="rigid-glue" {return rigid_glue(args);}
    if args[1]=="braid-restrict" {return braid_restrict(args);}
    if args[1]=="verify" {return verify(Path::new(&get_arg(args,"--certificate",None)?));}
    if args[1]=="spectrum" {return spectrum(args);}
    if args[1]=="research" {return research(args);}
    if args[1]=="rayleigh" {return rayleigh(args);}
    if args[1]=="position-spectrum" {return position_spectrum(args);}
    if args[1]=="braid-symmetry" {return braid_symmetry(args);}
    if args[1]=="mechanism-audit" {return mechanism_audit(args);}
    if args[1]=="s18-fiber-census" {return s18_fiber_census(args);}
    if args[1]=="analyze-nulls" {return analyze_nulls(args);}
    let d: usize = parse_num(get_arg(args, "--d", None)?, "--d")?;
    if !(4..=MAX_D).contains(&d) { return Err(format!("require 4 <= d <= {MAX_D}")); }
    let ambient = Ambient::parse(&get_arg(args, "--group", None)?)?;
    let cap: usize = parse_num(get_arg(args, "--cap", Some(if args[1] == "orbit-map" { "1024" } else { "200000" }))?, "--cap")?;
    let seconds: u64 = parse_num(get_arg(args, "--seconds", Some("60"))?, "--seconds")?;
    if cap == 0 || seconds == 0 { return Err("cap and seconds must be positive".into()); }
    if args[1] == "orbit-map" { return orbit_map(args, d, ambient, cap, seconds); }
    if args[1]=="smart-hunt" || args[1]=="passport-exact" {return smart_hunt(args,d,ambient,cap,seconds);}
    let t = parse_tuple(&get_arg(args, "--tuple", None)?, d)?;
    validate_tuple(d, ambient, &t)?;
    let result = explore(d, t, cap, Instant::now()+Duration::from_secs(seconds));
    println!("status={} quotient_states={} weak={} orders={:?}", result.status, result.states.len(), result.weak_witness.is_some(), result.orders);
    if result.status.starts_with("certified_") {
        let path = PathBuf::from(get_arg(args, "--out", None)?);
        write_certificate(&path, d, ambient, &t, &result).map_err(|e| e.to_string())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
