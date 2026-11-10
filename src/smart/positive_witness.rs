//! Exact positive words shared between completed jobs; never cache nonescape.
use super::*;

pub(super) struct Cache {
    words: Mutex<HashMap<Key, Vec<u8>>>,
    output: Mutex<fs::File>,
    hits: AtomicU64,
    capacity: usize,
}

fn check(d: usize, seed: [Perm; 4], word: &[u8]) -> bool {
    if !tuple_product_is_identity(d, &seed) || !generators_transitive(d, &seed[..3]) { return false; }
    let mut t = seed;
    for &mv in word {
        if mv >= 6 { return false; }
        t = hurwitz_move_dynamic(d, t, mv as usize);
    }
    pair_transitive(d, t[0], t[1])
}

pub(super) fn short_witness(d:usize,seed:[Perm;4])->Option<&'static [u8]> {
    // Common exported escape words. Each rejection is a literal positive
    // test on this tuple, independent of how these words were selected.
    // No full-group generation or component classification is inferred.
    const WORDS:&[&[u8]]=&[&[1,1],&[4,4],&[1,1,1],&[4,4,4],
        &[1,3,1],&[1,0,0,0,2,1]];
    for &word in WORDS {
        let mut t=seed;
        for &mv in word {t=hurwitz_move_dynamic(d,t,mv as usize);}
        if pair_transitive(d,t[0],t[1]) {return Some(word);}
    }
    None
}

impl Cache {
    pub(super) fn new(args: &[String], d: usize, out: &Path) -> Result<Self, String> {
        let capacity = 50000;
        let mut words = HashMap::new();
        if let Some(path) = optional_arg(args, "--positive-witnesses") {
            let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
            if text.lines().next() != Some("degree\tmoves\tseed") { return Err("invalid positive-witness header".into()); }
            for (i, line) in text.lines().skip(1).enumerate() {
                let fields = line.split('\t').collect::<Vec<_>>();
                if fields.len() != 3 || fields[0] != d.to_string() { return Err(format!("invalid positive-witness row {}", i + 2)); }
                let word = if fields[1].is_empty() { Vec::new() } else {
                    fields[1].split(',').map(|s| s.parse::<u8>().map_err(|_| "invalid braid move".to_string())).collect::<Result<Vec<_>, _>>()?
                };
                if word.len() > 4096 { return Err("positive word exceeds validation bound".into()); }
                let seed = parse_tuple(fields[2], d)?;
                if !check(d, seed, &word) { return Err(format!("false positive witness in row {}", i + 2)); }
                if words.len() < capacity { words.entry(canonical(&seed, d)).or_insert(word); }
            }
        }
        let mut output = OpenOptions::new().write(true).create_new(true).open(out.join("POSITIVE_WITNESSES.tsv")).map_err(|e| e.to_string())?;
        writeln!(output, "degree\tmoves\tseed").map_err(|e| e.to_string())?;
        Ok(Self { words: Mutex::new(words), output: Mutex::new(output), hits: AtomicU64::new(0), capacity })
    }

    pub(super) fn contains(&self, key: Key) -> bool {
        let hit = self.words.lock().unwrap().contains_key(&key);
        if hit { self.hits.fetch_add(1, Ordering::Relaxed); }
        hit
    }

    pub(super) fn insert(&self, d: usize, seed: [Perm; 4], word: &[u8]) -> io::Result<()> {
        if word.len() > 4096 { return Ok(()); }
        if !check(d, seed, word) { return Err(io::Error::other("invalid reconstructed positive braid word")); }
        let key = canonical(&seed, d);
        let mut words = self.words.lock().unwrap();
        if words.contains_key(&key) || words.len() >= self.capacity { return Ok(()); }
        words.insert(key, word.to_vec());
        let mut output = self.output.lock().unwrap();
        writeln!(output, "{d}\t{}\t{}", word.iter().map(u8::to_string).collect::<Vec<_>>().join(","), tuple_string_dynamic(d, &unpack(key, d)))?;
        Ok(())
    }

    pub(super) fn summary(&self, out: &mut fs::File) -> io::Result<()> {
        self.output.lock().unwrap().sync_all()?;
        writeln!(out, "verified positive words retained = {}\npositive witness cache hits = {}\npositive witness capacity = {}\nPositive words certify transitive first pairs; cache misses/caps certify nothing.", self.words.lock().unwrap().len(), self.hits.load(Ordering::Relaxed), self.capacity)
    }
}

#[cfg(test)] mod tests {
    use super::*;
    #[test] fn short_words_retain_every_known_weak_state() {
        for entry in fs::read_dir("evidence/certificates").unwrap() {
            let path=entry.unwrap().path();let name=path.file_name().unwrap().to_string_lossy();
            if !name.starts_with("CERTIFIED_") || !name.ends_with("weak.txt") {continue;}
            let case=parse_certificate_case(&path).unwrap();
            for key in case.states.unwrap() {assert!(short_witness(case.d,unpack(key,case.d)).is_none(),"{name}");}
        }
        let t=parse_tuple("[(2 3)(5 7)(6 8)(9 10), (3 4)(5 6)(7 9)(8 10), (1 6 10)(2 9)(3 8)(4 5), (1 10 7)(2 5)(3 9)(4 8)]",10).unwrap();
        assert!(check(10,t,short_witness(10,t).unwrap()));
    }
    #[test] fn traced_probe_preserves_decisions_and_rng() {
        let t=parse_tuple("[(2 3)(5 7)(6 8)(9 10), (3 4)(5 6)(7 9)(8 10), (1 6 10)(2 9)(3 8)(4 5), (1 10 7)(2 5)(3 9)(4 8)]",10).unwrap();
        for weak_only in [false,true] {
            for seed in 0..16 {
                let (mut plain_rng,mut traced_rng)=(SplitMix64::new(seed),SplitMix64::new(seed));
                let plain=probe(10,t,&mut plain_rng,64,weak_only);
                let mut word=vec![6];
                let traced=probe_with_word(10,t,&mut traced_rng,64,weak_only,Some(&mut word));
                assert_eq!(plain,traced);
                assert_eq!(plain_rng.below(10000),traced_rng.below(10000));
                assert!(word.len()<=64 && word.iter().all(|&m|m<6));
                if traced.0 {assert!(check(10,t,&word));}
            }
        }
    }
    #[test] fn closure_word_is_literal_and_conjugation_invariant() {
        // Use a guaranteed product-one transitive tuple, independent of notation convention.
        let mut a = Perm::identity(4); a.img[..4].copy_from_slice(&[1,0,3,2]);
        let mut b = Perm::identity(4); b.img[..4].copy_from_slice(&[0,2,1,3]);
        let t = [a,b,Perm::identity(4),product_one_fourth(4,a,b,Perm::identity(4))];
        assert!(check(4,t,&[])); assert!(!check(4,t,&[6]));
        let result = explore_goal(4,t,100,Instant::now()+Duration::from_secs(2),true);
        let word = result.weak_word.unwrap(); assert!(check(4,t,&word));
        let q = perm_of_type(&[4],4,None);
        assert!(check(4,t.map(|p|p.conjugate_by(q,4)),&word));
        let t=parse_tuple("[(2 3)(5 7)(6 8)(9 10), (3 4)(5 6)(7 9)(8 10), (1 6 10)(2 9)(3 8)(4 5), (1 10 7)(2 5)(3 9)(4 8)]",10).unwrap();
        assert!(!check(10,t,&[]));assert!(check(10,t,&[1,3,1]));
        let result=explore_goal(10,t,100000,Instant::now()+Duration::from_secs(3),true);
        let word=result.weak_word.unwrap();assert!(!word.is_empty());assert!(check(10,t,&word));
        let q=perm_of_type(&[10],10,None);assert!(check(10,t.map(|p|p.conjugate_by(q,10)),&word));
    }
}
