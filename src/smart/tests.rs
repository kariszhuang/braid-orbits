use super::*;
const A8: &str = "[(6 7 8), (3 4 6)(5 7 8), (1 3)(2 7)(4 6)(5 8), (1 3 4)(2 7 6)]";

#[test]
fn closed_component_profiles_check_all_three_generators() {
    for name in ["CERTIFIED_A25_weak.txt","CERTIFIED_A25_order3_weak.txt","CERTIFIED_A25_order6_genus12_weak.txt","CERTIFIED_A25_order12_genus12_weak.txt","CERTIFIED_S18_genus1_weak.txt"] {
        let case=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();
        let states=case.states.unwrap();let (fixed,masks,partitions)=component_profile(case.d,&states);
        assert_eq!(masks.values().sum::<usize>(),states.len());
        assert_eq!(partitions.values().sum::<usize>(),states.len());
        assert!(fixed.iter().any(|&n|n>0));
        assert!(partitions.keys().all(|p|p.len()>1 && p.iter().sum::<usize>()==case.d));
        assert!(genus_of(&case.seed,case.d)<=case.d as i64-5);
        for (i,&n) in fixed.iter().enumerate() {
            assert_eq!(n,masks.iter().filter(|(mask,_)|**mask&(1<<i)!=0).map(|(_,n)|n).sum());
        }
    }
}

#[test]
fn defect_decomposes_genus_on_known_component() {
    let case = parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt")).unwrap();
    let graph = graph_from_states(case.d, case.states.unwrap()).unwrap();
    for state in graph.states {
        let t = unpack(state, 18);
        let (defect, left, right) = defect_of(&t, 18);
        assert!(defect >= 0 && left >= 0 && right >= 0);
        assert_eq!(defect + left + right, genus_of(&t, 18));
    }
    let identity_case = parse_tuple("[(), (2 3 4), (1 2 3 4), (1 4 2 3)]", 4).unwrap();
    assert_eq!(genus_of(&identity_case, 4), 1);
    assert!(enumerate_spectrum_graph(4, identity_case).unwrap().states.into_iter()
        .all(|state| defect_of(&unpack(state, 4), 4).0 == 0));
}

#[test]
fn hexagon_centralizer_enumeration_is_exact() {
    for k in 2..=3 {
        let d=6*k;let (a,b)=matching_pair(&vec![6;k]);
        assert_eq!(point_orbits(d,&[a,b]).len(),k);
        assert_eq!(cycle_type(a.compose(b,d),d),vec![3;2*k]);
        let mut seen=HashSet::new();
        let blocks=point_orbits(d,&[a,b]);
        assert!(enumerate_pair_centralizer(&blocks,a,b,0,0,Perm::identity(d),&mut |g| {
            assert!(seen.insert(g));
            assert_eq!(g.compose(a,d),a.compose(g,d));
            assert_eq!(g.compose(b,d),b.compose(g,d));
            Ok(true)
        }).unwrap());
        assert_eq!(seen.len(),6usize.pow(k as u32)*factorial_usize(k));
        if k==3 {
            let mut rng=SplitMix64::new(20260926);
            let h0=pair_braid_conjugator(d,a,b,&mut rng).unwrap();
            let mut types=BTreeMap::new();
            for g in seen {
                let h=h0.compose(g,d);
                *types.entry(cycle_type(h,d)).or_insert(0usize)+=1;
            }
            assert_eq!(types.values().sum::<usize>(),1296);
            assert_eq!(types.get(&vec![3;6]),Some(&73));
        }
    }
    for sizes in [vec![16,2],vec![14,4],vec![12,6],vec![10,8],vec![8,6,4]] {
        let d:usize=sizes.iter().sum();let (a,b)=matching_pair(&sizes);
        assert_eq!(point_orbits(d,&[a,b]).iter().map(Vec::len).collect::<Vec<_>>(),sizes);
        let mut count=0;
        let blocks=point_orbits(d,&[a,b]);
        assert!(enumerate_pair_centralizer(&blocks,a,b,0,0,Perm::identity(d),&mut |g| {
            count+=1;assert_eq!(g.compose(a,d),a.compose(g,d));assert_eq!(g.compose(b,d),b.compose(g,d));Ok(true)
        }).unwrap());
        assert_eq!(count,sizes.iter().product::<usize>());
    }
}

#[test]
fn new_a25_order_three_weak_failure_is_distinct_and_exact() {
    let cert=Path::new("evidence/certificates/CERTIFIED_A25_order3_weak.txt");
    verify(cert).unwrap();
    let text=fs::read_to_string(cert).unwrap();
    let t=parse_tuple(text.lines().find_map(|s|s.strip_prefix("seed=")).unwrap(),25).unwrap();
    assert_eq!(genus_of(&t,25),0);
    let mut class=vec![1;10];class.extend(vec![3;5]);
    for a in &t[..3] {assert_eq!(cycle_type(*a,25),class);}
    let old=include_str!("../../evidence/certificates/CERTIFIED_A25_weak.txt");
    let old_states:HashSet<_>=old.lines().filter(|l|l.starts_with("STATE\t")).collect();
    assert!(text.lines().filter(|l|l.starts_with("STATE\t")).all(|s|!old_states.contains(s)));
    for s in text.lines().filter_map(|s|s.strip_prefix("STATE\t")) {
        let q=unpack(parse_hex_key(&s.split('\t').collect::<Vec<_>>()).unwrap(),25);
        let mut parts=point_orbits(25,&q[..2]).iter().map(Vec::len).collect::<Vec<_>>();
        parts.sort();
        assert!(parts==vec![1,1,1,3,3,16] || parts==vec![10,15]);
    }
}

#[test]
fn a25_order_six_genus_twelve_failure_is_distinct_and_exact() {
    let cert=Path::new("evidence/certificates/CERTIFIED_A25_order6_genus12_weak.txt");
    verify(cert).unwrap();
    let text=fs::read_to_string(cert).unwrap();
    let t=parse_tuple(text.lines().find_map(|s|s.strip_prefix("seed=")).unwrap(),25).unwrap();
    assert_eq!(genus_of(&t,25),12);
    let mut class=vec![2;2];class.extend(vec![3;3]);class.extend(vec![6;2]);
    for a in &t[..3] {assert_eq!(cycle_type(*a,25),class);}
    let mut exceptional=vec![3;3];exceptional.extend(vec![4;4]);
    assert_eq!(cycle_type(t[3],25),exceptional);
    assert_eq!(text.lines().filter(|l|l.starts_with("STATE\t")).count(),4);
    let mut partitions=BTreeMap::new();
    for s in text.lines().filter_map(|s|s.strip_prefix("STATE\t")) {
        let q=unpack(parse_hex_key(&s.split('\t').collect::<Vec<_>>()).unwrap(),25);
        let mut parts=point_orbits(25,&q[..2]).iter().map(Vec::len).collect::<Vec<_>>();
        parts.sort();*partitions.entry(parts).or_insert(0usize)+=1;
    }
    assert_eq!(partitions,BTreeMap::from([(vec![5,20],2),(vec![10,15],2)]));
    for old in [include_str!("../../evidence/certificates/CERTIFIED_A25_weak.txt"),
                include_str!("../../evidence/certificates/CERTIFIED_A25_order3_weak.txt")] {
        let old_states:HashSet<_>=old.lines().filter(|l|l.starts_with("STATE\t")).collect();
        assert!(text.lines().filter(|l|l.starts_with("STATE\t")).all(|s|!old_states.contains(s)));
    }
}

#[test]
fn a25_order_twelve_one_state_failure_is_distinct_and_exact() {
    let path=Path::new("evidence/certificates/CERTIFIED_A25_order12_genus12_weak.txt");
    verify(path).unwrap();
    let case=parse_certificate_case(path).unwrap();let t=case.seed;
    assert_eq!(genus_of(&t,25),12);
    for p in t {assert_eq!(cycle_type(p,25),vec![3,3,3,4,4,4,4]);}
    assert_eq!(case.states.as_ref().unwrap().len(),1);
    let key=canonical(&t,25);
    for mv in 0..6 {assert_eq!(canonical(&hurwitz_move_dynamic(25,t,mv),25),key);}
    assert!(!pair_transitive(25,t[0],t[1]));
    for name in ["CERTIFIED_A25_weak.txt","CERTIFIED_A25_order3_weak.txt","CERTIFIED_A25_order6_genus12_weak.txt"] {
        let old=parse_certificate_case(&Path::new("evidence/certificates").join(name)).unwrap();
        assert_ne!(passport_string(&t,25),passport_string(&old.seed,25));
    }
}

#[test]
fn s18_genus_one_weak_certificate_is_exact() {
    let path=Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt");
    verify(path).unwrap();
    let text=fs::read_to_string(path).unwrap();
    let t=parse_tuple(text.lines().find_map(|s|s.strip_prefix("seed=")).unwrap(),18).unwrap();
    assert_eq!(genus_of(&t,18),1);
    for a in &t[..3] {assert_eq!(cycle_type(*a,18),vec![2;9]);}
    let mut exceptional=vec![1;6];exceptional.extend(vec![4;3]);
    assert_eq!(cycle_type(t[3],18),exceptional);
    assert!(tuple_conjugator(18,t,hurwitz_move_dynamic(18,t,0)).is_ok());
    assert!(tuple_conjugator(18,t,hurwitz_move_dynamic(18,t,1)).is_err());
    let h=tuple_conjugator(18,t,hurwitz_move_dynamic(18,t,0)).unwrap();
    assert_eq!(cycle_type(h,18),vec![3;6]);
    assert_eq!(point_orbits(18,&t[..2]).len(),3);
    assert_eq!(point_orbits(18,&[t[0],t[1],h]).len(),1);
    assert_eq!(point_orbits(18,&t[2..]).len(),3);
    let mut right=point_orbits(18,&t[2..]);
    for orbit in &mut right {orbit.sort_unstable();}
    right.sort();
    assert_eq!(right,vec![vec![0,3,7,11,12,16],vec![1,2,5,9,13,15],vec![4,6,8,10,14,17]]);
    assert_eq!(cycle_type(t[0].compose(t[1],18),18),vec![3;6]);
    let mut partitions=BTreeMap::new();
    let mut by_position=BTreeMap::new();
    for s in text.lines().filter_map(|s|s.strip_prefix("STATE\t")) {
        let q=unpack(parse_hex_key(&s.split('\t').collect::<Vec<_>>()).unwrap(),18);
        let mut parts=point_orbits(18,&q[..2]).iter().map(Vec::len).collect::<Vec<_>>();
        parts.sort();
        assert!(parts.len()>1);
        *partitions.entry(parts.clone()).or_insert(0usize)+=1;
        let exceptional_position=q.iter().position(|&p|cycle_type(p,18)==exceptional).unwrap();
        *by_position.entry((exceptional_position,parts)).or_insert(0usize)+=1;
    }
    assert_eq!(partitions.values().sum::<usize>(),40);
    assert_eq!(partitions,BTreeMap::from([(vec![2,16],18),(vec![4,14],14),
        (vec![6,6,6],4),(vec![8,10],4)]));
    assert_eq!(by_position,BTreeMap::from([
        ((0,vec![2,16]),7),((0,vec![6,6,6]),1),((0,vec![8,10]),2),
        ((1,vec![2,16]),7),((1,vec![6,6,6]),1),((1,vec![8,10]),2),
        ((2,vec![2,16]),2),((2,vec![4,14]),7),((2,vec![6,6,6]),1),
        ((3,vec![2,16]),2),((3,vec![4,14]),7),((3,vec![6,6,6]),1),
    ]));
}

#[test]
fn s18_fiber_braid_census_is_reproducible() {
    let root=std::env::temp_dir().join(format!("hurwitz-fiber-census-{}-{}",std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    let args=vec!["hurwitz".to_string(),"s18-fiber-census".to_string(),
        "--certificate".to_string(),"evidence/certificates/CERTIFIED_S18_genus1_weak.txt".to_string(),
        "--out".to_string(),root.display().to_string()];
    s18_fiber_census(&args).unwrap();
    let summary=fs::read_to_string(root.join("SUMMARY.txt")).unwrap();
    assert!(summary.contains("braid-compatible v = 80\ncentralizer orbits = 9\nfull S10 centralizer orbits = 1"));
    assert!(summary.contains("auxiliary three-point genus = 0"));
    let matching=fs::read_to_string(root.join("MATCHING_SLICE.tsv")).unwrap();
    assert_eq!(matching.lines().count(),406);
    let counts=fs::read_to_string(root.join("MATCHING_SUMMARY.tsv")).unwrap();
    assert!(counts.contains("[1, 1, 1, 1, 1, 1, 4, 4, 4]\t1\ttrue\ttrue\t18"));
    let mut totals=[0usize;3];
    for row in counts.lines().skip(1) {
        let fields=row.split('\t').collect::<Vec<_>>();
        let n=fields[4].parse::<usize>().unwrap();
        totals[0]+=n;
        if fields[2]=="true" {totals[1]+=n;}
        if fields[3]=="true" {totals[2]+=n;}
    }
    assert_eq!(totals,[405,390,198]);
    let skeleton=fs::read_to_string(root.join("NINE_BLOCK_SUMMARY.tsv")).unwrap();
    assert!(skeleton.contains("1\t9\t243\t234\t198\t18"));
    assert!(skeleton.contains("3\t6\t162\t156\t0\t0"));
    let mut base_types=BTreeMap::new();
    for row in matching.lines().skip(1) {
        let fields=row.split('\t').collect::<Vec<_>>();
        *base_types.entry((fields[1],fields[2])).or_insert(0usize)+=1;
    }
    assert_eq!(base_types,BTreeMap::from([(("1","[1, 1, 4]"),243),
        (("3","[2, 2, 2]"),81),(("3","[6]"),81)]));
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn balanced_even_paths_and_checkpoint_ranges_are_exact() {
    let mut rng=SplitMix64::new(7131);
    for (cycles,paths,even) in [(vec![],vec![],vec![2]),(vec![],vec![1,1],vec![2]),(vec![2],vec![],vec![2])] {
        let d=cycles.iter().chain(&paths).sum::<usize>()+2*even.iter().sum::<usize>();
        let (a,b)=matching_pair_with_even_paths(&cycles,&paths,&even);
        assert_eq!(cycle_type(a,d),cycle_type(b,d));
        assert_eq!(index(a,d),(d-paths.len()-2*even.len())/2);
        let h0=pair_braid_conjugator(d,a,b,&mut rng).unwrap();
        let target=[b.conjugate_by(a,d),a];
        let mut expected=HashSet::new();let mut h=Perm::identity(d);
        loop {
            if [a,b].map(|p|p.conjugate_by(h,d))==target {expected.insert(h);}
            if !next_permutation(&mut h.img[..d]) {break;}
        }
        let mut actual=HashSet::new();
        enumerate_pair_centralizer(&point_orbits(d,&[a,b]),a,b,0,0,Perm::identity(d),&mut |g| {
            assert!(actual.insert(h0.compose(g,d)));Ok(true)
        }).unwrap();
        assert_eq!(actual,expected);
    }
    let root=std::env::temp_dir().join(format!("hurwitz-checkpoints-{}-{}",std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    let run=|label:&str,h:u64,c:u64,limit:u64| {
        let out=root.join(label);
        let args=["test","hexagon-search","--pair-blocks","2","--even-paths","2","--seconds","30","--out",out.to_str().unwrap(),
            "--start-h",&h.to_string(),"--start-third",&c.to_string(),"--max-thirds",&limit.to_string()].map(str::to_string);
        hexagon_search(&args).unwrap();
        fs::read_to_string(out.join("SUMMARY.txt")).unwrap().lines().filter_map(|l|l.split_once(" = "))
            .map(|(k,v)|(k.to_string(),v.to_string())).collect::<HashMap<_,_>>()
    };
    let whole=run("whole",0,0,u64::MAX);
    assert_eq!(whole["search complete"],"true");
    let (mut h,mut c,mut thirds,mut generating)=(0u64,0u64,0u64,0u64);
    for chunk in 0..1000 {
        let part=run(&format!("chunk{chunk}"),h,c,3);
        assert_eq!(part["scan range resolved"],"true");
        thirds+=part["compatible third involutions visited"].parse::<u64>().unwrap();
        generating+=part["full-group seeds"].parse::<u64>().unwrap();
        let next=(part["scan next h"].parse().unwrap(),part["scan next third"].parse().unwrap());
        if part["traversal exhausted"]=="true" {break;}
        assert!(next>(h,c));(h,c)=next;
        assert!(chunk<999);
    }
    assert_eq!(thirds,whole["compatible third involutions visited"].parse::<u64>().unwrap());
    assert_eq!(generating,whole["full-group seeds"].parse::<u64>().unwrap());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn exact_blocks_and_one_state_braid_relations() {
    let cyclic=cycle_perm(4,&[0,1,2,3]);
    assert_eq!(invariant_blocks(4,&[cyclic]).unwrap(),vec![vec![0,2],vec![1,3]]);
    assert!(invariant_blocks(4,&[cyclic,cycle_perm(4,&[0,1])]).is_none());
    let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_A25_order12_genus12_weak.txt")).unwrap();
    let h:[Perm;3]=std::array::from_fn(|i|tuple_conjugator(case.d,case.seed,hurwitz_move_dynamic(case.d,case.seed,i)).unwrap());
    let d=case.d;
    assert!(h.iter().all(|p|!p.is_even(d)));
    assert_eq!(subgroup_order(d,&h),factorial_u128(d));
    for i in 0..2 {assert_eq!(h[i].compose(h[i+1],d).compose(h[i],d),h[i+1].compose(h[i],d).compose(h[i+1],d));}
    assert_eq!(h[0].compose(h[2],d),h[2].compose(h[0],d));
    let delta=h[0].compose(h[1],d).compose(h[2],d);
    assert!((0..4).fold(Perm::identity(d),|p,_|p.compose(delta,d)).is_identity(d));
}

#[test]
fn odd_path_pair_relabelings_match_small_symmetric_groups() {
    let mut rng=SplitMix64::new(7123);
    for (cycles,paths) in [(vec![],vec![3,3]),(vec![2],vec![3]),(vec![4],vec![1]),(vec![],vec![5])] {
        let d:usize=cycles.iter().chain(&paths).sum();
        let (a,b)=matching_pair_with_paths(&cycles,&paths);
        assert!(a.compose(a,d).is_identity(d) && b.compose(b,d).is_identity(d));
        assert_eq!(index(a,d),(d-paths.len())/2);
        assert_eq!(cycle_type(a,d),cycle_type(b,d));
        let target=[b.conjugate_by(a,d),a];
        let h0=pair_braid_conjugator(d,a,b,&mut rng).unwrap();
        let mut expected=HashSet::new();let mut h=Perm::identity(d);
        loop {
            if [a,b].map(|p|p.conjugate_by(h,d))==target {expected.insert(h);}
            if !next_permutation(&mut h.img[..d]) {break;}
        }
        let mut actual=HashSet::new();
        enumerate_pair_centralizer(&point_orbits(d,&[a,b]),a,b,0,0,Perm::identity(d),&mut |g| {
            assert!(actual.insert(h0.compose(g,d)));Ok(true)
        }).unwrap();
        assert_eq!(actual,expected);
    }
}

#[test]
fn octahedral_lifts_and_s18_position_factor() {
    for m in [1,5] {
        let n=6*m;let (h,k)=octahedral_braid_block(m);
        assert_eq!(h.compose(k,n).compose(h,n),k.compose(h,n).compose(k,n));
        assert_eq!(point_orbits(n,&[h,k]).len(),1);
        let y=h.compose(k,n);let z=y.compose(y,n).compose(y,n);
        for i in 0..n {assert_eq!(z.img[i] as usize,(i/m)*m+(i%m+3)%m);}
    }
    let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt")).unwrap();
    let graph=graph_from_states(case.d,case.states.unwrap()).unwrap();
    let fibers=position_fibers(case.d,&graph).unwrap();
    assert!(fibers.iter().all(|f|f.len()==10));
    let mut v:Vec<_>=(0..40).map(|i|i as f64).collect();
    project_fibers(&mut v,&fibers);
    let av=graph_matvec(&graph,&v);
    for f in &fibers {
        assert!(f.iter().map(|&i|v[i]).sum::<f64>().abs()<1e-10);
        assert!(f.iter().map(|&i|av[i]).sum::<f64>().abs()<1e-10);
    }
    let mut integers=vec![0;40];
    for line in include_str!("../../evidence/spectral_census/s18_position_20260926/NEGATIVE_VECTOR.tsv").lines().skip(1) {
        let columns:Vec<_>=line.split('\t').collect();
        integers[columns[0].parse::<usize>().unwrap()]=columns[2].parse::<i64>().unwrap();
    }
    for f in &fibers {assert_eq!(f.iter().map(|&i|integers[i]).sum::<i64>(),0);}
    assert_eq!(exact_rayleigh(&graph,&integers).unwrap(),(-5569071013532289376,1000000002960119576));
}

#[test]
fn higher_order_rigid_constraints_and_index_bound_match_all_small_permutations() {
    let mut rng=SplitMix64::new(198674);
    for d in 4..=6 {
        for trial in 0..4 {
            let h=if trial==0 {cycle_perm(d,&[0,1])} else {random_perm(d,Ambient::Symmetric,&mut rng)};
            let k=if trial==0 {cycle_perm(d,&[1,2])} else {h};
            let r=h.inverse(d).compose(k.inverse(d),d);
            let mut lengths=[0;MAX_D];
            for c in permutation_cycles(d,k) {for &x in &c {lengths[x]=c.len();}}
            for order in [0,2,3,4,5,6] {
                let mut expected=HashSet::new(); let mut indices=0u64;
                let mut a=Perm::identity(d);
                loop {
                    if a.compose(k,d)==k.compose(a,d)
                        && (0..order).fold(Perm::identity(d),|p,_|p.compose(a,d)).is_identity(d) {
                        indices|=1<<index(a,d);
                        let q=a.compose(r,d);
                        if q.compose(q,d).compose(q,d).is_identity(d) {expected.insert(a);}
                    }
                    if !next_permutation(&mut a.img[..d]) {break;}
                }
                assert_eq!(commuting_entry_indices(d,k,order),indices);
                let mut observed=HashSet::new(); let mut nodes=0;
                assert!(solve_rigid(d,k,r,order,&lengths,[u8::MAX;MAX_D],&mut nodes,100000,
                    Instant::now()+Duration::from_secs(30),&mut |a| {
                        assert!(observed.insert(a));
                        let t=rigid_tuple(d,h,k,a);
                        assert_eq!(t.map(|p|p.conjugate_by(h,d)),hurwitz_move_dynamic(d,t,0));
                        assert_eq!(t.map(|p|p.conjugate_by(k,d)),hurwitz_move_dynamic(d,t,1));
                        Ok(true)
                    }).unwrap());
                assert_eq!(observed,expected,"degree={d} trial={trial} entry_order={order}");
            }
        }
    }
}

#[test]
fn rigid_propagation_matches_independent_small_enumeration() {
    let mut rng=SplitMix64::new(7483);
    for d in 4..=7 {
        for trial in 0..4 {
            let h=if trial==0 {cycle_perm(d,&[0,1])} else {random_perm(d,Ambient::Symmetric,&mut rng)};
            let k=if trial==0 {cycle_perm(d,&[1,2])} else {h};
            let r=h.inverse(d).compose(k.inverse(d),d);
            let mut lengths=[0;MAX_D];
            for c in permutation_cycles(d,k) {for &x in &c {lengths[x]=c.len();}}
            let mut expected=HashSet::new();
            for m in 0..=d/2 {
                let mut class=vec![1;d-2*m]; class.extend(vec![2;m]);
                for a in conjugacy_class(d,&class) {
                    let q=a.compose(r,d);
                    if a.compose(k,d)==k.compose(a,d) && q.compose(q,d).compose(q,d).is_identity(d) {expected.insert(a);}
                }
            }
            let mut observed=HashSet::new(); let mut nodes=0;
            assert!(solve_rigid(d,k,r,2,&lengths,[u8::MAX;MAX_D],&mut nodes,100000,Instant::now()+Duration::from_secs(30),&mut |a| {
                assert!(observed.insert(a));
                let t=rigid_tuple(d,h,k,a);
                assert_eq!(t.map(|p|p.conjugate_by(h,d)),hurwitz_move_dynamic(d,t,0));
                assert_eq!(t.map(|p|p.conjugate_by(k,d)),hurwitz_move_dynamic(d,t,1));
                Ok(true)
            }).unwrap());
            assert_eq!(observed,expected);
            assert!(!solve_rigid(d,k,r,2,&lengths,[u8::MAX;MAX_D],&mut 0,0,Instant::now()+Duration::from_secs(30),&mut |_|panic!("budget exhausted")).unwrap());
        }
    }
    let cert=include_str!("../../evidence/certificates/CERTIFIED_A25_weak.txt");
    let t=parse_tuple(cert.lines().find_map(|s|s.strip_prefix("seed=")).unwrap(),25).unwrap();
    let h=tuple_conjugator(25,t,hurwitz_move_dynamic(25,t,0)).unwrap();
    let k=tuple_conjugator(25,t,hurwitz_move_dynamic(25,t,1)).unwrap();
    let r=h.inverse(25).compose(k.inverse(25),25);
    let q=t[0].compose(r,25);
    assert!(q.compose(q,25).compose(q,25).is_identity(25));
    let mut lengths=[0;MAX_D];
    for c in permutation_cycles(25,k) {for &x in &c {lengths[x]=c.len();}}
    let mut found=false; let mut nodes=0; let mut leaves=0;
    let complete=solve_rigid(25,k,r,2,&lengths,[u8::MAX;MAX_D],&mut nodes,1000000,Instant::now()+Duration::from_secs(30),&mut |a| {
        leaves+=1; found|=a==t[0]; Ok(!found)
    }).unwrap();
    assert!(found && !complete);
    println!("rigid positive control found known A25 a after {nodes} nodes and {leaves} solutions");
    let (d,hh,kk)=glue_braid_blocks(25,h,k,&["9","9","16"].map(str::to_string)).unwrap();
    assert_eq!(d,34);
    assert_eq!(hh.compose(kk,d).compose(hh,d),kk.compose(hh,d).compose(kk,d));
    assert!(glue_braid_blocks(25,h,k,&["16","16","9"].map(str::to_string)).is_err());
    assert!(glue_braid_blocks(25,h,k,&["t20:5".into()]).is_err());
    assert!(glue_braid_blocks(25,h,k,&["t21:0".into()]).is_err());
    let library=braid_block_library(25,h,k);
    assert_eq!(library.keys().copied().collect::<Vec<_>>(),vec![3,4,6,8,9,12,15,16,18,21,24,27,30,33]);
    for (&n,pairs) in &library {
        assert_eq!(pairs.len(),1);
        let (h,k)=pairs[0];
        assert_eq!(h.compose(k,n).compose(h,n),k.compose(h,n).compose(k,n));
        assert_eq!(point_orbits(n,&[h,k]).len(),1);
    }
    for m in 1..=MAX_D/4 {
      for phase in 0..m {
        let n=4*m; let (h,k)=tetrahedral_braid_block(m,phase);
        assert_eq!(h.compose(k,n).compose(h,n),k.compose(h,n).compose(k,n));
        assert_eq!(point_orbits(n,&[h,k]).len(),1);
        let hk=h.compose(k,n); let z=hk.compose(hk,n).compose(hk,n);
        for i in 0..n {assert_eq!(z.img[i] as usize,(i/m)*m+(i%m+3)%m);}
        let key=canonical(&[h,k,Perm::identity(n),product_one_fourth(n,h,k,Perm::identity(n))],n);
        let old=library.get(&n).map(|pairs|pairs.iter().any(|&(a,b)|
            key==canonical(&[a,b,Perm::identity(n),product_one_fourth(n,a,b,Perm::identity(n))],n))).unwrap_or(false);
        assert_eq!(old,m==1 || (m,phase)==(2,0) || (m,phase)==(4,2));
      }
    }
}

#[test]
fn larger_degree_encoding_counts_and_actions_are_exact() {
    let mut rng=SplitMix64::new(83054);
    for d in 26..=MAX_D {
        let a=cycle_perm(d,&(0..d).collect::<Vec<_>>());
        let b=cycle_perm(d,&[0,1]);
        let t=[a,b,Perm::identity(d),product_one_fourth(d,a,b,Perm::identity(d))];
        let key=canonical(&t,d);
        assert!(key.iter().all(|&w|w<factorial_u128(d)));
        assert_eq!(canonical(&unpack(key,d),d),key);
        validate_tuple(d,Ambient::Symmetric,&unpack(key,d)).unwrap();
        assert!(pair_transitive(d,a,b));
        let h=random_perm(d,Ambient::Symmetric,&mut rng);
        assert_eq!(canonical(&t.map(|p|p.conjugate_by(h,d)),d),key);
        for m in 0..=d/2 {
            assert_eq!(commuting_counts(&vec![1;d],m,&mut HashMap::new()),factorial_u128(d)/(factorial_u128(d-2*m)*factorial_u128(m)*(1u128<<m)));
        }
        let c=commuting_involution(d,Perm::identity(d),d/2,&mut rng).unwrap();
        assert_eq!(index(c,d),d/2);
        let invalid=[u128::MAX;3];
        assert_ne!(canonical(&t,d),invalid);
        assert!(unpack(invalid,d)[..3].iter().all(|p|p.img[..d].iter().all(|&x|(x as usize)<d)));
    }
    let quotient=SpectrumGraph{states:vec![],transitions:vec![[0;6]]};
    let h=cycle_perm(34,&[0,33]);
    let edges:Vec<_>=(0..6).map(|i|(0,i,0,h,h)).collect();
    let subsets=subset_action_graph(34,&quotient,&edges,1);
    let ordered=ordered_action_graph(34,&quotient,&edges,1);
    assert_eq!(subsets.transitions[0],[33;6]);
    assert_eq!(ordered.transitions[0],[33;6]);
}

#[test]
fn braid_fixed_construction_matches_brute_force_and_known_a25() {
    let mut rng = SplitMix64::new(20260926);
    for d in 4..=7 {
        for _ in 0..6 {
            let h = random_perm(d,Ambient::Symmetric,&mut rng);
            for m in 0..=d/2 {
                let mut lengths = vec![1;d-2*m]; lengths.extend(vec![2;m]);
                let class = conjugacy_class(d,&lengths);
                let expected = class.iter().filter(|&&c|c.compose(h,d)==h.compose(c,d)).count();
                assert_eq!(commuting_counts(&cycle_type(h,d),m,&mut HashMap::new()),expected as u128);
                let mut enumerated=HashSet::new();
                assert!(enumerate_commuting_involutions(&permutation_cycles(d,h),m,Perm::identity(d),&mut |c| {
                    assert!(enumerated.insert(c)); Ok(true)
                }).unwrap());
                let expected_set:HashSet<_>=class.iter().copied().filter(|c|c.compose(h,d)==h.compose(*c,d)).collect();
                assert_eq!(enumerated,expected_set);
                for _ in 0..4 {
                    let got = commuting_involution(d,h,m,&mut rng);
                    assert_eq!(got.is_some(),expected>0);
                    if let Some(c)=got {
                        assert_eq!(cycle_type(c,d),lengths);
                        assert_eq!(c.compose(h,d),h.compose(c,d));
                    }
                }
                let a = class[rng.below(class.len())];
                let b = class[rng.below(class.len())];
                let found = pair_braid_conjugator(d,a,b,&mut rng);
                if let Some(h)=found {
                    assert_eq!(a.conjugate_by(h,d),b.conjugate_by(a,d));
                    assert_eq!(b.conjugate_by(h,d),a);
                }
                if d<=6 {
                    let mut image:Vec<_>=(0..d as u8).collect();
                    let mut exists=false;
                    loop {
                        let mut h=Perm::identity(d); h.img[..d].copy_from_slice(&image);
                        exists |= a.conjugate_by(h,d)==b.conjugate_by(a,d) && b.conjugate_by(h,d)==a;
                        if !next_permutation(&mut image) {break;}
                    }
                    assert_eq!(found.is_some(),exists);
                }
            }
        }
    }
    // Known weak failure must lie in the constructed family, regardless
    // of its probability under the biased sampler.
    let cert=include_str!("../../evidence/certificates/CERTIFIED_A25_weak.txt");
    let t=parse_tuple(cert.lines().find_map(|s|s.strip_prefix("seed=")).unwrap(),25).unwrap();
    let h=tuple_conjugator(25,t,hurwitz_move_dynamic(25,t,0)).unwrap();
    assert_eq!(t[2].compose(h,25),h.compose(t[2],25));
    assert_eq!(t[3].compose(h,25),h.compose(t[3],25));
    assert_eq!(point_orbit_partition(25,&t[2..]),"15+10");
    let h_cycles=permutation_cycles(25,h);
    let mut cycle_lengths=[0;MAX_D];
    for cycle in &h_cycles { for &x in cycle {cycle_lengths[x]=cycle.len();} }
    for &p in &t[2..] {
        assert!((0..25).all(|x|cycle_lengths[x]==cycle_lengths[p.img[x] as usize]));
    }
    let k=tuple_conjugator(25,t,hurwitz_move_dynamic(25,t,1)).unwrap();
    assert_eq!(h.compose(k,25).compose(h,25),k.compose(h,25).compose(k,25));
    let hk=h.compose(k,25);
    let center=hk.compose(hk,25).compose(hk,25);
    println!("A25 symmetry group order={} orbits={:?} center={} exceptional={} inverse={}",subgroup_order(25,&[h,k]),point_orbits(25,&[h,k]),center.cycle_string(25),center==t[3],center==t[3].inverse(25));
    for block in point_orbits(25,&[h,k]) {
        let mut labels=[u8::MAX;MAX_D];
        for (i,&x) in block.iter().enumerate(){labels[x]=i as u8;}
        let restrict=|p:Perm| {let mut q=Perm::identity(block.len()); for (i,&x) in block.iter().enumerate(){q.img[i]=labels[p.img[x] as usize];} q};
        println!("symmetry block degree={} group_order={} h={:?} k={:?}",block.len(),subgroup_order(block.len(),&[restrict(h),restrict(k)]),cycle_type(restrict(h),block.len()),cycle_type(restrict(k),block.len()));
    }
    let count=commuting_counts(&cycle_type(h,25),10,&mut HashMap::new());
    assert!(count>0);
    println!("A25 <a,b,h> point orbits: {:?}",point_orbits(25,&[t[0],t[1],h]));
    println!("A25 sigma1 relabeling={} type={:?}; commuting m10 involutions={count}; sigma2 relabeling={} type={:?}",h.cycle_string(25),cycle_type(h,25),k.cycle_string(25),cycle_type(k,25));
    assert!((0..10000).any(|_|pair_braid_conjugator(25,t[0],t[1],&mut rng)==Some(h)));
    // Upper bound used by the u64 DP and 64-bit sampling conversion.
    for m in 0..=12 {
        let count=commuting_counts(&vec![1;25],m,&mut HashMap::new());
        assert_eq!(count as u128,factorial_u128(25)/(factorial_u128(25-2*m)*factorial_u128(m)*(1u128<<m)));
        assert!(count<usize::MAX as u128);
    }
    let mut visits=0;
    assert!(!enumerate_commuting_involutions(&permutation_cycles(25,h),10,Perm::identity(25),&mut |_| {
        visits+=1; Ok(false)
    }).unwrap());
    assert_eq!(visits,1);
    let counts=std::array::from_fn(|_|AtomicU64::new(0));
    let seed=(0..10).find_map(|_|braid_fixed_seed(5,Ambient::Alternating,2,0,None,&mut rng,&counts)).unwrap();
    validate_tuple(5,Ambient::Alternating,&seed).unwrap();
    assert_eq!(genus_of(&seed,5),0);
    assert_eq!(canonical(&seed,5),canonical(&hurwitz_move_dynamic(5,seed,0),5));
    let profile=cycle_type(seed[3],5);
    let seed=(0..100).find_map(|_|braid_fixed_seed(5,Ambient::Alternating,2,0,Some(&profile),&mut rng,&counts)).unwrap();
    assert_eq!(cycle_type(seed[3],5),profile);
}

fn brute_order(d: usize, gens: &[Perm]) -> usize {
    let mut seen = HashSet::from([Perm::identity(d)]);
    let mut q = vec![Perm::identity(d)];
    let mut head = 0;
    while head < q.len() {
        let x = q[head]; head += 1;
        for g in gens { let y = g.compose(x, d); if seen.insert(y) { q.push(y); } }
    }
    seen.len()
}

#[test]
fn voltage_lift_and_loop_conditioning_agree_on_a5() {
    let d = 5;
    let seed = parse_tuple("[(), (1 2 3 4 5), (1 2 4 5 3), (1 4 5 2 3)]", d).unwrap();
    let graph = enumerate_spectrum_graph(d, seed).unwrap();
    assert_eq!(graph.transitions.len(), 4);
    let profile = lift_profile(d, &graph).unwrap();
    let lift = labeled_lift_graph(d, unpack(graph.states[0], d), 1000).unwrap();
    assert_eq!(lift.transitions.len() as u128, profile.holonomy_order * 4);
    let standard = standard_point_lambda_max(d, 4, &profile.edges).unwrap();
    let full = eigen_for_graph(&lift, 1000, 128, 819);
    assert!(standard <= full.lambda2 + 1e-8);
    let point = subset_action_graph(d, &graph, &profile.edges, 1);
    let point_lambda2 = eigen_for_graph(&point, 1000, 128, 819).lambda2;
    assert!((point_lambda2 - standard.max(eigen_for_graph(&graph, 1000, 128, 819).lambda2)).abs() < 1e-8);
    let pairs = subset_action_graph(d, &graph, &profile.edges, 2);
    assert!(eigen_for_graph(&pairs, 1000, 128, 819).lambda2 <= full.lambda2 + 1e-8);
    let ordered_triples = ordered_action_graph(d, &graph, &profile.edges, 3);
    assert_eq!(ordered_triples.transitions.len(), lift.transitions.len());
    assert!((eigen_for_graph(&ordered_triples, 1000, 128, 819).lambda2 - full.lambda2).abs() < 1e-8);
    let mut rng = SplitMix64::new(1873);
    for _ in 0..16 {
        let null = loop_matched_schreier_graph(&graph, &mut rng);
        for (i, row) in graph.transitions.iter().enumerate() {
            for mv in 0..6 { assert_eq!(row[mv] as usize == i, null.transitions[i][mv] as usize == i); }
        }
    }
}

#[test]
fn exact_orders_match_independent_enumeration_and_old_engine() {
    let mut rng = SplitMix64::new(88123);
    for d in 4..=10 {
        for _ in 0..30 {
            let gens: Vec<_> = (0..3).map(|_| random_perm(d, Ambient::Symmetric, &mut rng)).collect();
            for k in 1..=3 {
                let got = subgroup_order(d, &gens[..k]);
                assert_eq!(got, subgroup_order_reference(d, &gens[..k]));
                if d <= 7 { assert_eq!(got, brute_order(d, &gens[..k]) as u128); }
            }
        }
    }
    for d in 4..=MAX_D {
        let cycle: Vec<_> = (0..d).collect();
        assert_eq!(subgroup_order(d, &[cycle_perm(d,&cycle),cycle_perm(d,&[0,1])]), factorial_u128(d));
        let alt_cycle: Vec<_> = (0..if d%2 == 1 {d} else {d-1}).collect();
        assert_eq!(subgroup_order(d, &[cycle_perm(d,&alt_cycle),cycle_perm(d,&[0,1,d-1])]), factorial_u128(d)/2);
    }
}

#[test]
fn canonicalization_and_known_counterexample() {
    let seed = parse_tuple(A8, 8).unwrap();
    validate_tuple(8, Ambient::Alternating, &seed).unwrap();
    let key = canonical(&seed, 8);
    let mut rng = SplitMix64::new(73452);
    for _ in 0..100 {
        let h = random_perm(8, Ambient::Symmetric, &mut rng);
        assert_eq!(key, canonical(&seed.map(|g| g.conjugate_by(h,8)),8));
        assert_eq!(key, canonical(&unpack(key,8),8));
        for mv in 0..6 {
            let moved = hurwitz_move_dynamic(8, seed, mv);
            let conjugated = seed.map(|g| g.conjugate_by(h,8));
            assert_eq!(canonical(&moved,8),canonical(&hurwitz_move_dynamic(8,conjugated,mv),8));
        }
    }
    let result = explore(8,seed,10000,Instant::now()+Duration::from_secs(30));
    assert_eq!(result.status, "certified_strong_failure");
    // Quotienting by S8 instead of A8 halves this example again.
    assert_eq!(result.states.len(),144);
    assert!(result.weak_witness.is_some());
    let expected = BTreeMap::from([(12,12),(18,18),(24,16),(36,12),(168,44),(360,8),(2520,34)]);
    assert_eq!(result.orders, expected);
    assert_eq!(explore(8,seed,1,Instant::now()+Duration::from_secs(30)).status,"capped");
    assert_eq!(explore(8,seed,1000,Instant::now()).status,"timeout");
    assert!(parse_tuple("[(0 1), (), (), ()]",8).is_err());
    assert!(parse_tuple("[(1 2)(2 3), (), (), ()]",8).is_err());
}

#[test]
fn a9_saved_candidate_independent_check() {
    for (text,size) in [
        ("[(1 9)(2 8)(3 7)(4 6), (1 6)(2 4)(3 5)(7 8), (1 9)(2 6)(3 8)(4 7), (2 4 9)(5 8 7)]",144),
        ("[(1 2)(3 4)(5 8)(6 7), (1 4 3)(6 7 8), (2 6 3)(5 9 7), (1 3)(2 6)(5 9)(7 8)]",216),
        ("[(1 2)(3 6)(5 8)(7 9), (3 6 5)(4 7 8), (1 5 7)(2 4 8), (1 8 6)(2 7 9)]",972),
        ("[(1 2 3)(6 9 7), (2 5 3)(6 7 9), (2 8 5)(3 4 6), (1 8 2)(3 6 4)]",360),
        ("[(1 6 3)(5 8 7), (2 7 3)(4 9 6), (1 7 6)(5 9 8), (2 3 5)(4 7 9)]",36),
        ("[(1 7 5)(2 9 6)(3 8 4), (1 2 3)(4 6 5)(7 8 9), (1 8 4)(2 3 7)(5 6 9), (1 5 6)(2 8 7)(3 9 4)]",36),
    ] {
    let seed=parse_tuple(text,9).unwrap();
    assert_eq!(brute_order(9,&seed[..3]),181440);
    let result=explore(9,seed,10000,Instant::now()+Duration::from_secs(30));
    assert_eq!(result.status,"certified_strong_failure");
    assert_eq!(result.states.len(),size);
    assert!(result.weak_witness.is_some());
    let closed:HashSet<_>=result.states.iter().copied().collect();
    for &key in &result.states {
        let t=unpack(key,9);
        let brute=brute_order(9,&t[..2]) as u128;
        assert_eq!(subgroup_order(9,&t[..2]),brute);
        assert!(brute<181440);
        for mv in 0..6 { assert!(closed.contains(&canonical(&hurwitz_move_dynamic(9,t,mv),9))); }
    }
    }
}

#[test]
fn certificate_roundtrip_rejects_corrupted_closure() {
    let d=8;
    let seed=parse_tuple(A8,d).unwrap();
    let result=explore(d,seed,10000,Instant::now()+Duration::from_secs(30));
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let path=std::env::temp_dir().join(format!("hurwitz-test-{}-{stamp}.txt",std::process::id()));
    write_certificate(&path,d,Ambient::Alternating,&seed,&result).unwrap();
    verify(&path).unwrap();
    let data=fs::read_to_string(&path).unwrap();
    let last=data.lines().last().unwrap();
    fs::write(&path,data.replace(&format!("{last}\n"),"")).unwrap();
    assert!(verify(&path).is_err());
    fs::remove_file(path).unwrap();
}

#[test]
fn spectral_graph_and_jacobi_sanity() {
    let labeled = SpectrumGraph { states: Vec::new(), transitions: vec![
        [1, 2, 1, 2, 1, 2], [2, 0, 2, 0, 2, 0], [0, 1, 0, 1, 0, 1],
    ] };
    let permutation = [2, 0, 1];
    let mut relabeled = labeled.clone();
    for (old, row) in labeled.transitions.iter().enumerate() {
        relabeled.transitions[permutation[old]] = row.map(|next| permutation[next as usize] as u32);
    }
    assert_eq!(canonical_labeled_code(&labeled, 3), canonical_labeled_code(&relabeled, 3));
    assert_eq!(exact_rayleigh(&labeled, &[-1, 0, 1]).unwrap(), (-6, 2));
    assert!(exact_rayleigh(&labeled, &[1, 0, 1]).is_err());
    assert_eq!(canonical_labeled_code(&labeled, 2), "");
    assert_eq!(infer_group_degree(Path::new("/tmp/A7_shards/shard_0/generating_components.csv")).unwrap(),(Ambient::Alternating,7));
    let a=Perm::identity(5);
    let b=perm_of_type(&[3,1,1],5,None);
    let c=perm_of_type(&[5],5,None);
    let seed=[a,b,c,product_one_fourth(5,a,b,c)];
    validate_tuple(5,Ambient::Alternating,&seed).unwrap();
    let graph=enumerate_spectrum_graph(5,seed).unwrap();
    assert!(!graph.states.is_empty());
    assert_eq!(graph.states.len(),graph.transitions.len());
    assert_eq!(graph_components(&graph),1);
    assert!(graph.transitions.iter().all(|row| row.len()==6));
    let values=jacobi_eigenvalues(vec![2.0,1.0,1.0,2.0],2);
    assert!((values[0]-3.0).abs()<1e-9 && (values[1]-1.0).abs()<1e-9);
    let exact=eigen_for_graph(&graph,usize::MAX,graph.transitions.len(),7);
    let (lambda2,lambda_min)=lanczos_extremes(&graph,graph.transitions.len(),7);
    assert!((lambda2-exact.lambda2).abs()<1e-8);
    assert!((lambda_min-exact.lambda_min).abs()<1e-8);
    let mut rng=SplitMix64::new(42);
    let diagnostic=graph_diagnostics(&graph,&exact,8,4,3,0,11);
    assert_eq!((diagnostic.braid_12_defects,diagnostic.braid_23_defects,diagnostic.commute_13_defects),(0,0,0));
    assert_eq!(diagnostic.moments.len(),8);
    assert!(diagnostic.maxcut_conflict_fraction.is_finite());
    let matched=cycle_matched_schreier_graph(&graph,&mut rng);
    let matched_diagnostic=graph_diagnostics(&matched,&eigen_for_graph(&matched,usize::MAX,matched.transitions.len(),8),0,1,0,0,12);
    assert_eq!(diagnostic.generator_profiles,matched_diagnostic.generator_profiles);
    for n in [7,8,31] {
        let regular=random_regular_graph(n,6,20,&mut rng).unwrap();
        for (i,row) in regular.transitions.iter().enumerate() {
            assert_eq!(row.iter().copied().collect::<HashSet<_>>().len(),6);
            assert!(!row.contains(&(i as u32)));
            for &j in row { assert!(regular.transitions[j as usize].contains(&(i as u32))); }
        }
    }
    assert_eq!(finite_stats([1.0, f64::NAN, 3.0]), (2, 2.0, 1.0, 3.0));
    assert!((spectral_wasserstein(&[6.0, 2.0, -2.0], &[6.0, 1.0, -1.0])-1.0).abs()<1e-12);
}

#[test]
fn a25_weak_certificate_seed_is_supported() {
    const A25: &str = "[(1 10)(3 5)(4 24)(6 13)(8 21)(9 19)(11 16)(17 22)(18 20)(23 25), (2 6)(3 7)(8 16)(9 17)(10 12)(11 18)(14 15)(19 25)(20 22)(21 23), (1 11)(2 23)(3 17)(4 22)(5 13)(6 10)(7 12)(9 15)(14 21)(16 24), (1 7 17 20)(2 10 11 8)(3 13 23 19)(4 16 18)(5 12 6)(9 25 14)(15 21 24 22)]";
    let seed = parse_tuple(A25, 25).unwrap();
    validate_tuple(25, Ambient::Alternating, &seed).unwrap();
    assert_eq!(point_orbit_count_perms(25, seed[0], seed[1]), 6);
    assert_eq!(subgroup_order(25, &seed[..3]), factorial_u128(25) / 2);
    for i in 0..4 {
        for j in i+1..4 {
            assert_eq!(cycle_type(seed[i].compose(seed[j], 25), 25), [2, 2, 3, 3, 3, 6, 6]);
        }
    }
    let root = canonical(&seed, 25);
    let relabeled = unpack(root, 25);
    validate_tuple(25, Ambient::Alternating, &relabeled).unwrap();
    assert_eq!(canonical(&relabeled, 25), root);
}

#[test]
fn degree_eleven_involution_class_has_expected_size() {
    let lengths = [1, 1, 1, 2, 2, 2, 2];
    let class = conjugacy_class(11, &lengths);
    // 11! / (2^4 * 4! * 3!), the centralizer class-size formula.
    assert_eq!(class.len(), 17325);
    assert_eq!(class.iter().copied().collect::<HashSet<_>>().len(), class.len());
    assert!(class.iter().all(|&p| cycle_type(p, 11) == lengths));
}

#[test]
fn weak_search_and_low_index_samples() {
    let seed=parse_tuple(A8,8).unwrap();
    let result=explore_goal(8,seed,10000,Instant::now()+Duration::from_secs(30),true);
    assert_eq!(result.status,"weak_found");
    let mut rng=SplitMix64::new(82112);
    for d in [8,9,13,16] {
        for ambient in [Ambient::Alternating,Ambient::Symmetric] {
            for _ in 0..5 {
                let seed=low_index_seed(d,ambient,false,&mut rng).expect("low index sample");
                validate_tuple(d,ambient,&seed).unwrap();
                let total:usize=seed.iter().map(|&p|index(p,d)).sum();
                assert!((2*d-2..=2*d+2).contains(&total));
            }
        }
    }
}

#[test]
fn targeted_families_and_probe_accounting() {
    let mut rng=SplitMix64::new(91357);
    for d in 4..=MAX_D {
        for _ in 0..100 {
            let a=random_perm(d,Ambient::Symmetric,&mut rng);
            let b=random_perm(d,Ambient::Symmetric,&mut rng);
            assert_eq!(pair_transitive(d,a,b),point_orbit_count_perms(d,a,b)==1);
        }
        assert!(!pair_transitive(d,Perm::identity(d),Perm::identity(d)));
    }
    // These fixed braid words move each chosen pair to slots 0,1 up to
    // simultaneous conjugacy, proving the all-six-pair rejection rule.
    let words=[((0,1),&[][..]),((0,2),&[4][..]),((0,3),&[5,4][..]),
               ((1,2),&[0,1][..]),((1,3),&[5,0,1][..]),((2,3),&[1,2,0,1][..])];
    for d in 4..=10 {
        for _ in 0..20 {
            let a=random_perm(d,Ambient::Symmetric,&mut rng);
            let b=random_perm(d,Ambient::Symmetric,&mut rng);
            let c=random_perm(d,Ambient::Symmetric,&mut rng);
            let t=[a,b,c,product_one_fourth(d,a,b,c)];
            for &((i,j),word) in &words {
                let front=word.iter().fold(t,|u,&mv|hurwitz_move_dynamic(d,u,mv));
                assert_eq!(pair_transitive(d,t[i],t[j]),pair_transitive(d,front[0],front[1]));
            }
        }
    }
    let passport=parse_passport("3.3.3/3.3.3/3.3.3/3.3.3",9,Ambient::Alternating).unwrap();
    let small=parse_passport("3.3/3.3/3.3/3.3",9,Ambient::Alternating).unwrap();
    assert_eq!(conjugacy_class(9,&small[1]).len(),3360);
    for _ in 0..5 {
        let t=passport_seed(9,Ambient::Alternating,&passport,false,&mut rng).unwrap();
        validate_tuple(9,Ambient::Alternating,&t).unwrap();
        for i in 0..4 { assert_eq!(cycle_type(t[i],9),passport[i]); }
    }
    assert!(parse_passport("3/3/3/3",9,Ambient::Alternating).is_err());
    assert!(parse_passport("3.3.3/3.3.3/3.3.3/3.3.3",9,Ambient::Symmetric).is_err());
    assert!(parse_passport("0/1/1/1",9,Ambient::Alternating).is_err());
    let parents=vec![(8,parse_tuple(A8,8).unwrap())];
    for d in 8..=10 {
        for ambient in [Ambient::Alternating,Ambient::Symmetric] {
            for _ in 0..3 {
                let t=deform_seed(d,ambient,&parents,&mut rng).unwrap();
                validate_tuple(d,ambient,&t).unwrap();
            }
        }
    }
    let t=parse_tuple(A8,8).unwrap();
    let (success,_,step)=probe(8,t,&mut rng,64,false);
    assert!(!success); assert_eq!(step,65);
    let (success,weak,step)=probe(8,t,&mut rng,64,true);
    assert!(success && weak && step<=64);
}

#[test]
fn local_switch_preserves_entry_types_and_filters_candidates() {
    let case=parse_certificate_case(Path::new("evidence/certificates/CERTIFIED_S18_genus1_weak.txt")).unwrap();
    let states=case.states.unwrap().into_iter().map(|key|unpack(key,18)).collect::<Vec<_>>();
    let parent_types=states.iter().map(|s|s[..3].iter().map(|&p|cycle_type(p,18)).collect::<Vec<_>>()).collect::<HashSet<_>>();
    let mut rng=SplitMix64::new(20260929);
    let counts: [AtomicU64;6]=std::array::from_fn(|_|AtomicU64::new(0));
    for _ in 0..200 {
        if let Some(t)=local_switch_seed(18,Ambient::Symmetric,&states,2,18,false,false,&mut rng,&counts) {
            validate_tuple(18,Ambient::Symmetric,&t).unwrap();
            assert!(!pair_transitive(18,t[0],t[1]));
            assert!(parent_types.contains(&t[..3].iter().map(|&p|cycle_type(p,18)).collect::<Vec<_>>()));
        }
    }
    assert_eq!(counts[0].load(Ordering::Relaxed),200);
    assert!(counts[1].load(Ordering::Relaxed)>0);
    for i in 1..6 {assert!(counts[i].load(Ordering::Relaxed)<=counts[i-1].load(Ordering::Relaxed));}
}

#[test]
fn orbit_map_distinguishes_closed_failure_from_partial_search() {
    let certificate = include_str!("../../evidence/certificates/CERTIFIED_A25_weak.txt");
    let seed_text = certificate.lines().find_map(|line| line.strip_prefix("seed=")).unwrap();
    let seed = parse_tuple(seed_text, 25).unwrap();
    let empty = HashMap::new();
    let closed = map_orbit(25, seed, 16, Instant::now() + Duration::from_secs(30), &empty);
    assert_eq!(closed.status, "closed_weak_failure");
    assert_eq!(closed.states.len(), 4);
    assert!(!closed.weak);
    let capped = map_orbit(25, seed, 1, Instant::now() + Duration::from_secs(30), &empty);
    assert_eq!(capped.status, "capped");
    let known: HashMap<_, _> = closed.states.iter().copied().map(|key| (key, 1)).collect();
    let neighbor = hurwitz_move_dynamic(25, seed, 1);
    assert!(known.contains_key(&canonical(&neighbor, 25)));
}

#[test]
fn commuting_powers_match_independent_permutations_and_resume() {
    for d in 4..=7 {
        let mut p=Perm::identity(d);let mut all=Vec::new();let mut types=BTreeMap::new();
        loop {
            all.push(p);types.entry(cycle_type(p,d)).or_insert(p);
            if !next_permutation(&mut p.img[..d]) {break;}
        }
        for (lengths,h) in types {for q in 2..=6 {
            let expected:HashSet<_>=all.iter().copied().filter(|&c|c.compose(h,d)==h.compose(c,d) && power_identity(d,c,q)).collect();
            let cycles=permutation_cycles(d,h);let mut got=Vec::new();
            let deadline=Instant::now()+Duration::from_secs(30);
            assert!(enumerate_commuting_power(&cycles,0,q,Perm::identity(d),deadline,&mut 0,&mut HashMap::new(),&mut |c| {got.push(c);Ok(true)}).unwrap());
            assert_eq!(got.len(),expected.len());assert_eq!(got.iter().copied().collect::<HashSet<_>>(),expected);
            assert_eq!(commuting_power_count(&lengths,q,&mut HashMap::new()),got.len() as u128);
            for skip in [0,got.len()/2,got.len()] {
                let mut left=skip as u128;let mut tail=Vec::new();
                assert!(enumerate_commuting_power(&cycles,0,q,Perm::identity(d),deadline,&mut left,&mut HashMap::new(),&mut |c| {tail.push(c);Ok(true)}).unwrap());
                assert_eq!(left,0);assert_eq!(tail,got[skip..]);
            }
        }}
    }
}

#[test]
fn new_modular_representations_and_rigid_branch_resume_are_exact() {
    for (n,m) in [(8,2),(9,3),(12,2),(15,5),(18,3),(20,5),(25,5),(30,5),(32,2)] {
        for seed in 0..4 {
            let (h,k)=modular_braid_block(n,m,seed).unwrap();
            assert_eq!(h.compose(k,n).compose(h,n),k.compose(h,n).compose(k,n));
            assert!(generators_transitive(n,&[h,k]));
            let hk=h.compose(k,n);let center=hk.compose(hk,n).compose(hk,n);
            let mut expected=Perm::identity(n);
            for x in 0..n {expected.img[x]=(x/m*m+(x%m+1)%m) as u8;}
            assert_eq!(center,expected);
        }
    }
    let d=4;let (h,k)=tetrahedral_braid_block(1,0);
    let mut lengths=[0;MAX_D];
    for c in permutation_cycles(d,k) {for &x in &c {lengths[x]=c.len();}}
    let r=h.inverse(d).compose(k.inverse(d),d);let deadline=Instant::now()+Duration::from_secs(30);
    let mut whole=Vec::new();let mut nodes=0;
    assert!(solve_rigid(d,k,r,0,&lengths,[u8::MAX;MAX_D],&mut nodes,u64::MAX,deadline,&mut |a| {whole.push(a);Ok(true)}).unwrap());
    let mut chunked=Vec::new();let mut cursor=0;
    loop {
        let (complete,next,_)=rigid_branches(d,k,r,0,&lengths,cursor,&mut 0,u64::MAX,deadline,&mut |a| {
            chunked.push(a);Ok(true)
        }).unwrap();
        cursor=next;if complete {break;}
    }
    assert_eq!(whole,chunked);
    assert_eq!(rigid_branches(d,k,r,0,&lengths,cursor,&mut 0,u64::MAX,deadline,&mut |_|panic!("finished")).unwrap().1,cursor);
}

#[test]
fn quotient_relabelings_preserve_every_transitive_completion() {
    let d=6;let (a,b)=matching_pair(&[2,4]);
    let h0=pair_braid_conjugator(d,a,b,&mut SplitMix64::new(1)).unwrap();
    for q in 2..=6 {
        let mut all=HashSet::new();let mut reduced=HashSet::new();let mut h_seen=HashSet::new();
        enumerate_pair_centralizer(&point_orbits(d,&[a,b]),a,b,0,0,Perm::identity(d),&mut |g| {
            let h=h0.compose(g,d);let keep=h_seen.insert(colored_representation_key(d,[a,b,h]));
            enumerate_commuting_power(&permutation_cycles(d,h),0,q,Perm::identity(d),Instant::now()+Duration::from_secs(30),&mut 0,&mut HashMap::new(),&mut |c| {
                let t=[a,b,c,product_one_fourth(d,a,b,c)];
                if generators_transitive(d,&t[..3]) {let key=canonical(&t,d);all.insert(key);if keep {reduced.insert(key);}}
                Ok(true)
            })?;
            Ok(true)
        }).unwrap();
        if q<=4 {assert!(!all.is_empty());}assert_eq!(all,reduced);
    }
}

#[test]
fn higher_order_search_checkpoint_ranges_are_exact() {
    let root=std::env::temp_dir().join(format!("hurwitz-power-cursors-{}-{}",std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    for q in [3,4,5,6] {
        let run=|label:&str,h:u64,c:u64,limit:u64| {
            let out=root.join(format!("q{q}_{label}"));
            let args=["test","hexagon-search","--pair-blocks","2","--even-paths","2","--third-order",&q.to_string(),
                "--seconds","30","--out",out.to_str().unwrap(),"--start-h",&h.to_string(),
                "--start-third",&c.to_string(),"--max-thirds",&limit.to_string()].map(str::to_string);
            hexagon_search(&args).unwrap();
            fs::read_to_string(out.join("SUMMARY.txt")).unwrap().lines().filter_map(|l|l.split_once(" = "))
                .map(|(k,v)|(k.to_string(),v.to_string())).collect::<HashMap<_,_>>()
        };
        let whole=run("whole",0,0,u64::MAX);
        assert_eq!(whole["search complete"],"true");
        assert_eq!(whole["group"],if q%2==1 {"A_6"} else {"auto_6"});
        assert!(!whole.contains_key("third involution transpositions"));
        let (mut h,mut c,mut thirds,mut generating)=(0u64,0u64,0u64,0u64);
        let mut done=false;
        for chunk in 0..1000 {
            let part=run(&format!("chunk{chunk}"),h,c,3);
            assert_eq!(part["scan range resolved"],"true");
            thirds+=part["compatible commuting third entries visited"].parse::<u64>().unwrap();
            generating+=part["full-group seeds"].parse::<u64>().unwrap();
            let next=(part["scan next h"].parse().unwrap(),part["scan next third"].parse().unwrap());
            if part["traversal exhausted"]=="true" {done=true;break;}
            assert!(next>(h,c));(h,c)=next;
        }
        assert!(done);
        assert_eq!(thirds,whole["compatible commuting third entries visited"].parse::<u64>().unwrap());
        assert_eq!(generating,whole["full-group seeds"].parse::<u64>().unwrap());
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn unequal_third_involutions_preserve_braid_symmetry_and_exhaust_class() {
    let d=8;
    let (a,b)=matching_pair_with_even_paths(&[2,2],&[1,1,1,1],&[]);
    let mut rng=SplitMix64::new(20260929);
    let h0=pair_braid_conjugator(d,a,b,&mut rng).unwrap();
    let class=conjugacy_class(d,&[2,2,2,2]);
    let mut count=0;
    enumerate_pair_centralizer(&point_orbits(d,&[a,b]),a,b,0,0,Perm::identity(d),&mut |g| {
        let h=h0.compose(g,d);
        let expected:HashSet<_>=class.iter().copied().filter(|&c|c.compose(h,d)==h.compose(c,d)).collect();
        let mut actual=HashSet::new();
        enumerate_commuting_involutions(&permutation_cycles(d,h),4,Perm::identity(d),&mut |c| {
            let t=[a,b,c,product_one_fourth(d,a,b,c)];
            assert_eq!(t.map(|x|x.conjugate_by(h,d)),hurwitz_move_dynamic(d,t,0));
            assert!(actual.insert(c));count+=1;Ok(true)
        })?;
        assert_eq!(actual,expected);Ok(true)
    }).unwrap();
    assert!(count>0);
    let mut a=Perm::identity(4);
    loop {
        let mut b=Perm::identity(4);
        loop {
            let mut c=Perm::identity(4);
            loop {
                let t=[a,b,c,product_one_fourth(4,a,b,c)];
                if point_orbits(4,&t).len()==1 && t.iter().any(|&p|cycle_type(p,4).len()<=2) {
                    assert!(has_transitive_pair(4,&t));
                }
                if !next_permutation(&mut c.img[..4]) {break;}
            }
            if !next_permutation(&mut b.img[..4]) {break;}
        }
        if !next_permutation(&mut a.img[..4]) {break;}
    }
}

#[test]
fn direct_classes_and_centralizer_rows_match_independent_enumeration() {
    for d in 4..=6 {
        let mut p=Perm::identity(d);
        let mut classes=HashMap::<Vec<usize>,Vec<Perm>>::new();
        loop {
            classes.entry(cycle_type(p,d)).or_default().push(p);
            if !next_permutation(&mut p.img[..d]) {break;}
        }
        let all:Vec<_>=classes.values().flatten().copied().collect();
        for (lengths,rows) in &classes {
            assert_eq!(&conjugacy_class(d,lengths),rows);
            let signature=lengths.iter().fold(0u64,|s,&l|s+(1u64<<(4*l)));
            assert!(rows.iter().all(|&p|small_cycle_signature(d,p)==signature));
            let a=perm_of_type(lengths,d,None);
            let centralizer:Vec<_>=all.iter().copied().filter(|&g|g.compose(a,d)==a.compose(g,d)).collect();
            for bs in classes.values() {
                let representatives=centralizer_rows(d,a,bs.clone());
                let mut covered=HashSet::new();
                for b in representatives {
                    let orbit:HashSet<_>=centralizer.iter().map(|&g|b.conjugate_by(g,d)).collect();
                    assert!(orbit.is_disjoint(&covered));
                    covered.extend(orbit);
                }
                assert_eq!(covered,bs.iter().copied().collect());
            }
        }
    }
}

#[test]
fn exact_passport_chunks_match_brute_a4() {
    let passport=parse_passport("3/3/3/3",4,Ambient::Alternating).unwrap();
    let class=conjugacy_class(4,&passport[1]);
    assert_eq!(class.len(),8);
    assert_eq!(class.iter().copied().collect::<HashSet<_>>().len(),8);
    let a=perm_of_type(&passport[0],4,None);
    let accepted:Vec<_>=(0..64).map(|i| {
        let b=class[i/8]; let c=class[i%8];
        cycle_type(product_one_fourth(4,a,b,c),4)==passport[3] && brute_order(4,&[a,b,c])==12
    }).collect();
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let root=std::env::temp_dir().join(format!("hurwitz-scan-test-{}-{stamp}",std::process::id()));
    for (start,count) in [(0,64),(0,13),(13,51)] {
        let out=root.join(format!("{start}_{count}"));
        let args:Vec<String>=["hurwitz","passport-exact","--group","A","--d","4","--passport","3/3/3/3","--threads","3","--start",&start.to_string(),"--count",&count.to_string(),"--out",out.to_str().unwrap()].iter().map(|s|s.to_string()).collect();
        run(&args).unwrap();
        let summary=fs::read_to_string(out.join("SUMMARY.txt")).unwrap();
        let expected=accepted[start..start+count].iter().filter(|&&b|b).count();
        assert!(summary.contains(&format!("accepted generating seeds = {expected}\n")));
        assert!(summary.contains(&format!("scan next = {}\n",start+count)));
        assert!(summary.contains("scan range resolved = true\n"));
    }
    let candidate_b:Vec<_>=class.iter().copied().filter(|&b|!pair_transitive(4,a,b)).collect();
    let out=root.join("weak");
    let args:Vec<String>=["hurwitz","passport-exact","--group","A","--d","4","--passport","3/3/3/3","--threads","3","--start","0","--count","64","--weak-only","--out",out.to_str().unwrap()].iter().map(|s|s.to_string()).collect();
    run(&args).unwrap();
    let summary=fs::read_to_string(out.join("SUMMARY.txt")).unwrap();
    let candidate_total=candidate_b.len()*candidate_b.len();
    assert!(summary.contains("scan raw class pairs = 64\n"));
    assert!(summary.contains(&format!("scan next = {candidate_total}\n")));
    assert!(summary.contains(&format!("fixed-entry transitive pairs pruned = {}\n",64-candidate_total)));
    assert!(summary.contains("scan range resolved = true\n"));
    let out=root.join("centralizer");
    let args:Vec<String>=["hurwitz","passport-exact","--group","A","--d","4","--passport","3/3/3/3","--threads","3","--centralizer-rows","--weak-only","--out",out.to_str().unwrap()].iter().map(|s|s.to_string()).collect();
    run(&args).unwrap();
    let summary=fs::read_to_string(out.join("SUMMARY.txt")).unwrap();
    assert!(summary.contains("scan raw class pairs = 64\n"));
    assert!(summary.contains("scan enumeration = centralizer-v1\n"));
    assert!(summary.contains("scan range resolved = true\n"));
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn replay_conjugacy_duplicates_are_not_unresolved() {
    let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let root=std::env::temp_dir().join(format!("hurwitz-promoted-repeats-{}-{stamp}",std::process::id()));
    fs::create_dir_all(&root).unwrap();let input=root.join("replay.txt");
    let tuple="[(1 2 3 4), (1 2), (1 2), (1 4 3 2)]";
    let t=parse_tuple(tuple,4).unwrap();assert!(tuple_product_is_identity(4,&t));
    assert!(tuple_generates_dynamic(4,Ambient::Symmetric,&t));
    let mut h=Perm::identity(4);h.img.swap(0,2);
    let other=tuple_string_dynamic(4,&t.map(|p|p.conjugate_by(h,4)));
    fs::write(&input,format!("group=S_4 seed={tuple}\ngroup=S_4 seed={other}\n")).unwrap();
    let args=["hurwitz".into(),"smart-hunt".into(),"--group".into(),"S".into(),"--d".into(),"4".into(),
        "--threads".into(),"1".into(),"--weak-only".into(),"--replay-only".into(),"--replay".into(),
        input.to_string_lossy().into_owned(),"--out".into(),root.join("out").to_string_lossy().into_owned()];
    run(&args).unwrap();let summary=fs::read_to_string(root.join("out/SUMMARY.txt")).unwrap();
    for expected in ["replay inputs = 1\n","accepted generating seeds = 1\n","solved seeds = 1\n","duplicate promoted starts = 0\n","hard unsolved candidates = 0\n"] {
        assert!(summary.contains(expected),"{summary}");
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
#[ignore = "explicit reproducible microbenchmark"]
fn benchmark_orders() {
    for d in [9,13,16] {
        let mut rng = SplitMix64::new(76123);
        let inputs: Vec<Vec<Perm>> = (0..16).map(|_| (0..3).map(|_| random_perm(d,Ambient::Alternating,&mut rng)).collect()).collect();
        for (name, f) in [("reference",subgroup_order_reference as fn(usize,&[Perm])->u128),("filtered",subgroup_order)] {
            let start=Instant::now();
            let mut sum=0;
            for gens in &inputs { sum+=f(d,gens); }
            println!("d={d} engine={name} seconds={:.6} checksum={sum}",start.elapsed().as_secs_f64());
        }
    }
}
