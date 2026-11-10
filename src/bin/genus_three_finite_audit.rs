//! Independent compact-certificate replay. NO producer/constructor imports.
//! Python separately supplies exact connection counts and frame inequivalence.
use std::{
    collections::HashSet,
    fs,
    io::{BufRead, BufReader},
    path::Path,
    time::Instant,
};
const N: usize = 15;
type P = [u8; N];
type T = [P; 4];
fn identity() -> P {
    std::array::from_fn(|i| i as u8)
}
fn multiply(a: P, b: P) -> P {
    std::array::from_fn(|i| a[b[i] as usize])
}
fn inverse(a: P) -> P {
    let mut b = identity();
    for i in 0..N {
        b[a[i] as usize] = i as u8;
    }
    b
}
fn cycles(a: P) -> Vec<usize> {
    let mut seen = [false; N];
    let mut lengths = vec![];
    for root in 0..N {
        if seen[root] {
            continue;
        }
        let (mut x, mut len) = (root, 0);
        while !seen[x] {
            seen[x] = true;
            len += 1;
            x = a[x] as usize;
        }
        lengths.push(len);
    }
    lengths.sort_unstable();
    lengths
}
fn index(a: P) -> usize {
    let mut seen = 0u16;
    let mut count = 0;
    for root in 0..N {
        if seen & (1 << root) != 0 {
            continue;
        }
        count += 1;
        let mut x = root;
        while seen & (1 << x) == 0 {
            seen |= 1 << x;
            x = a[x] as usize;
        }
    }
    N - count
}
fn smooth(t: T) -> bool {
    let k = index(multiply(t[0], t[1]));
    k == index(t[0]) + index(t[1]) || k == index(t[2]) + index(t[3])
}
fn move_tuple(mut t: T, m: u8) -> T {
    let i = m as usize % 3;
    let (a, b) = (t[i], t[i + 1]);
    if m < 3 {
        t[i] = multiply(multiply(a, b), inverse(a));
        t[i + 1] = a;
    } else {
        t[i] = b;
        t[i + 1] = multiply(multiply(inverse(b), a), b);
    }
    t
}
fn parse_p(raw: &str) -> Result<P, String> {
    let values = raw
        .split(',')
        .map(|s| s.parse::<u8>().map_err(|e| e.to_string()))
        .collect::<Result<Vec<_>, _>>()?;
    let p: P = values.try_into().map_err(|_| "permutation length")?;
    let mut q = p;
    q.sort_unstable();
    if q != identity() {
        return Err("invalid images".into());
    }
    Ok(p)
}
fn packed(p: P) -> u64 {
    p.iter()
        .enumerate()
        .fold(0, |s, (i, &x)| s | ((x as u64) << (4 * i)))
}
fn lines(path: &Path) -> Result<Vec<String>, String> {
    BufReader::new(fs::File::open(path).map_err(|e| e.to_string())?)
        .lines()
        .collect::<Result<_, _>>()
        .map_err(|e| e.to_string())
}
fn partition(raw: &str) -> Result<Vec<usize>, String> {
    let mut p = raw
        .split(',')
        .map(|s| s.parse::<usize>().map_err(|e| e.to_string()))
        .collect::<Result<Vec<_>, _>>()?;
    if p.iter().sum::<usize>() != N || p.contains(&0) {
        return Err("partition degree".into());
    }
    p.sort_unstable();
    Ok(p)
}
fn audit(folder: &Path, mu: Vec<usize>, nu: Vec<usize>, workers: usize) -> Result<String, String> {
    let began = Instant::now();
    let mut frames = vec![];
    let source = lines(&folder.join("FRAMES.tsv"))?;
    if source.first().map(String::as_str) != Some("frame\ta\tb") {
        return Err("frame header".into());
    }
    let special = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4];
    let p: P = std::array::from_fn(|i| ((i + 1) % N) as u8);
    for line in &source[1..] {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 3 || fields[0].parse::<usize>().ok() != Some(frames.len()) {
            return Err("frame rank/columns".into());
        }
        let (a, b) = (parse_p(fields[1])?, parse_p(fields[2])?);
        if cycles(a) != special || cycles(b).len() != 5 || multiply(a, b) != p {
            return Err("frame type/product".into());
        }
        frames.push([a, b]);
    }
    let mut completions = vec![];
    let mut distinct = HashSet::new();
    let r = inverse(p);
    let source = lines(&folder.join("COMPLETIONS.tsv"))?;
    if source.first().map(String::as_str) != Some("completion\tc") {
        return Err("completion header".into());
    }
    for line in &source[1..] {
        let (rank, raw) = line.split_once('\t').ok_or("completion columns")?;
        if rank.parse::<usize>().ok() != Some(completions.len()) {
            return Err("completion rank".into());
        }
        let c = parse_p(raw)?;
        let d = multiply(inverse(c), r);
        if cycles(c) != mu || cycles(d) != nu || !distinct.insert(packed(c)) {
            return Err("completion type/product/duplicate".into());
        }
        completions.push([c, d]);
    }
    drop(source);
    drop(distinct);
    let source = lines(&folder.join("WORDS.tsv"))?;
    let mut words = vec![vec![]];
    if source.first().map(String::as_str) != Some("token\tword") {
        return Err("dictionary header".into());
    }
    for line in &source[1..] {
        let (rank, raw) = line.split_once('\t').ok_or("dictionary columns")?;
        if rank.parse::<usize>().ok() != Some(words.len()) {
            return Err("dictionary numbering".into());
        }
        let word = raw
            .split(',')
            .map(|s| s.parse::<u8>().map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        if word.is_empty() || word.len() > 32 || word.iter().any(|&m| m >= 6) {
            return Err("invalid witness word".into());
        }
        words.push(word);
    }
    let tokens = fs::read(folder.join("WITNESSES.u16le")).map_err(|e| e.to_string())?;
    let states = frames
        .len()
        .checked_mul(completions.len())
        .ok_or("state count overflow")?;
    if states == 0 || tokens.len() != states.checked_mul(2).ok_or("byte count overflow")? {
        return Err("witness byte length".into());
    }
    let parse_seconds = began.elapsed().as_secs_f64();
    let now = Instant::now();
    let chunk = frames.len().div_ceil(workers).max(1) * completions.len() * 2;
    let hist = std::thread::scope(|scope| -> Result<Vec<u64>, String> {
        let mut handles = vec![];
        for (part, bytes) in tokens.chunks(chunk).enumerate() {
            let frames = &frames;
            let cs = &completions;
            let words = &words;
            handles.push(scope.spawn(move || -> Result<Vec<u64>, String> {
                let mut counts = vec![0u64; words.len()];
                let base = part * chunk / 2;
                for (frame_offset, block) in bytes.chunks(cs.len() * 2).enumerate() {
                    let pair = frames[base / cs.len() + frame_offset];
                    for (i, (token, cd)) in block.chunks_exact(2).zip(cs).enumerate() {
                        let token = u16::from_le_bytes([token[0], token[1]]) as usize;
                        if token >= words.len() {
                            return Err("invalid proof dictionary token".into());
                        }
                        counts[token] += 1;
                        if token == 0 {
                            continue;
                        }
                        let mut t = [pair[0], pair[1], cd[0], cd[1]];
                        for &m in &words[token] {
                            t = move_tuple(t, m);
                        }
                        if smooth(t) {
                            return Err(format!(
                                "false negative witness at state {}",
                                base + frame_offset * cs.len() + i
                            ));
                        }
                    }
                }
                Ok(counts)
            }));
        }
        let mut total = vec![0u64; words.len()];
        for handle in handles {
            let h = handle.join().map_err(|_| "audit worker panic")??;
            for (i, c) in h.into_iter().enumerate() {
                total[i] += c;
            }
        }
        Ok(total)
    })?;
    let seconds = now.elapsed().as_secs_f64();
    Ok(format!("{{\"frames\":{},\"completions\":{},\"states\":{states},\"token_histogram\":{hist:?},\"workers\":{workers},\"parse_seconds\":{parse_seconds},\"replay_seconds\":{seconds}}}",frames.len(),completions.len()))
}
fn main() {
    let a = std::env::args().collect::<Vec<_>>();
    let result = (|| {
        if a.len() != 5 {
            return Err("usage: finite_audit FOLDER MU_CSV NU_CSV WORKERS".into());
        }
        let workers = a[4].parse::<usize>().map_err(|e| e.to_string())?;
        if workers == 0 || workers > 10 {
            return Err("workers1..10".into());
        }
        audit(
            Path::new(&a[1]),
            partition(&a[2])?,
            partition(&a[3])?,
            workers,
        )
    })();
    match result {
        Ok(text) => println!("{text}"),
        Err(e) => {
            eprintln!("finite-audit: {e}");
            std::process::exit(1);
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_moves_preserve_product_and_invert() {
        let mut a = identity();
        a.swap(0, 1);
        let mut b = identity();
        b[0] = 1;
        b[1] = 2;
        b[2] = 0;
        let mut c = identity();
        c.swap(2, 3);
        let d = inverse(multiply(multiply(a, b), c));
        let t = [a, b, c, d];
        for m in 0..6 {
            let u = move_tuple(t, m);
            assert_eq!(move_tuple(u, (m + 3) % 6), t);
            assert_eq!(
                multiply(multiply(u[0], u[1]), multiply(u[2], u[3])),
                identity()
            );
        }
    }
}
