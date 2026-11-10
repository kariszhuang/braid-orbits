use std::collections::{HashMap, HashSet};

fn factorial(n: usize) -> u128 {
    (1..=n as u128).product()
}

fn partitions(n: usize, limit: usize, current: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if n == 0 {
        out.push(current.clone());
        return;
    }
    for part in (1..=n.min(limit)).rev() {
        current.push(part);
        partitions(n - part, part, current, out);
        current.pop();
    }
}

fn dimension(shape: &[usize]) -> u128 {
    let n: usize = shape.iter().sum();
    let hooks: u128 = shape
        .iter()
        .enumerate()
        .flat_map(|(r, &length)| {
            (0..length).map(move |c| {
                (length - c + shape.iter().skip(r + 1).filter(|&&other| other > c).count()) as u128
            })
        })
        .product();
    factorial(n) / hooks
}

fn subshapes(
    shape: &[usize],
    target: usize,
    row: usize,
    max: usize,
    current: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if row == shape.len() {
        if target == 0 {
            out.push(current.clone());
        }
        return;
    }
    for length in (0..=shape[row].min(max).min(target)).rev() {
        if length > 0 {
            current.push(length);
        }
        subshapes(shape, target - length, row + 1, length, current, out);
        if length > 0 {
            current.pop();
        }
    }
}

fn rim_height(shape: &[usize], sub: &[usize]) -> Option<usize> {
    let mut cells = HashSet::new();
    for (row, &length) in shape.iter().enumerate() {
        for col in sub.get(row).copied().unwrap_or(0)..length {
            cells.insert((row, col));
        }
    }
    if cells.is_empty() {
        return None;
    }
    if cells.iter().any(|&(r, c)| {
        cells.contains(&(r + 1, c))
            && cells.contains(&(r, c + 1))
            && cells.contains(&(r + 1, c + 1))
    }) {
        return None;
    }
    let start = *cells.iter().next().unwrap();
    let mut reached = HashSet::from([start]);
    let mut frontier = vec![start];
    while let Some((r, c)) = frontier.pop() {
        for next in [
            (r.wrapping_sub(1), c),
            (r + 1, c),
            (r, c.wrapping_sub(1)),
            (r, c + 1),
        ] {
            if cells.contains(&next) && reached.insert(next) {
                frontier.push(next);
            }
        }
    }
    if reached.len() != cells.len() {
        return None;
    }
    let rows: HashSet<_> = cells.iter().map(|&(r, _)| r).collect();
    Some(rows.len() - 1)
}

fn character(
    shape: &[usize],
    cycles: &[usize],
    index: usize,
    memo: &mut HashMap<(Vec<usize>, usize), i128>,
) -> i128 {
    if index == cycles.len() {
        return i128::from(shape.is_empty());
    }
    if cycles[index..].iter().all(|&c| c == 1) {
        return dimension(shape) as i128;
    }
    let key = (shape.to_vec(), index);
    if let Some(&answer) = memo.get(&key) {
        return answer;
    }
    let target = shape.iter().sum::<usize>() - cycles[index];
    let mut subs = Vec::new();
    subshapes(shape, target, 0, shape[0], &mut Vec::new(), &mut subs);
    let answer = subs
        .into_iter()
        .filter_map(|sub| {
            rim_height(shape, &sub).map(|height| {
                let sign = if height % 2 == 0 { 1 } else { -1 };
                sign * character(&sub, cycles, index + 1, memo)
            })
        })
        .sum();
    memo.insert(key, answer);
    answer
}

fn power_mod(mut base: u128, mut exponent: u128, modulus: u128) -> u128 {
    let mut answer = 1;
    while exponent > 0 {
        if exponent & 1 == 1 {
            answer = answer * base % modulus;
        }
        base = base * base % modulus;
        exponent >>= 1;
    }
    answer
}

fn product_one_count(n: usize, involutions: usize, fourth: &[usize]) -> u128 {
    let mut shapes = Vec::new();
    partitions(n, n, &mut Vec::new(), &mut shapes);
    let alpha = vec![2; involutions];
    let mut memo_alpha = HashMap::new();
    let mut memo_beta = HashMap::new();
    let terms: Vec<_> = shapes
        .iter()
        .map(|shape| {
            (
                character(shape, &alpha, 0, &mut memo_alpha),
                character(shape, fourth, 0, &mut memo_beta),
                dimension(shape),
            )
        })
        .collect();
    let mut multiplicities = HashMap::<usize, usize>::new();
    for &length in fourth {
        *multiplicities.entry(length).or_default() += 1;
    }
    let centralizer_beta: u128 = multiplicities
        .iter()
        .map(|(&length, &count)| (length as u128).pow(count as u32) * factorial(count))
        .product();
    let class_alpha = factorial(n) / ((2_u128).pow(involutions as u32) * factorial(involutions));
    let class_beta = factorial(n) / centralizer_beta;
    // Frobenius's class-algebra formula, evaluated exactly by three prime moduli.
    let sample_space = class_alpha.checked_pow(3).unwrap();
    let (mut total, mut modulus) = (0_u128, 1_u128);
    for prime in [1_000_000_007_u128, 1_000_000_009, 998_244_353] {
        let character_sum = terms.iter().fold(0_u128, |sum, &(a, b, dim)| {
            assert_ne!(dim % prime, 0);
            let a = a.rem_euclid(prime as i128) as u128;
            let b = b.rem_euclid(prime as i128) as u128;
            let denominator = power_mod(dim % prime, 2, prime);
            (sum + power_mod(a, 3, prime) * b % prime * power_mod(denominator, prime - 2, prime)
                % prime)
                % prime
        });
        let residue = power_mod(class_alpha % prime, 3, prime) * (class_beta % prime) % prime
            * power_mod(factorial(n) % prime, prime - 2, prime)
            % prime
            * character_sum
            % prime;
        let step = ((residue + prime - total % prime) % prime)
            * power_mod(modulus % prime, prime - 2, prime)
            % prime;
        total += modulus * step;
        modulus *= prime;
    }
    assert!(modulus > sample_space && total <= sample_space);
    total
}

fn main() {
    let class_alpha = factorial(18) / (2_u128.pow(9) * factorial(9));
    let total = product_one_count(18, 9, &[4, 4, 4, 1, 1, 1, 1, 1, 1]);
    let component = 10 * factorial(18);
    assert_eq!(32 * total, 351 * component);
    println!("ordered passport product-one tuples (all groups) = {total}");
    println!(
        "passport tuples / 18! = {:.12}",
        total as f64 / factorial(18) as f64
    );
    println!("known 40-state component with exceptional entry last = {component}");
    println!(
        "probability of known component among product-one passport tuples = {:.12}",
        component as f64 / total as f64
    );
    println!(
        "passport probability from three random involutions = {:.12e}",
        total as f64 / class_alpha.pow(3) as f64
    );
    println!(
        "known-component probability from three random involutions = {:.12e}",
        component as f64 / class_alpha.pow(3) as f64
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matchings(points: &[usize], permutation: Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if points.is_empty() {
            out.push(permutation);
            return;
        }
        let a = points[0];
        for &b in &points[1..] {
            let mut next = permutation.clone();
            next[a] = b;
            next[b] = a;
            let rest: Vec<_> = points[1..].iter().copied().filter(|&x| x != b).collect();
            matchings(&rest, next, out);
        }
    }

    fn cycle_type(permutation: &[usize]) -> Vec<usize> {
        let mut seen = vec![false; permutation.len()];
        let mut lengths = Vec::new();
        for start in 0..permutation.len() {
            if seen[start] {
                continue;
            }
            let (mut at, mut length) = (start, 0);
            while !seen[at] {
                seen[at] = true;
                length += 1;
                at = permutation[at];
            }
            lengths.push(length);
        }
        lengths.sort_unstable_by(|a, b| b.cmp(a));
        lengths
    }

    #[test]
    fn small_character_and_class_count() {
        let mut memo = HashMap::new();
        assert_eq!(character(&[2, 1], &[3], 0, &mut memo), -1);
        assert_eq!(character(&[1, 1, 1], &[2, 1], 0, &mut memo), -1);
        // In S4 there are three fixed-point-free involutions. For triples,
        // every product is another such involution, never a 4-cycle.
        assert_eq!(product_one_count(4, 2, &[4]), 0);

        let mut small = Vec::new();
        matchings(&(0..6).collect::<Vec<_>>(), (0..6).collect(), &mut small);
        let mut brute = 0;
        for a in &small {
            for b in &small {
                for c in &small {
                    let product: Vec<_> = (0..6).map(|i| a[b[c[i]]]).collect();
                    brute += usize::from(cycle_type(&product) == [4, 1, 1]);
                }
            }
        }
        assert_eq!(product_one_count(6, 3, &[4, 1, 1]), brute as u128);

        let mut shapes = Vec::new();
        partitions(18, 18, &mut Vec::new(), &mut shapes);
        let alpha = vec![2; 9];
        let beta = [4, 4, 4, 1, 1, 1, 1, 1, 1];
        let mut a_memo = HashMap::new();
        let mut b_memo = HashMap::new();
        let a_norm: i128 = shapes
            .iter()
            .map(|s| character(s, &alpha, 0, &mut a_memo).pow(2))
            .sum();
        let b_norm: i128 = shapes
            .iter()
            .map(|s| character(s, &beta, 0, &mut b_memo).pow(2))
            .sum();
        let regular: i128 = shapes.iter().map(|s| (dimension(s) as i128).pow(2)).sum();
        assert_eq!(a_norm, (2_u128.pow(9) * factorial(9)) as i128);
        assert_eq!(
            b_norm,
            (4_u128.pow(3) * factorial(3) * factorial(6)) as i128
        );
        assert_eq!(regular, factorial(18) as i128);
        assert_eq!(
            product_one_count(18, 9, &beta),
            702_260_365_847_040_000
        );
    }
}
