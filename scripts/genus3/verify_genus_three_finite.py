#!/usr/bin/env python3
"""Independent literal-word and character-count audit of the finite S15 targets.

No producer imports, no inferred negative from a cap/geometry test. A full
negative witness for EVERY distinct necessary state certifies the passport.
Safe mixed cores require a separate full-B4 positive audit and are rejected here.
"""
import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import struct
import subprocess
import time

if not __debug__:
    raise RuntimeError('Do not run the independent verifier with Python -O')

N = 15
IDENTITY = tuple(range(N))
A = (3, 3, 3, 3, 3)
B = (2, 3, 3, 3, 4)
SPECIAL = (1,) * 9 + (2, 4)


def mul(a, b):
    return tuple(a[x] for x in b)


def inv(a):
    q = [0] * len(a)
    for i, x in enumerate(a):
        q[x] = i
    return tuple(q)


def profile(a):
    seen = set()
    lengths = []
    for root in range(len(a)):
        if root in seen:
            continue
        x, length = root, 0
        while x not in seen:
            seen.add(x)
            length += 1
            x = a[x]
        lengths.append(length)
    return tuple(sorted(lengths))


def smooth(t):
    indices = [N - len(profile(p)) for p in t]
    product_index = N - len(profile(mul(t[0], t[1])))
    return product_index in (sum(indices[:2]), sum(indices[2:]))


def braid(t, move):
    t = list(t)
    i = move % 3
    a, b = t[i:i + 2]
    if move < 3:
        t[i:i + 2] = mul(mul(a, b), inv(a)), a
    else:
        t[i:i + 2] = b, mul(mul(inv(b), a), b)
    return tuple(t)


def hooks(partition):
    # det(1+t P)/(1+t), from the natural permutation representation.
    n = sum(partition)
    poly = [1] + [0] * n
    for r in partition:
        for k in range(n, r - 1, -1):
            poly[k] -= (-1) ** r * poly[k - r]
    out = [poly[0]]
    for k in range(1, n):
        out.append(poly[k] - out[-1])
    assert poly[n] == out[-1]
    return out


def connection(mu, nu):
    n = sum(mu)
    assert n == sum(nu)
    def class_size(partition):
        z = 1
        for length, count in Counter(partition).items():
            z *= length ** count * factorial(count)
        return factorial(n) // z
    h, k = hooks(mu), hooks(nu)
    value = Fraction(class_size(mu) * class_size(nu), factorial(n)) * sum(
        Fraction((-1) ** j * h[j] * k[j], comb(n - 1, j)) for j in range(n))
    assert value.denominator == 1 and value >= 0
    return value.numerator


def pair_key(a, b):
    # Independent root traversal along ab, not the producer's frame lookup.
    p = mul(a, b)
    candidates = []
    for root in range(N):
        order = [root]
        for _ in range(N - 1):
            order.append(p[order[-1]])
        assert len(set(order)) == N and p[order[-1]] == root
        labels = {x: i for i, x in enumerate(order)}
        candidates.append(tuple(labels[q[x]] for q in (a, b) for x in order))
    # The unique 4-cycle cannot be preserved by a nonidentity power of p.
    assert len(set(candidates)) == N, 'nontrivial first-pair centralizer'
    return min(candidates)


def audit(folder, native=None, workers=1):
    start = time.monotonic()
    folder = Path(folder)
    summary = json.loads((folder / 'SUMMARY.json').read_text())
    if summary['schema'] == 2:
        return audit_slice(folder, summary, native, workers)
    mixed = summary['population'].startswith('degree15 mixed ')
    assert mixed or summary['population'].startswith('degree15 pure ')
    assert summary['safe_mixed_states'] == 0, 'positive core needs independent full braid closure'
    pairs, frame_keys, sets = {}, set(), {}
    counts, words = Counter(), Counter()
    block_counts, total, max_word = Counter(), 0, 0
    digest = hashlib.sha256()
    with (folder / 'STATES.tsv').open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            assert int(row['state']) == total, 'nonsequential state'
            t = tuple(tuple(map(int, row[k].split(','))) for k in ('a', 'b', 'c', 'd'))
            assert all(len(p) == N and sorted(p) == list(IDENTITY) for p in t), 'invalid images'
            assert mul(mul(t[0], t[1]), mul(t[2], t[3])) == IDENTITY, 'product'
            profiles = tuple(profile(p) for p in t)
            assert profiles[0] == SPECIAL
            assert sorted(profiles[1:]) == sorted((A, A, B) if mixed else (A, A, A)), 'passport'
            assert profile(mul(t[0], t[1])) == (N,), 'not necessary full-product state'
            assert (sum(N - len(p) for p in profiles) - 2*N + 2) // 2 == 3
            frame = int(row['frame'])
            pair = t[:2]
            if frame not in pairs:
                key = pair_key(*pair)
                assert key not in frame_keys, 'conjugate first-pair frames'
                frame_keys.add(key)
                pairs[frame] = pair
                sets[frame] = set()
            assert pairs[frame] == pair
            assert t[2] not in sets[frame], 'duplicate state in centralizer-free frame'
            sets[frame].add(t[2])
            counts[(frame, profiles[2], profiles[3])] += 1
            assert row['status'] == 'nonsmooth'
            word = tuple(map(int, row['word'].split(','))) if row['word'] else ()
            assert word and all(0 <= m < 6 for m in word), 'invalid or empty negative word'
            u = t
            for m in word:
                u = braid(u, m)
            assert not smooth(u), f'false negative witness at state {total}'
            max_word = max(max_word, len(word))
            words[word] += 1
            # ab is a full cycle, so these are ALL possible block systems.
            p = mul(t[0], t[1])
            order = [0]
            for _ in range(N - 1):
                order.append(p[order[-1]])
            labels = {x: i for i, x in enumerate(order)}
            block = None
            for k in (3, 5):
                if all(all(labels[q[x]] % k == labels[q[order[labels[x] % k]]] % k
                           for x in order) for q in t):
                    block = k
                    break
            block_counts[block] += 1
            digest.update((';'.join(','.join(map(str, p)) for p in t) + '\n').encode())
            total += 1
    # Exact character counts + distinct centralizer-free states prove coverage;
    # arithmetic totals are not being used as a substitute for word replay.
    expected_frames = {A: connection(SPECIAL, A) // N}
    if mixed:
        expected_frames[B] = connection(SPECIAL, B) // N
    assert Counter(profile(pair[1]) for pair in pairs.values()) == Counter(expected_frames)
    assert sorted(pairs) == list(range(len(pairs)))
    for frame, pair in pairs.items():
        remaining = list((A, A, B) if mixed else (A, A, A))
        remaining.remove(profile(pair[1]))
        orientations = set((tuple(remaining), tuple(reversed(remaining))))
        for mu, nu in orientations:
            assert counts[(frame, mu, nu)] == connection(mu, nu), 'incomplete factorization population'
        assert sum(counts[(frame, mu, nu)] for mu, nu in orientations) == len(sets[frame])
    assert total == summary['candidate_states'] == summary['nonsmooth']
    assert block_counts[None] == summary['primitive_states']
    assert total - block_counts[None] == summary['imprimitive']
    assert max_word == summary['maximum_negative_word_length']
    with (folder/'STATES.tsv').open('rb') as file:
        states_sha = hashlib.file_digest(file,'sha256').hexdigest()
    return dict(status='verified_complete_negative_passport', population=summary['population'],
                quotient='simultaneous S15 conjugation', states=total,
                first_pair_frames=len(pairs), primitive=block_counts[None],
                three_block_states=block_counts[3], five_block_states=block_counts[5],
                unresolved=0, max_negative_word=max_word,
                word_histogram=[dict(word=w, states=n) for w, n in sorted(words.items())],
                tuple_stream_sha256=digest.hexdigest(),
                states_sha256=states_sha,
                seconds=time.monotonic()-start)


def audit_slice(folder, summary, native=None, workers=1):
    start = time.monotonic()
    ordinary = tuple(sorted(summary['first_ordinary_profile']))
    last = tuple(sorted(summary.get('remaining_last_profile', A)))
    first = tuple(sorted(summary.get('remaining_first_profile', A)))
    assert len(first) == len(last) == 5 and sum(first) == sum(last) == N
    frames = []
    frame_keys = set()
    with (folder/'FRAMES.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            assert int(row['frame']) == len(frames)
            pair = tuple(tuple(map(int, row[k].split(','))) for k in ('a', 'b'))
            assert all(len(p) == N and sorted(p) == list(IDENTITY) for p in pair)
            assert profile(pair[0]) == SPECIAL and profile(pair[1]) == ordinary
            assert mul(*pair) == tuple((i+1) % N for i in range(N)), 'not standard full product'
            key = pair_key(*pair)
            assert key not in frame_keys
            frame_keys.add(key)
            frames.append(pair)
    assert len(frames) == summary['frames'] == connection(SPECIAL, ordinary)//N
    if native:
        return audit_slice_native(folder, summary, frames, ordinary, first, last, native, workers, start)
    completions = []
    r = tuple((i-1) % N for i in range(N))
    with (folder/'COMPLETIONS.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            assert int(row['completion']) == len(completions)
            c = tuple(map(int, row['c'].split(',')))
            assert len(c) == N and sorted(c) == list(IDENTITY)
            d = mul(inv(c), r)
            assert profile(c) == first and profile(d) == last
            completions.append((c, d))
    assert len(set(completions)) == len(completions) == connection(first, last) == summary['completions']
    words = [()]
    with (folder/'WORDS.tsv').open() as f:
        for row in csv.DictReader(f, delimiter='\t'):
            assert int(row['token']) == len(words)
            w = tuple(map(int, row['word'].split(',')))
            assert w and all(0 <= m < 6 for m in w)
            words.append(w)
    data = (folder/'WITNESSES.u16le').read_bytes()
    assert len(data) == len(frames)*len(completions)*2
    histogram = Counter()
    unresolved = 0
    for state, (token,) in enumerate(struct.iter_unpack('<H', data)):
        assert token < len(words), 'invalid proof dictionary token'
        if not token:
            unresolved += 1
            continue
        frame, ci = divmod(state, len(completions))
        a, b = frames[frame]
        c, d = completions[ci]
        t = (a, b, c, d)
        for m in words[token]:
            t = braid(t, m)
        assert not smooth(t), f'false negative witness at compact state {state}'
        histogram[words[token]] += 1
    total = len(frames)*len(completions)
    assert total == summary['candidate_states']
    assert unresolved == summary['unresolved'] and total-unresolved == summary['nonsmooth']
    assert max(map(len, histogram), default=0) == summary['maximum_negative_word_length']
    files = ('FRAMES.tsv', 'COMPLETIONS.tsv', 'WORDS.tsv', 'WITNESSES.u16le', 'SUMMARY.json')
    hashes = {}
    for name in files:
        with (folder/name).open('rb') as f:
            hashes[name] = hashlib.file_digest(f, 'sha256').hexdigest()
    return dict(status='verified_complete_negative_passport' if not unresolved else 'verified_partial_no_exclusion',
                population=summary['population'], first_ordinary_profile=ordinary,
                quotient='simultaneous S15 conjugation', selected_ordered_layout_states=total,
                negative=total-unresolved, unresolved=unresolved,
                max_negative_word=summary['maximum_negative_word_length'],
                word_histogram=[dict(word=w, states=n) for w, n in sorted(histogram.items())],
                hashes=hashes, seconds=time.monotonic()-start)


def audit_slice_native(folder, summary, frames, ordinary, first, last, native, workers, start):
    assert 1 <= workers <= 10
    command = [str(native), str(folder), ','.join(map(str,first)), ','.join(map(str,last)), str(workers)]
    process = subprocess.run(command, capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    result = json.loads(process.stdout)
    assert result['frames'] == len(frames)
    assert result['completions'] == summary['completions'] == connection(first,last)
    total = len(frames)*result['completions']
    assert result['states'] == total == summary['candidate_states']
    words = [()]
    with (folder/'WORDS.tsv').open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            assert int(row['token']) == len(words)
            word = tuple(map(int,row['word'].split(',')))
            assert word and all(0 <= m < 6 for m in word)
            words.append(word)
    histogram = result['token_histogram']
    assert len(histogram) == len(words) and sum(histogram) == total
    unresolved = histogram[0]
    assert unresolved == summary['unresolved'] and total-unresolved == summary['nonsmooth']
    maximum = max((len(w) for w,c in zip(words,histogram) if c),default=0)
    assert maximum == summary['maximum_negative_word_length']
    hashes = {}
    for name in ('FRAMES.tsv','COMPLETIONS.tsv','WORDS.tsv','WITNESSES.u16le','SUMMARY.json'):
        with (folder/name).open('rb') as f:
            hashes[name] = hashlib.file_digest(f,'sha256').hexdigest()
    with Path(native).open('rb') as f:
        native_sha = hashlib.file_digest(f,'sha256').hexdigest()
    return dict(status='verified_complete_negative_passport' if not unresolved else 'verified_partial_no_exclusion',
                population=summary['population'], first_ordinary_profile=ordinary,
                quotient='simultaneous S15 conjugation', selected_ordered_layout_states=total,
                negative=total-unresolved, unresolved=unresolved, max_negative_word=maximum,
                word_histogram=[dict(word=w,states=c) for w,c in zip(words[1:],histogram[1:]) if c],
                independent_native=result, native_sha256=native_sha, hashes=hashes,
                seconds=time.monotonic()-start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--native', type=Path, help='separately implemented literal replay; no producer imports')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    result = audit(args.folder,args.native,args.workers)
    text = json.dumps(result, indent=2) + '\n'
    if args.out:
        with args.out.open('x') as f:
            f.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
