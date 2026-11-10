"""Apply Proposition 1 of the supplied proof, never its converse.

This checks exact finite hypotheses. The geometric deduction over C remains
theorem-based, not a proof-assistant or numerical moduli computation.
"""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
from verify_smooth_geometry import parse_tuple, cycles

PROOF_SHA256 = '6b1e2a3d96ec76e6f62d0ff51be1d02c10b453339b9ce0933e7f782352ba380d'


def classify(entry):
    n = entry['degree']
    t = parse_tuple('[' + ','.join(entry['tuple']) + ']', n)
    witnesses = []
    failed_systems = []
    point_cycles = [len(cycles(p)) for p in t]
    for bi, b in enumerate(entry['block_quotients']):
        m, r = b['block_size'], b['block_count']
        q = parse_tuple('[' + ','.join(b['tuple']) + ']', r)
        labels = {x - 1: i for i, block in enumerate(b['blocks']) for x in block}
        quotient_cycles = [len(cycles(p)) for p in q]
        a_passing_entries = []
        candidate_failures = []
        for k in range(4):
            if q[k] != tuple(range(r)):
                continue
            a_passing_entries.append(k + 1)
            active = sorted({labels[x] for x in range(n) if t[k][x] != x})
            cycle_failures = [dict(entry=j + 1, actual_cycles=point_cycles[j],
                quotient_cycles=quotient_cycles[j], required_cycles=m * quotient_cycles[j],
                ramification_excess=m * quotient_cycles[j] - point_cycles[j])
                for j in range(4) if j != k and point_cycles[j] != m * quotient_cycles[j]]
            failed = (['b'] if len(active) != 1 else []) + (['c'] if cycle_failures else [])
            if failed:
                candidate_failures.append(dict(special_entry=k + 1,
                    failed_conditions=failed, active_blocks=[x + 1 for x in active],
                    ordinary_cycle_count_failures=cycle_failures))
                continue
            h = b['genus']
            if h < 1:
                raise ValueError('Proposition hypotheses contradict quotient genus >= 1')
            witnesses.append(dict(block_quotient_index=bi, special_entry=k + 1,
                active_block=active[0] + 1, block_count=r, block_size=m,
                quotient_genus=h, claim=6 if h == 1 else 7,
                hypotheses=dict(quotient_entry_identity=True, active_blocks=1,
                    ordinary_cycle_count_equality=True)))
        failed_systems.append(dict(block_quotient_index=bi, block_count=r, block_size=m,
            quotient_cycle_counts=quotient_cycles, a_passing_entries=a_passing_entries,
            candidate_failures=candidate_failures))
    decisions = {w['quotient_genus'] == 1 for w in witnesses}
    if len(decisions) > 1:
        raise ValueError('contradictory geometric witnesses: review proof/data')
    value = next(iter(decisions)) if decisions else None
    explanation = ('Proposition 1(a)-(c) holds for a recorded invariant block system. '
        + ('Claim 6: the fixed quotient has genus 1; translating the unique moving branch point gives an isotrivial source family.'
           if value else 'Claim 7: the fixed quotient has genus at least 2; the unique branch point moves, so de Franchis finiteness rules out an isotrivial source family.')) if decisions else (
        'No recorded block system satisfies all of Proposition 1(a)-(c). '
        'The supplied sufficient criterion does not decide isotriviality; this is not evidence of non-isotriviality.')
    result = dict(isotrivial=value, isotriviality=dict(
        status='isotrivial' if value is True else 'non-isotrivial' if value is False else 'unknown',
        method='proposition-1-single-active-block-v1',
        source='/sources/isotrivial-proof.pdf',
        source_sha256=PROOF_SHA256,
        explanation=explanation, witnesses=witnesses))
    if value is None:
        result['isotriviality']['failed_block_systems'] = failed_systems
    return result
