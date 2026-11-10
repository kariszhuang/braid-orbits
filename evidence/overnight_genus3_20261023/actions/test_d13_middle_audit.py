#!/usr/bin/env python3
"""Regression for the actual producer/auditor marking mismatch."""
from itertools import combinations
import json
from pathlib import Path
import unittest
import audit_d13_double_middle as audit

class FrameMarkingTest(unittest.TestCase):
    def test_independent_marking_matches_frozen_producer_frame_set(self):
        raw=Path('/Volumes/T7/Research/genus3_d13_double_middle_20261023_0256/RESULT.json')
        expected={tuple(row['a'])for row in json.loads(raw.read_text())['frames']}
        p=tuple(range(1,13))+(0,)
        actual=set()
        for edges in combinations(combinations(range(13),2),2):
            if len(set(edges[0]+edges[1]))!=4:continue
            cc=audit.cycles(audit.mul(audit.match(edges),p))
            if len(cc)==3 and min(map(len,cc))>=3:
                actual.add(audit.match(audit.canonical_edges(edges)))
        self.assertEqual(actual,expected)

if __name__=='__main__':unittest.main()
