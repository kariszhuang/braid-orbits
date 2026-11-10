"""Independent tree--cotree orientation and exact small gauge comparisons."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import itertools
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from genus_three_intrinsic import IntrinsicFrame,PERMS,MUL,INV,IDENTITY,TRANS,cycles,compose,uniform3_quotient


def direct_keys(beta,gamma,active):
    """All gamma-pure voltages, before residual beta-cycle gauge reduction."""
    q=len(beta);frame=IntrinsicFrame(beta,gamma);white=cycles(gamma)
    roots={c[0]for c in white};free=[j for j in range(q)if j not in roots]
    labels={x:j for j,row in enumerate(cycles(beta))for x in row}
    gauges=list(itertools.product(range(6),repeat=len(cycles(beta))))
    def key(local,voltage):
        result=None
        for h in gauges:
            value=tuple(MUL[MUL[h[labels[j]]][local[j]]][INV[h[labels[j]]]]for j in range(q))
            value+=tuple(MUL[MUL[h[labels[gamma[j]]]][voltage[j]]][INV[h[labels[j]]]]for j in range(q))
            if result is None or value<result:result=value
        return result
    direct=set();intrinsic=set();raw_direct=0;raw_intrinsic=0
    for selected in itertools.product(TRANS,repeat=len(active)):
        local=[IDENTITY]*q
        for j,t in zip(active,selected):local[j]=t
        for values in itertools.product(range(6),repeat=len(free)):
            z=[IDENTITY]*q
            for j,t in zip(free,values):z[j]=t
            voltage=[MUL[INV[z[gamma[j]]]][z[j]]for j in range(q)]
            good=True
            for face in frame.faces:
                r=IDENTITY
                for j in face:r=MUL[MUL[local[frame.delta[j]]][voltage[j]]][r]
                if r!=IDENTITY:good=False;break
            if good:raw_direct+=1;direct.add(key(local,voltage))
        for handles in itertools.product(range(6),repeat=2):
            lifted=frame.solve(local,handles)
            if lifted is None:continue
            raw_intrinsic+=1
            voltage=[PERMS.index(tuple(lifted[2][3*j+x]%3 for x in range(3)))for j in range(q)]
            intrinsic.add(key(local,voltage))
    return direct,intrinsic,raw_direct,raw_intrinsic


class IntrinsicTests(unittest.TestCase):
    def test_known_one_branch_positive_genus_two_has_nine_handle_lifts(self):
        beta=uniform3_quotient(3);frame=IntrinsicFrame(beta,beta)
        local=[PERMS.index((1,2,0)),IDENTITY,IDENTITY]
        found=[t for h in itertools.product(range(6),repeat=2)if(t:=frame.solve(local,h))is not None]
        self.assertEqual(len(found),9)
        for t in found:
            self.assertEqual(sum(9-len(cycles(p))for p in t),20)
            self.assertEqual(compose(compose(t[0],t[1]),compose(t[2],t[3])),tuple(range(9)))

    def test_all_six_point_pure_voltage_classes_equal_intrinsic_handles(self):
        beta=uniform3_quotient(6);gamma=(1,3,4,0,5,2)
        # Complement of one opposite pair in the A4 edge action.
        active=(0,1,3,5)
        direct,intrinsic,raw_direct,raw_intrinsic=direct_keys(beta,gamma,active)
        self.assertEqual(direct,intrinsic);self.assertTrue(direct)
        self.assertEqual(raw_direct,6*raw_intrinsic)
        print('INTRINSIC_GAUGE',len(direct),raw_direct,raw_intrinsic)

    def test_unsupported_surface_and_disconnected_data_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'genus one'):IntrinsicFrame((1,0,2),(0,2,1))
        beta=uniform3_quotient(6)
        with self.assertRaisesRegex(ValueError,'disconnected'):IntrinsicFrame(beta,beta)

    def test_all_profile_beta_reduction_has_exact_small_lists(self):
        from verify_genus_three_intrinsic import elliptic_betas
        profile=lambda p:tuple(sorted(map(len,cycles(p))))
        self.assertEqual([profile(p)for p in elliptic_betas(7)],[])
        self.assertEqual([profile(p)for p in elliptic_betas(8)],[(4,4)])
        self.assertEqual([profile(p)for p in elliptic_betas(9)],[(1,4,4),(2,3,4),(3,3,3)])
        self.assertEqual(len(elliptic_betas(10)),5)
        self.assertEqual(len(elliptic_betas(11)),4)


@unittest.skipUnless(os.environ.get('G3_INTRINSIC_CONTROL'),'requires exact finite intrinsic control')
class IntrinsicAuditTests(unittest.TestCase):
    def test_independent_complete_coordinates_and_adversarial_markings(self):
        from verify_genus_three_intrinsic import verify
        source=Path(os.environ['G3_INTRINSIC_CONTROL'])
        self.assertTrue(verify(source)['passed'])
        for name,old,new in [('COORDINATES.tsv','\t1,2\t','\t1,4\t'),('QUOTIENTS.tsv','\t1,2\n','\t2,1\n'),('POPULATION.tsv','\ttrue\t','\tfalse\t')]:
            with tempfile.TemporaryDirectory()as tmp:
                target=Path(tmp)/'control';shutil.copytree(source,target)
                path=target/name;value=path.read_text();changed=value.replace(old,new,1)
                self.assertNotEqual(value,changed,'corruption fixture must actually change')
                path.write_text(changed)
                with self.assertRaises((ValueError,AssertionError)):verify(target)


if __name__=='__main__':unittest.main()
