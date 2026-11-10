"""Full quotient factor maps and distinctness for the new ramified lift pilots."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path

from audit_guided_native import audit_guided_native
from profile_mechanism_breaker import install_bases,ramified_double_profile
from run_smooth_degree9 import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--job',type=Path,required=True)
    p.add_argument('--base-catalog',type=Path,required=True);p.add_argument('--auditor',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args();args.out.mkdir()
    prior=json.loads(args.base_catalog.read_text())['components'];install_bases(prior)
    old={(e['degree'],tuple(tuple(p) for p in e['component_identity_key'])) for e in prior if 'component_identity_key' in e}
    report=audit_guided_native(args.job,args.auditor,collect_records=True,native_components=True)
    report.pop('_closed_records');profiles=[];identities=set()
    for c,item in zip([c for c in report['components'] if c['summary']['status']=='closed_smooth'],report['closed_identities']):
        identity=(item['degree'],tuple(tuple(p) for p in item['key']))
        if identity in old or identity in identities:raise ValueError('not an additional distinct closed identity')
        identities.add(identity);n=item['degree'];g=int(c['summary']['genus'])
        # Original curated entries have exact genus; the stopped root entries
        # were compared by complete independent minimum-state identities above.
        if any(e['degree']==n and e.get('genus')==g for e in prior if 'component_identity_key' not in e):
            raise ValueError('old curated same-genus population needs full identity comparison')
        folder=args.job/item['folder']
        entry=dict(degree=n,states=int(c['summary']['states']),genus=g,
            representative=c['summary']['seed'],raw_folder=str(folder),table_sha256=c['table_sha256'],
            summary_sha256=sha(folder/'SUMMARY.tsv'),component_identity_key=item['key'])
        value=ramified_double_profile(entry);value.update(summary_sha256=entry['summary_sha256'],component_identity_key=item['key'])
        if value['group_order']!=int(c['summary']['group_order']):raise ValueError('independent point order differs')
        profiles.append(value);print('RAMIFIED',n,g,value['states'],value['group_order'],value['quotient_projection_fiber'],flush=True)
    (args.out/'CATALOG.json').write_text(json.dumps(dict(components=profiles,source_checks=report,
        distinct_from_prior_closed_identities=True,base_catalog_sha256=sha(args.base_catalog),
        components_added=len(profiles),states_added=sum(e['states'] for e in profiles),
        source_job=str(args.job),verifier_sha256=sha(args.auditor),
        novelty='New relative ramification pattern in this project, not a literature novelty or primitive/weak/geometric claim.'),indent=2)+'\n')

if __name__=='__main__':main()
