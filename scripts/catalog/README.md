# Public catalogue

Builds the public catalogue: classify isotriviality by Proposition 1, give each component a stable id, write versioned releases with deltas, package the website and preview it locally.

| file | what it does |
|---|---|
| `audit_mathematical_patterns.py` | Bounded, local-only audit of a frozen smooth catalog; no discovery/deployment. |
| `bind_public_component_ids.py` | Bootstrap stable IDs for the unversioned catalog from its complete certificates. |
| `build_public_smooth_catalog.py` | Publish the approved math-only catalog from a verified export; never starts searches. |
| `check_smooth_release.py` | Verify a versioned smooth publication locally, optionally against public URLs. |
| `classify_smooth_isotrivial.py` | Apply Proposition 1 of the supplied proof, never its converse. |
| `clean_smooth_catalog.py` | Produce a compact mathematical catalog for positive-genus smooth components. |
| `count_isotrivial_failures.py` | Count failures of Nitiwit's sufficient criterion; reads JSON, writes nothing. |
| `export_catalog_growth_public.py` | Publish a frozen selection of independently certified catalog-growth records. |
| `export_research_catalog.py` | Export positive-genus everywhere-smooth components as filterable math data. |
| `extend_public_smooth_catalog.py` | Merge independently certified new closures into an existing public catalog. |
| `package_current_math_release.py` | Package one current snapshot plus selected deltas, with no historical snapshots. |
| `package_math_ui.py` | Attach the checked static UI and pinned proof to a local release (no deploy). |
| `preview_math_catalog.py` | Local QA only: reproduce the public read-only URL mappings. |
| `summarize_smooth_unknowns.py` | Report exact failures of Proposition 1 for an annotated local catalog. |

Tests (8): `test_build_public_smooth_catalog.py`, `test_clean_smooth_catalog.py`, `test_export_catalog_growth_public.py`, `test_export_research_catalog.py`, `test_extend_public_smooth_catalog.py`, `test_mathematical_patterns.py`, `test_package_current_math_release.py`, `test_smooth_isotrivial.py`
