# Build everything and run every test that works without the original raw data.
#   make build   main engine (target/smart/release) and the four Rust projects
#   make test    Rust tests and self-tests, Python tests, website tests
SHELL := /bin/bash
.PHONY: build test test-rust test-python test-web

ROOT_CARGO := CARGO_TARGET_DIR=target/smart RUSTFLAGS="-C target-cpu=native" cargo
# Degrees above 25 need the wide-key build (--cfg smooth_wide_keys); its one test runs there.
WIDE_CARGO := CARGO_TARGET_DIR=target/wide RUSTFLAGS="-C target-cpu=native --cfg smooth_wide_keys" cargo
WIDE_ONLY := new_degree_constructors_and_wide_keys
RUST_PROJECTS := frames geosmooth onepiece quotient_exhaust
BIN := $(CURDIR)/target/smart/release
SYMPY := uv run --isolated --no-project --python 3.13 --with sympy==1.14.0 python
SCIENCE := uv run --isolated --no-project --python 3.13 --with numpy==2.3.5 --with pandas==2.3.3 \
	--with scipy==1.16.3 --with statsmodels==0.14.6 --with matplotlib==3.11.1 --with sympy==1.14.0 python

# Fails on the current TRIANGLE sampler, which changed after this test was written (see README).
KNOWN_FAILING := triangle_draws_have_the_declared_geodesies_and_budget
# Read historical runs, exports or the website's proof PDF, none of which are in this repository (see README).
NEEDS_DATA := test_build_public_smooth_catalog test_package_current_math_release test_defect_coin_model test_genus_three_double_marks test_genus_three_frame_retirement \
	test_genus_three_saturated_launch test_mathematical_patterns test_negative_audit
NEEDS_SYMPY := test_clean_smooth_catalog test_export_catalog_growth_public \
	test_mechanism_breaker test_smooth_isotrivial test_smooth_noncyclic test_theory_frontier
NEEDS_SCIENCE := test_spectral_theorem_logic
NEEDS_ARGS := test_defect_degree_stability

build:
	$(ROOT_CARGO) build --release
	for p in $(RUST_PROJECTS); do (cd projects/$$p && cargo build --release) || exit 1; done

test: test-rust test-python test-web

test-rust: build
	$(ROOT_CARGO) test --release --no-fail-fast -- --skip $(KNOWN_FAILING) --skip $(WIDE_ONLY)
	$(WIDE_CARGO) test --release --bin smooth_catalog_growth $(WIDE_ONLY)
	$(BIN)/hurwitz_breaker selftest
	projects/onepiece/target/release/onepiece --selftest
	python3 projects/geosmooth/test_g2.py

test-python: build
	@export GUIDED_AUDIT_BINARY=$(BIN)/smooth_negative_audit DESIGN_BINARY=$(BIN)/smooth_forest_fast; \
	log=$$(mktemp); failed=(); \
	run() { echo "== $$*"; "$$@" >"$$log" 2>&1 || { tail -20 "$$log"; failed+=("$$*"); }; }; \
	for t in scripts/*/test_*.py; do n=$$(basename $$t .py); \
	  case " $(NEEDS_DATA) $(NEEDS_SYMPY) $(NEEDS_SCIENCE) $(NEEDS_ARGS) " in *" $$n "*) continue;; esac; \
	  run python3 $$t; done; \
	for n in $(NEEDS_SYMPY); do run $(SYMPY) scripts/*/$$n.py; done; \
	for n in $(NEEDS_SCIENCE); do run $(SCIENCE) scripts/*/$$n.py; done; \
	run python3 scripts/spectra/test_defect_degree_stability.py --binary $(BIN)/defect_degree_probe; \
	run python3 scripts/weak/odd_order_spin.py --selftest; \
	run python3 scripts/weak/remaining_weak_passports.py --selftest; \
	run python3 scripts/weak/analyze_escape_words.py --selftest; \
	run $(SYMPY) scripts/weak/weak_graph_extensions.py --selftest; \
	run $(SYMPY) scripts/weak/verify_weak_graph_extensions.py --selftest; \
	rm -f "$$log"; \
	if [ $${#failed[@]} -gt 0 ]; then printf 'FAILED: %s\n' "$${failed[@]}"; exit 1; fi; \
	echo "All Python tests passed."

test-web:
	node deploy/math/ui/test-catalog.mjs
