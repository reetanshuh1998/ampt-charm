# Open-charm AMPT framework: one-command entry points (see README.md)
.PHONY: all deps build test-results test-repro test-quick test clean
all: deps build test

deps:           ## download + verify AMPT, LHAPDF 6.5.1, CTEQ6.1M, EPS09s; build LHAPDF (deps/)
	scripts/setup_deps.sh
build:          ## patch the official AMPT and compile bin/ampt
	scripts/build.sh
test-results:   ## regenerate all published tables and figures from the shipped analysis outputs (~30 s)
	tests/test_results.sh
test-repro:     ## regenerate events of the paper production and compare byte for byte (~2 min)
	tests/test_exact_reproduction.sh
test-quick:     ## small new production + analysis + independent re-analysis + comparison with published (~5-10 min)
	tests/test_quick_production.sh
test: test-results test-repro test-quick
clean:
	rm -rf build bin quick prod
