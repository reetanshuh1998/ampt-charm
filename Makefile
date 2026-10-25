# Open-charm AMPT framework: one-command entry points (see README.md)
.PHONY: all deps build eps-tool test-results test-repro test-quick test clean
all: deps build test

deps:           ## download + verify AMPT, LHAPDF 6.5.1, CTEQ6.1M, EPS09s; build LHAPDF (deps/)
	scripts/setup_deps.sh
build:          ## patch the official AMPT and compile bin/ampt
	scripts/build.sh
eps-tool:       ## offline EPS09s error-set evaluator bin/eps09s_rg (README Sec. 2.3; needs make deps)
	mkdir -p bin && g++ -O2 -I deps/eps09s analysis/eps09s/eps09s_rg.cpp -o bin/eps09s_rg
test-results:   ## regenerate all published tables and figures from the shipped analysis outputs (~30 s)
	tests/test_results.sh
test-repro:     ## regenerate events of the paper production and compare byte for byte (~2 min)
	tests/test_exact_reproduction.sh
test-quick:     ## small new production + analysis + independent re-analysis + comparison with published (~10-20 min)
	tests/test_quick_production.sh
test: test-results test-repro test-quick
clean:
	rm -rf build bin quick prod regenerated
