# Build and run the interpolation experiments in several flag configurations.
#   make            build + run the "strict" configuration
#   make run CFG=x  build + run configuration x
#   make matrix     run every configuration and compare the CSVs with strict

SHELL    := /bin/bash
CXX      ?= g++
BASE     := -std=c++20 -Wall -Wextra -Iinclude
LDLIBS   := -lquadmath

CFG      ?= strict
FLAGS_strict   := -O2 -ffp-contract=off
FLAGS_contract := -O2 -march=native -ffp-contract=fast
FLAGS_fast     := -O2 -ffast-math
FLAGS_O0       := -O0 -ffp-contract=off
FLAGS_native   := -O3 -march=native -ffp-contract=off
FLAGS_noinline := -O2 -ffp-contract=off -fno-inline -fno-inline-small-functions
FLAGS_x87      := -O2 -m32 -mfpmath=387 -ffp-contract=off

CONFIGS  := strict native contract fast O0 noinline $(shell echo 'int main(){}' | $(CXX) -m32 -x c++ - -o /dev/null 2>/dev/null && echo x87)

TESTS    := $(basename $(notdir $(wildcard tests/*.cpp)))
HEADERS  := $(wildcard include/*.hpp)
BIN      := build/$(CFG)

.PHONY: all build run matrix clean
all: run

build: $(addprefix $(BIN)/,$(TESTS))

$(BIN)/%: tests/%.cpp $(HEADERS) | $(BIN)
	$(CXX) $(BASE) $(FLAGS_$(CFG)) $< -o $@ $(LDLIBS)

$(BIN):
	mkdir -p $@

# Tests of guaranteed properties exit non-zero on failure; keep going so the
# whole report is produced, and fail at the end.
run: build
	@mkdir -p results/$(CFG); status=0; \
	for t in $(TESTS); do \
	  echo "##### $(CFG) / $$t"; \
	  RESULTS_DIR=results/$(CFG) ./$(BIN)/$$t | tee results/$(CFG)/$$t.txt; \
	  [ $${PIPESTATUS[0]} -eq 0 ] || { echo "!!! $$t failed in $(CFG)"; status=1; }; \
	done; exit $$status

matrix:
	@for c in $(CONFIGS); do $(MAKE) --no-print-directory run CFG=$$c > results-$$c.log 2>&1 \
	  && echo "$$c: all checks passed" || echo "$$c: some checks FAILED (see results-$$c.log)"; done
	@echo; echo "CSV files differing from strict:"
	@for c in $(filter-out strict,$(CONFIGS)); do \
	  for f in results/strict/*.csv; do \
	    cmp -s $$f results/$$c/$$(basename $$f) || echo "  $$c: $$(basename $$f)"; \
	  done; done

clean:
	rm -rf build results results-*.log
