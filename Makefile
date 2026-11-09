CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
PYTHON ?= python3

.PHONY: verify figures reproduce
verify: build/experiments
	./build/experiments verify

build/experiments: code/experiments.cpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@

figures:
	$(PYTHON) code/figures.py

reproduce: build/experiments
	mkdir -p build/reproduced results
	./build/experiments verify > build/reproduced/verification.txt
	./build/experiments one > build/reproduced/single_link.csv
	./build/experiments finite > build/reproduced/finite_budgets.csv
	./build/experiments limit > build/reproduced/limiting_constants.csv
	cp build/reproduced/verification.txt build/reproduced/*.csv results/
	$(MAKE) figures
