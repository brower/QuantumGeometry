CXX=g++
CXXFLAGS=-std=c++17 -O3 -Wall -Wno-unknown-warning-option \
	-Wno-deprecated-declarations -Wno-sign-compare
# CXXFLAGS=-std=c++14 -O3 -Wall -Xpreprocessor -fopenmp

# directory where executables can access symmetry group data files
GRP_DIR=$(shell pwd)/../grp
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""

SRC=$(wildcard *.c*)
BIN_DIR=../bin
BIN=$(basename $(SRC))
INC_DIR=../include
INC=$(wildcard $(INC_DIR)/*.h)

.PHONY: clean

all: $(BIN)

%: %.cc $(INC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(INC_DIR) $< -o $(BIN_DIR)/$@

clean:
	$(RM) -r $(addprefix $(BIN_DIR)/, $(BIN))
