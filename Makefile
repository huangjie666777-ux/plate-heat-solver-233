CXX      ?= g++
EIGEN    ?= third_party/eigen
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
CPPFLAGS += -Iinclude -I$(EIGEN)

BUILD    := build
LIB      := $(BUILD)/libplate_heat233.a
LIBOBJ   := $(BUILD)/mesh.o $(BUILD)/solve.o $(BUILD)/vtk.o

.PHONY: all test run-example clean

all: $(LIB) $(BUILD)/dual_material $(BUILD)/selftest

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.cpp include/plate_heat233/plate_heat233.hpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIB): $(LIBOBJ)
	$(AR) rcs $@ $^

$(BUILD)/dual_material: examples/dual_material.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L$(BUILD) -lplate_heat233 -o $@

$(BUILD)/selftest: tests/selftest.cpp $(LIB)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L$(BUILD) -lplate_heat233 -o $@

test: $(BUILD)/selftest
	./$(BUILD)/selftest

run-example: $(BUILD)/dual_material
	cd $(BUILD) && ./dual_material

clean:
	rm -rf $(BUILD)
