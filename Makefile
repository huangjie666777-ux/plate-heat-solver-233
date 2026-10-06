CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude -Ithird_party/eigen3

LIB_SRCS := src/mesh.cpp src/solver.cpp src/vtk.cpp
LIB_OBJS := $(LIB_SRCS:.cpp=.o)

.PHONY: all example test clean

all: libplate_heat233.a examples/two_material

libplate_heat233.a: $(LIB_OBJS)
	$(AR) rcs $@ $^

examples/two_material: examples/two_material.o libplate_heat233.a
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

example: examples/two_material
	cd examples && ./two_material

test: tests/run_tests
	./tests/run_tests

tests/run_tests: tests/run_tests.o libplate_heat233.a
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(LIB_OBJS) examples/two_material.o examples/two_material \
	    tests/run_tests.o tests/run_tests libplate_heat233.a \
	    examples/two_material.msh examples/two_material.vtk

