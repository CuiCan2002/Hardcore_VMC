CXX ?= g++
CPPFLAGS += -Iinclude
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

ifeq ($(OS),Windows_NT)
EXE := .exe
else
EXE :=
endif

.PHONY: all test clean
all: hardcore_vmc$(EXE)

MPICXX ?= mpicxx
.PHONY: mpi
mpi: hardcore_vmc_mpi$(EXE)
hardcore_vmc_mpi$(EXE): src/main.cpp src/vmc.cpp include/vmc.hpp
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -DVMC_MPI src/main.cpp src/vmc.cpp -o $@ $(LDFLAGS)

hardcore_vmc$(EXE): src/main.cpp src/vmc.cpp include/vmc.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/main.cpp src/vmc.cpp -o $@ $(LDFLAGS)

vmc_tests$(EXE): tests/test_vmc.cpp src/vmc.cpp include/vmc.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/test_vmc.cpp src/vmc.cpp -o $@ $(LDFLAGS)

test: vmc_tests$(EXE)
	./vmc_tests$(EXE)

clean:
ifeq ($(OS),Windows_NT)
	-if exist hardcore_vmc.exe del hardcore_vmc.exe
	-if exist hardcore_vmc_mpi.exe del hardcore_vmc_mpi.exe
	-if exist vmc_tests.exe del vmc_tests.exe
else
	rm -f hardcore_vmc hardcore_vmc_mpi vmc_tests
endif
