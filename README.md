#modules 
gcc/9.2.0 openmpi/3.1.6

# C++ Reservoir Model for the Columbia River Basin (CRB)

This repository contains a C++ reservoir model designed to optimize operations for the Columbia River Basin (CRB). The model includes historical reservoir rules and provides an interface between the reservoir model and an optimization solver.

## Project Structure

The project consists of three main files:

- **HYSSR_CRB_model.cpp**: Contains the core functionality for optimizing reservoir operations.
- **HYSSR_CRB_model_hyssr.cpp**: Implements historical reservoir rules.
- **runCRBparallelJared.cpp**: Provides the interface between the reservoir model and the optimization solver.

## Compilation

The project uses a Makefile to compile the necessary files. The Makefile assumes the use of `mpicxx` for parallel compilation. 

### Makefile

```Makefile
####### Compiler, tools and options

# CC            = gcc
CXX           = mpicxx
CXXFLAGS      = -c -O3 -Wall 

####### Compile
all: CRB_FCRPS

CRB_FCRPS: runCRBparallelJared.o borgmm.o  mt19937ar.o
	$(CXX) runCRBparallelJared.o borgmm.o mt19937ar.o -o CRB_FCRPS

runCRBparallelJared.o: runCRBparallelJared.cpp HYSSR_CRB_model.cpp \
	../borg_jared/borgmm.h
	$(CXX) $(CXXFLAGS) runCRBparallelJared.cpp

borgmm.o: ../borg_jared/borgmm.c ../borg_jared/borgmm.h \
	../borg_jared/borg.h
	$(CXX) $(CXXFLAGS) ../borg_jared/borgmm.c

mt19937ar.o: ../borg_jared/mt19937ar.c ../borg_jared/mt19937ar.h
	$(CXX) $(CXXFLAGS) ../borg_jared/mt19937ar.c

clean:
	rm -rf *.o 
	rm CRB_FCRPS
