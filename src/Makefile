

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

