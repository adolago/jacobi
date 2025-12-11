


EXE=JACOBI2D

.DEFAULT_GOAL=JACOBI2D

OBJECTS=$(addprefix build/,main.o functions.o)

build/functions.o : functions.h

CXX=g++
CXX_CFLAGS_MANDATORY=-std=c++14

CXX_CFLAGS_OPTIONAL=-O3


CXX_CFLAGS=$(CXX_CFLAGS_MANDATORY) $(CXX_CFLAGS_OPTIONAL)

CXX_LFLAGS=


build/%.o:  ./%.cpp
	[ -d $(@D) ] || mkdir $(@D) && true
	$(CXX) $(CXX_CFLAGS) -c $< -o $@

$(EXE):  $(OBJECTS)
	$(CXX) $(CXX_LFLAGS) $^ -o $@
