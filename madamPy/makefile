.PHONY: default
default: build test

PYTHON = /Users/dekelsch/homebrew/Cellar/python\@3.7/3.7.12/bin/python3.7
#PYTHON = python3
NUMPY_INCLUDE  = ${shell ${PYTHON} -c 'import numpy; print( numpy.get_include() )'}
MPI4PY_INCLUDE = ${shell ${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'}
PYTHON_INCLUDE = /Users/dekelsch/homebrew/Cellar/python\@3.7/3.7.12/Frameworks/Python.framework/Versions/3.7/include/python3.7m/
PYTHON_LIB = /Users/dekelsch/homebrew/Cellar/python\@3.7/3.7.12/Frameworks/Python.framework/Versions/3.7/lib



include ../module.mk

MPICC = mpic++
BREW = $(HOME)/homebrew/Cellar
CXXFLAGS += -fPIC
LDFLAGS += -L${PYTHON_LIB} -shared -ldl -framework CoreFoundation -lpython3.7m 
SOURCES=energy.cpp ../src/*.cpp
.PHONY: build 
build: PartiClass.so
PartiClass.so: PyPartiClass.cpp
	${MPICC} ${CXXFLAGS} -I${NUMPY_INCLUDE} -I${MPI4PY_INCLUDE} -I${PYTHON_INCLUDE} -o $@ $< $(SOURCES) ${LDFLAGS} ${LDLIBS}

MPIEXEC = mpiexec
NP_FLAG = -n
NP = 4
.PHONY: test
test: build
	${MPIEXEC} -f hostfile ${NP_FLAG} ${NP} ${PYTHON} run_partiClass.py


.PHONY: clean
clean:
	${RM} *.png*
	${RM} *.dat*
	${RM} PartiClass.so
