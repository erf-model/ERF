#---------------------------
# Makefile for AMReX project
#---------------------------

CXX = mpicxx
CXXFLAGS = -O2 -std=c++17 -fopenmp -Wall

# AMReX paths
AMREX_INC = /pscratch/sd/n/nataraj2/AIModCon/GitHub/amrex_repos/amrex/include
AMREX_LIB = /pscratch/sd/n/nataraj2/AIModCon/GitHub/amrex_repos/amrex/lib

HDRS = ReadPlotFile.H

LIBS = -lamrex -lpthread -lgfortran -ldl -lm -lstdc++

# Executable
EXE = out

SRCS = main.cpp ReadPlotFile.cpp
OBJS = $(SRCS:.cpp=.o)


# Default target
all: $(EXE)


# Link executable
$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(AMREX_LIB) $(LIBS)


# Compile object files
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -I$(AMREX_INC) -c $< -o $@


# Clean
clean:
	rm -f *.o $(EXE)


.PHONY: all clean
