#---------------------------
# Makefile for AMReX project
#---------------------------

CXX = mpicxx
CXXFLAGS = -O2 -std=c++17 -fopenmp -Wall

# Replace with your AMReX paths
AMREX_INC = /pscratch/sd/n/nataraj2/AIModCon/GitHub/amrex_repos/amrex/include
AMREX_LIB = /pscratch/sd/n/nataraj2/AIModCon/GitHub/amrex_repos/amrex/lib

SRCS = main.cpp ReadPlotFile.cpp
HDRS = ReadPlotFile.H

OBJS = $(SRCS:.cpp=.o)
EXE  = out

# Add math, dl, pthread as required by AMReX
LIBS = -lamrex -lpthread -lgfortran -ldl -lm

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(AMREX_LIB) $(LIBS)

# Compile rule with header dependency
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -I$(AMREX_INC) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXE)

.PHONY: all clean
