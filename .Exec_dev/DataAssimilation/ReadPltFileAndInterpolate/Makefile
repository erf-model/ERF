

#---------------------------
# Makefile for AMReX project
#---------------------------
CXX = mpicxx
CXXFLAGS = -O2 -std=c++17 -fopenmp -Wall

# Replace with your AMReX paths
AMREX_INC = /pscratch/sd/n/nataraj2/AIModCon/GitHub/amrex_repos/amrex/include
AMREX_LIB = /pscratch/sd/n/nataraj2/AIModCon/GitHub/amrex_repos/amrex/lib

HDRS = ReadPlotFile.H
LIBS = -lamrex -lpthread -lgfortran -ldl -lm -lstdc++

# Default: do nothing if 'make' is called without target
all:
	@echo "Use 'make serial' or 'make parallel'"

# -------------------------
# Serial build
# -------------------------
SRCS_SERIAL = main_serial.cpp ReadPlotFile.cpp
OBJS_SERIAL = $(SRCS_SERIAL:.cpp=.o)
EXE_SERIAL  = out_serial

serial: $(EXE_SERIAL)

$(EXE_SERIAL): $(OBJS_SERIAL)
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(AMREX_LIB) $(LIBS)

# -------------------------
# Parallel build
# -------------------------
SRCS_PARALLEL = main_parallel.cpp ReadPlotFile.cpp
OBJS_PARALLEL = $(SRCS_PARALLEL:.cpp=.o)
EXE_PARALLEL  = out_parallel

parallel: $(EXE_PARALLEL)

$(EXE_PARALLEL): $(OBJS_PARALLEL)
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(AMREX_LIB) $(LIBS)

# -------------------------
# Compile rules for object files
# -------------------------
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -I$(AMREX_INC) -c $< -o $@

# -------------------------
# Clean
# -------------------------
clean:
	rm -f *.o out_serial out_parallel

.PHONY: all serial parallel clean


