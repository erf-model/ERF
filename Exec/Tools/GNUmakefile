ERF_HOME   ?= ../..
AMREX_HOME ?= $(ERF_HOME)/Submodules/AMReX

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE

BL_NO_FORT = TRUE

USE_NETCDF = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

DEFINES += -DERF_USE_NETCDF

EBASE = PlotfileToNetCDF

#Turn on NetCDF macro define
has_netcdf_mpi := $(shell pkg-config --cflags netcdf-mpi > /dev/null 2>&1; echo $$?)
ifeq ($(has_netcdf_mpi),0)
 includes  += $(shell pkg-config --cflags netcdf-mpi)
 LIBRARIES += $(shell pkg-config --libs netcdf-mpi)
else
 includes  += $(shell pkg-config --cflags netcdf)
 LIBRARIES += $(shell pkg-config --libs netcdf)
endif

CEXE_sources += ERF_PlotfileToNetCDF.cpp
CEXE_sources += ERF_NCPlotFile.cpp
CEXE_sources += ERF_NCInterface.cpp
CEXE_headers += ERF_NCInterface.H

INCLUDE_LOCATIONS += $(ERF_HOME)/Source/Tools
VPATH_LOCATIONS   += $(ERF_HOME)/Source/Tools

INCLUDE_LOCATIONS += $(ERF_HOME)/Source/IO
VPATH_LOCATIONS   += $(ERF_HOME)/Source/IO

# ERF_NCPlotFile.{H,cpp} need ERF_Constants.H and ERF_DataStruct.H (and the
# structs the latter pulls in)
INCLUDE_LOCATIONS += $(ERF_HOME)/Source
VPATH_LOCATIONS   += $(ERF_HOME)/Source

INCLUDE_LOCATIONS += $(ERF_HOME)/Source/DataStructs
VPATH_LOCATIONS   += $(ERF_HOME)/Source/DataStructs

INCLUDE_LOCATIONS += $(ERF_HOME)/Source/PBL
VPATH_LOCATIONS   += $(ERF_HOME)/Source/PBL

INCLUDE_LOCATIONS += $(ERF_HOME)/Source/Utils
VPATH_LOCATIONS   += $(ERF_HOME)/Source/Utils

include $(AMREX_HOME)/Src/Base/Make.package
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(EBASE).ex
