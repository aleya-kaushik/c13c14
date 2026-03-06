# Makefile for SiB4
#
# You will need GNU's version of Make, which is the default on Linux, OS X
# but not necessarily other Unixes.
#
# If you want to change the compiler or optimization without editing the
# makefile, set the COMPILER or OPT environment variables on the command
# line like so:
# user@host:~$ make COMPILER=ifort OPT=opt
#
# Values for COMPILER: gcc*, pgi, ifort
# Values for OPT: debug, opt*
# Values for NC: 3, 4*
# * = default
#
# You can tell the Makefile where to find NetCDF and LAPACK/BLAS libraries
# by setting the NETCDF_ROOT and LAPACK_ROOT environment variables as well:
# user@host:~$ export NETCDF_ROOT=/usr/local # bash
# user@host:~$ setenv NETCDF_ROOT /usr/local # c shell
#
#The includes will be searched for in 
#     $(NETCDF_ROOT)/include and $(LAPACK_ROOT)/include
#The libraries will be searched for in
#     $(NETCDF_ROOT)/lib and $(LAPACK_ROOT)/lib
#
#
#Please note, when linking/running you need to have
#  the libraries loaded into the library path.
#  For bash, this can be done two ways:
#  Option 1) Include the paths in the .bashrc file
#  Option 2) Set the LD_LIBRARY_PATH using the appropriate
#              paths in the files found in the scripts directory
#              (exportgcc, exportifort, exportpgi)
#            In the scripts directory, type:
#              >source ../scripts/exportgcc
#
#=======================================================================
# Set Defaults
# -------------------------------------------------------------------------
ifndef COMPILER
  COMPILER = ifort
endif

ifndef OPT
  OPT = opt
endif

ifndef NC
  NC = 4
endif


# Grab System Information
# -------------------------------------------------------------------------
#SHELL := /bin/sh
#OS = $(shell uname -s)
#PROC = $(shell uname -p)
#NNAME = $(shell uname -n)
#VERSION = $(shell uname -r)

SHELL := /bin/bash
OS = $(uname -s)
PROC = $(uname -p)
NNAME = $(uname -n)
VERSION = $(uname -r)


# Set Compiler Variables
# -------------------------------------------------------------------------
ifeq ($(COMPILER),ifort)
     SUFFIX = intel
     F90 = ifort
endif

ifeq ($(COMPILER),gcc)
   SUFFIX = gcc
   F90 = gfortran
endif

ifeq ($(COMPILER),pgi)
   SUFFIX = pgi
   F90 = pgf90
endif


# Set Library Directories
# -------------------------------------------------------------------------
# BLAS (Basic Linear Algebra System) and LAPACK (Linear Algebra Package)
ifndef LAPACK_ROOT
  ifeq ($(OS),Linux)
       LAPACK_ROOT = $MKLROOT
#       ifeq ($(NNAME),pampas)
#            LAPACK_ROOT = /usr/local/versions/atlas-kathy
#       endif
#       ifeq ($(NNAME),bamboo)
#            LAPACK_ROOT = /usr/local/atlas
#       endif
  endif
endif

# NETCDF and HDF
ifndef NETCDF_ROOT
  ifeq ($(NC),3)
       NETCDF_ROOT = /usr
  else
#       NETCDF_ROOT = /apps/intel-2020.2/netcdf-4.7.4
#       H5ROOT = /apps/intel-2020.2/hdf5-1.10.6
       NETCDF_ROOT = $NETCDF_FORTRAN_ROOT
       H5ROOT = $HDF5_ROOT
  endif

  DIREXIST = $(shell if test -d $(NETCDF_ROOT); then echo yes; else echo no; fi)
  ifeq ($(DIREXIST),no)
#       NETCDF_ROOT = /apps/intel-2020.2/netcdf-4.7.4
#       H5ROOT = /apps/intel-2020.2/hdf5-1.10.6
        NETCDF_ROOT = $NETCDF_FORTRAN_ROOT
        H5ROOT = $HDF5_ROOT
  endif
     
endif


# Set Library Flags
# -------------------------------------------------------------------------
ifeq ($(OS),Darwin)
  ifeq ($(VERSION),10.8.0)
       LALIB = -framework vecLib
  else
       LALIB = -framework Accelerate  
  endif
else
       LALIB = -qmkl=sequential
endif

ifeq ($(NC),3)
#  INCLUDES = -I/apps/intel-2020.2/netcdf-4.7.4/include
#  LIBS = $(LALIB) -L/apps/intel-2020.2/netcdf-4.7.4/lib -lnetcdf
  INCLUDES = -I$(NETCDF_FORTRAN_ROOT)/include
  LIBS = $(LALIB) -L$(NETCDF_FORTRAN_ROOT)/lib -lnetcdf
else
  INCLUDES = -I$(NETCDF_FORTRAN_ROOT)/include -I$(HDF5_ROOT)/include
  LIBS = $(LALIB) -L$(NETCDF_FORTRAN_ROOT)/lib -L$(NETCDF_C_ROOT)/lib \
         -L$(HDF5_ROOT)/lib -lnetcdff -lnetcdf \
          -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lz
#  INCLUDES = -I/apps/intel-2020.2/netcdf-4.7.4/include -I/apps/intel-2020.2/hdf5-1.10.6/include
#  LIBS = $(LALIB) -L/apps/intel-2020.2/netcdf-4.7.4/lib -L/apps/intel-2020.2/hdf5-1.10.6/lib -lnetcdff -lnetcdf \
          -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lz
endif


# Set Compiler Option Flags
# -------------------------------------------------------------------------
ifeq ($(COMPILER),gcc)
  ifeq ($(OPT),opt)
    F90FLAGS = -O2
  else ifeq ($(OPT),debug)
    F90FLAGS = -g -Wall -Wuninitialized -fbounds-check \
         -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=nan
  endif
  ifeq ($(PROC),powerpc)
    F90FLAGS += -fconvert=little-endian
  endif
  F90FLAGS += -fimplicit-none -Wsurprising $(INCLUDES)
  LFLAGS = $(LIBS)
  ifdef PROF
	F90FLAGS += -pg
	LFLAGS += -pg
  endif
endif

ifeq ($(COMPILER),pgi)
  ifeq ($(OPT),opt)
    F90FLAGS = -fast -Mnoframe
  else ifeq ($(OPT),debug)
    F90FLAGS = -g -Mbounds -Ktrap=fp
  endif
  F90FLAGS += -Kieee -Minfo=loop,inline -Minform=inform $(INCLUDES)
  LFLAGS   = -v -Minform=inform $(LIBS)
  ifdef PROF
	F90FLAGS += -pg
	LFLAGS += -pg
  endif
endif

ifeq ($(COMPILER),ifort)
  ifeq ($(OPT),opt)
    F90FLAGS = -O2
  else ifeq ($(OPT),debug)
    F90FLAGS = -fp-stack-check -g -check bounds -fpe0 -traceback -warn interface
  endif
  ifeq ($(PROC),powerpc)
    F90FLAGS += -convert little_endian
  endif
  F90FLAGS += -xHost -ip -no-prec-div -static-intel -implicitnone -heap-arrays 12000000 -warn all $(INCLUDES)
  LFLAGS = $(LIBS)
endif


# Objects (DO NOT EDIT - Unless You Add or Remove Files)
# -------------------------------------------------------------------------
# Variable objects

VAR_OBJS  = \
       kinds.o \
       module_fractsib.o \
       module_oparams.o \
       module_pparams.o \
       module_sibconst.o \
       module_io.o \
       module_isodata.o \
       module_local.o \
       module_param.o \
       module_pftinfo.o \
       module_phosib.o \
       module_poolinfo.o \
       module_sib.o \
       module_sibvs.o \
       module_time.o \
       fcos_solver_Ogee.o \
       cos_soil_production.o

# Scientific objects
SCI_OBJS  = \
        addinc.o \
        balan_carbon.o \
        balan_eh2o.o \
        c13_iso_calc.o \
        cas_update.o \
        cfrax_calc.o \
        cos_calc.o \
        delef.o \
        delhf.o \
        delwf.o \
        dtess_eau.o \
        equipools_control.o \
        flux_rbrd.o \
        flux_update.o \
        flux_vmf.o \
        flux_vrbrd.o \
        hydro_canopy.o \
        hydro_nveg.o \
        hydro_sets.o \
        hydro_setv.o \
        hydro_snow.o \
        hydro_soil.o \
        mapper.o \
        phen_defined.o \
        phen_dynamic.o \
        phen_gss.o \
        phen_update.o \
        phocycalc.o \
        phonveg.o \
        phosib.o \
        phosort.o \
        phostress.o \
        pool_assim.o \
        pool_auto_resp.o \
        pool_auto_tran.o \
        pool_graze.o \
        pool_het_retr.o \
        pool_update.o \
        radabs.o \
        radfac.o \
        radnet.o \
        sibcontrol.o \
        sibmain.o \
        sibslv.o \
        sif_calc.o \
        snow_compact.o \
        snow_combine.o \
        snow_subdivide.o \
        veg_update.o \
        zenith.o 

# Netcdf Objects
NCDF_OBJS = nc_util.o

# SiB Drive objects
DRV_OBJS = \
        co2_module.o \
        cos_module.o \
        diagnostic_create.o \
        diagnostic_save.o \
        diagnostic_write.o \
        driver_init.o \
        driver_interp.o \
        driver_read_global.o \
        driver_read_single.o \
        driver_update.o \
        tm5_init.o \
        tm5_interp.o \
        tm5_read.o \
        tm5_update.o \
	fire_init.o \
	fire_interp.o \
	fire_read_global.o \
	fire_read_single.o \
	fire_update.o \
        grid_init.o \
        local_set.o \
        output_closer.o \
        output_control.o \
        output_init.o \
        read_aero.o \
        read_ciso.o \
        read_namel.o \
        read_outopts.o \
        read_pftinfo.o \
        read_phys.o \
        read_pgdd.o \
        read_pstg.o \
        read_pool.o \
        read_poolinfo.o \
        read_routopts.o \
        read_sibvs.o \
        restart_read.o \
        restart_write.o \
        set_continue_spinup.o \
        sibtype_init.o \
        sibtype_setup.o \
        time_init.o \
        time_manager.o \
        SiBDRV.o

# Compilation (DO NOT EDIT BELOW THIS LINE)
# -------------------------------------------------------------------------
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(F90FLAGS) -c $< -o $@

all: SiB4D

remake: clean all

SiB4D: $(VAR_OBJS) $(SCI_OBJS) $(NCDF_OBJS) $(DRV_OBJS)
	$(F90) $^ -o $@ $(LFLAGS)
	@echo ""
	@date
	@echo Finished building == SiB4D
	@echo ""

clean:
	rm -f SiB4D *.o *.mod *.stb *~ *__genmod*

distclean: clean
	rm -rf SiB4D-* #*#

help:
	@echo ""
	@echo "Supported Rules:"
	@echo "1) [nothing], SiB4D: compile SiB4 model, " \
	"recompiling dependencies as necessary."
	@echo "2) SiB4D: compile SiB4 model only, recompiling dependencies " \
	"as necessary."
	@echo "3) <object name>.o: compile specific object. matches the file" \
	"name."
	@echo "4) remake: delete all objects, and re-compile the sib4 model " \
	"and sibmerge."
	@echo "5) clean: remove all compiler-generated files and executable."
	@echo "6) help: print this message."
	@echo ""

# Dependencies

$(VAR_OBJS):

$(SCI_OBJS): $(VAR_OBJS)

$(NCDF_OBJS): $(VAR_OBJS)

$(DRV_OBJS): $(VAR_OBJS) $(NCDF_OBJS)
