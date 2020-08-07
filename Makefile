#----------------------------------------------------------------------!
#                                                                      !
#                Makefile for HFBTHO                                   !
#                                                                      !
# Authors: N. Schunck, LLNL                                            !
#                                                                      !
# This compiles HFBTHO (version > 139)                                 !
#                                                                      !
#                                   First included in release v135     !
#                                                                      !
#----------------------------------------------------------------------!

#======================================================================#
#                                                                      #
#              HFBTHO - Makefile (version >= 135)                      #
#                                                                      #
#  NOTA BENE: When running this Makefile, make sure:                   #
#             - all files are in UNIX format (this includes            #
#               all sources .f and .f90, but also this Makefile)       #
#             - all indentations after 'Beginning of the action' are   #
#               proper tabulations, not spaces...                      #
#                                                                      #
#  Pre-defined compilers are Portland Fortran (PGI), Pathscale         #
#  (PATHSCALE), Lahey (LAHEY), Intel Fortran compiler (IFORT), GNU     #
#  fortran (GFORTRAN) and CRAY fortran (CRAY).                         #
#                                                                      #
#  Makefile Options                                                    #
#  ================                                                    #
#                                                                      #
#   - COMPILER .....: Presets compiler's options. Choices are:         #
#                     PATHSCALE, LAHEY, PGI, IFORT, GFORTRAN, CRAY     #
#   - FORTRAN ......: Set compiler's name                              #
#   - DEBUG ........: Turns on all optimizations on (FALSE) or off     #
#                     (TRUE). In the latter case, debugging options    #
#                     are set.                                         #
#   - ACML_LIBRARY .: If TRUE, uses an external library for BLAS and   #
#                     LAPACK (you must set the links to the library    #
#                     under keyword ACML below), if FALSE, the library #
#                     is built.                                        #
#   - STATIC_LIB ...: TRUE or FALSE, allows static linking             #
#   - hide_openmp ..: If 1, no OpenMP; if 0, OpenMP pragmas are active #
#                                                                      #
#======================================================================#

SHELL := /bin/bash

# Names
VERSION        = 200d
HFBTHO_EXE     = main
HFBTHO_SOURCE  = main_v$(VERSION).f90
HFBTHO_OBJ     = main_v$(VERSION).o

# Makefile Options
COMPILER      = GFORTRAN
FORTRAN       = gfortran
DEBUG         = FALSE
ACML_LIBRARY  = FALSE
STATIC_LIB    = FALSE

# Preprocessor options
hide_openmp   = 1

# BLAS and LACPACK libraries
ifeq ($(ACML_LIBRARY),TRUE)
     LAPACK_OBJ =
     BLAS_OBJ   =
     ACML       = -L$(HOME)/lib/lapack-3.4.1 -llapack_LINUX -lblas_LINUX
else
     LAPACK_OBJ = lapack.o
     BLAS_OBJ   = blas.o
     ACML       =
endif

# Module names
HFBTHO_SOLVER_MOD = hfbtho.mod
HFBTHO_SOLVER_OBJ = hfbtho_v$(VERSION).o

# Defining compiler options for: GNU FORTRAN COMPILER (gfortran)
ifeq ($(COMPILER),GFORTRAN)

      FORMAT_F77 = -ffixed-form
      FORMAT_F90 = -ffree-form -ffree-line-length-none
      STATIC     =

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = -static
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = -cpp -Dhide_openmp=$(hide_openmp) $(STATIC) -O3
      else
            OPTIONS_FC = -cpp -Dhide_openmp=$(hide_openmp) \
                              $(STATIC) -Warray-bounds -Wall -Wunderflow -Warray-temporaries \
                              -Wcharacter-truncation -Wtabs -Wintrinsic-shadow -Walign-commons \
                              -frange-check -Wconversion
      endif

      ifeq ($(hide_openmp),0)
            OPTIONS = $(OPTIONS_FC) -fopenmp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: IFORT FORTRAN COMPILER (ifort)
ifeq ($(COMPILER),IFORT)

      FORMAT_F77 = -fixed -80
      FORMAT_F90 = -free -extend_source
      STATIC     =

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = -static
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = -cpp -Dhide_openmp=$(hide_openmp) $(STATIC) -O3
      else
            OPTIONS_FC = -cpp -Dhide_openmp=$(hide_openmp) $(STATIC) -check all
      endif

      ifeq ($(hide_openmp),0)
            OPTIONS = $(OPTIONS_FC) -openmp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: PORTLAND FORTRAN COMPILER (pgf90)
ifeq ($(COMPILER),PGI)

      FORMAT_F77 = -fpic -m64 -tp=shanghai-64 -Mfixed
      FORMAT_F90 = -fpic -m64 -tp=shanghai-64
      STATIC     =

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = -Bstatic
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = -Mpreprocess -Dhide_openmp=$(hide_openmp) $(STATIC) -fast
      else
            OPTIONS_FC = -Mpreprocess -Dhide_openmp=$(hide_openmp) \
                                      $(STATIC) -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Minfo=all
      endif

      ifeq ($(hide_openmp),0)
            OPTIONS = $(OPTIONS_FC) -mp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: PATHSCALE COMPILER (pathf90)
ifeq ($(COMPILER),PATHSCALE)

      FORMAT_F77 = -fixedform
      FORMAT_F90 =
      STATIC     =

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = -static
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = -cpp -Dhide_openmp=$(hide_openmp) $(STATIC) -woff -O3
      else
            OPTIONS_FC = -cpp -Dhide_openmp=$(hide_openmp) -fullwarn -C
      endif

      ifeq ($(hide_openmp),0)
            OPTIONS = $(OPTIONS_FC) -mp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: LAHEY COMPILER (lf95)
ifeq ($(COMPILER),LAHEY)

      FORMAT_F77 = --fix
      FORMAT_F90 =
      STATIC     =

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = --staticlink
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = -Cpp -Dhide_openmp=$(hide_openmp) $(STATIC) -O3
      else
            OPTIONS_FC = -Cpp -Dhide_openmp=$(hide_openmp) \
                              $(STATIC) --chk --chkglobal
      endif

      ifeq ($(hide_openmp),0)
            OPTIONS = $(OPTIONS_FC) --openmp
      else
            OPTIONS = $(OPTIONS_FC)
      endif

endif

# Defining compiler options for: CRAY COMPILER
ifeq ($(COMPILER),CRAY)

      FORMAT_F77 = -f fixed
      FORMAT_F90 = -f free
      STATIC     =

      ifeq ($(STATIC_LIB),TRUE)
              STATIC = --staticlink
      endif

      ifeq ($(DEBUG),FALSE)
            OPTIONS_FC = -e Z -Dhide_openmp=$(hide_openmp) $(STATIC) -O3
      else
            OPTIONS_FC = -e Z -Dhide_openmp=$(hide_openmp) $(STATIC) -e c -e D
      endif

      ifeq ($(USE_OPENMP),1)
            OPTIONS = $(OPTIONS_FC)
      else
            OPTIONS = $(OPTIONS_FC) -h noomp
      endif

endif

#=========================#
# Beginning of the action #
#=========================#

all: $(HFBTHO_EXE)

# Libraries
$(LAPACK_OBJ) : lapack.f
	$(FORTRAN) $(FORMAT_F77) $(OPTIONS) -c $<

$(BLAS_OBJ) : blas.f
	$(FORTRAN) $(FORMAT_F77) $(OPTIONS) -c $<

LIBRARIES = $(BLAS_OBJ) $(LAPACK_OBJ)

# HFBTHO core
$(HFBTHO_SOLVER_MOD) : hfbtho_v$(VERSION).f90
	$(FORTRAN) $(FORMAT_F90) $(OPTIONS) -c $<

# HFBTHO core
$(HFBTHO_OBJ) : $(HFBTHO_SOURCE) $(HFBTHO_SOLVER_MOD)
	$(FORTRAN) $(FORMAT_F90) $(OPTIONS) -c $<

$(HFBTHO_EXE) : $(HFBTHO_OBJ) $(HFBTHO_SOLVER_OBJ) $(LIBRARIES)
	$(FORTRAN) $(FORMAT_F90) $(OPTIONS) -o $@ $^ $(ACML)


clean ::
	-rm -f *.o *.oo *.ipo *.mod


