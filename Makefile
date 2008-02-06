#===============================================================================
#-------------------------------------------------------------------------------
#
#  MAKEFILE FOR SHEBA
#
#-------------------------------------------------------------------------------
#===============================================================================
#
#  (C) James Wookey, December 2003 - February 2008
#  Department of Earth Sciences, University of Bristol
#  Wills Memorial Building, Queen's Road, Bristol, BR8 1RJ, UK
#  j.wookey@bristol.ac.uk
#
#-------------------------------------------------------------------------------
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
#
#-------------------------------------------------------------------------------
#   CVS: $Revision: 1.10 $ $Date: 2008/02/06 09:51:48 $
#-------------------------------------------------------------------------------
#
#   SHEBA requires a FORTRAN 90 compiler. Compilers known to work are:    
#      Linux : ifc/ifort, g95, gfortran
#      Mac (PowerPC) : g95, gfortran, XLF
#      Mac (i386) : gfortran, ifort (maybe)
#      Sun (sparc) : Solaris f90
#
#-------------------------------------------------------------------------------
# Set compiler and options here:

# -- flags which work for ifort (version>8) on MacOSX
FC= ifort 
OPT77 = -w95 -cm 
OPT90 =
OPT =  -axT -Vaxlib -assume byterecl    

# -- flags which work for f90 on Solaris
#FC = f90
#OPT90 =
#OPT77 =
#OPT = 

# -- flags which work for gfortran on Mac (and probably Linux)
#FC = gfortran -O2
#OPT90 =
#OPT77 = -w
#OPT = 

#-------------------------------------------------------------------------------
## Uncomment for sac/sac2000 (binaries sac files are native endian)
#F90SAC_FLAGS = -DDISABLE_C_IO
#-------------------------------------------------------------------------------
## Uncomment for MacSAC (binary sac files are always big-endian)
F90SAC_FLAGS = -DDISABLE_C_IO -DFORCE_BIGENDIAN_SACFILES 
#-------------------------------------------------------------------------------
# Note: f90sac C routines are disabled for sheba; we don't need them.

EXECDIR=/usr/local/sac/macros

#-------------------------------------------------------------------------------
# No editing *should* be required below here ...
#-------------------------------------------------------------------------------
#
#	Code Objects
#
MODULES = sheba_config.o event_info.o array_sizes.o
SUBROUTINES = sheba.o misc.o input.o desplit.o output.o teanby.o \
                  rumpker.o scs.o cluster.o split.o\
                  traceops.o particle.o hilbert.o scsslow.o\
                  calcsnr.o reorient.o

F90SAC = f90sac_distrib.o 
#
#	Executable info
#
# 
#

all:$(EXECDIR)/sheba \
      $(EXECDIR)/sheba_plot_errclu.gmt \
      $(EXECDIR)/sheba_plot_error.gmt \
      $(EXECDIR)/sheba_combine_plots.csh \
      $(EXECDIR)/cleansheba \
      $(EXECDIR)/split_sheba\
      $(EXECDIR)/sheba_stack 
      
#
#     SHEBA EXECUTABLE
#
$(EXECDIR)/sheba:${F90SAC} ${MODULES} sheba_main.o ${SUBROUTINES}
	$(FC) $(OPT) -o $(EXECDIR)/sheba ${F90SAC} ${MODULES} sheba_main.o ${SUBROUTINES}
#
$(EXECDIR)/sheba_stack:${F90SAC} ${MODULES} sheba_stack.o ${SUBROUTINES}
	$(FC) $(OPT) -o $(EXECDIR)/sheba_stack ${F90SAC} ${MODULES} sheba_stack.o ${SUBROUTINES}

#	F90SAC requires special options to compile ...
f90sac_distrib.o: f90sac_distrib.F90
	$(FC) $(OPT) $(OPT90) $(F90SAC_FLAGS) -c $*.F90 

#
#     GMT PLOTTING SCRIPTS + OTHER SHELL SCRIPTS
#
$(EXECDIR)/sheba_plot_errclu.gmt:sheba_plot_errclu.gmt
	chmod +x sheba_plot_errclu.gmt; \cp sheba_plot_errclu.gmt $(EXECDIR)

$(EXECDIR)/sheba_plot_error.gmt:sheba_plot_error.gmt
	chmod +x sheba_plot_error.gmt; \cp sheba_plot_error.gmt $(EXECDIR)

$(EXECDIR)/sheba_combine_plots.csh:sheba_combine_plots.csh
	chmod +x sheba_combine_plots.csh; \cp sheba_combine_plots.csh $(EXECDIR)

$(EXECDIR)/cleansheba:cleansheba
	chmod +x cleansheba; \cp cleansheba $(EXECDIR)

#
#     SAC Macro
#
$(EXECDIR)/split_sheba:split_sheba
	cp -f split_sheba $(EXECDIR)

distrib:
	rm -rf ../SHEBA_distrib
	mkdir ../SHEBA_distrib
	cp -f *.f90 ../SHEBA_distrib
	cp -f *.f ../SHEBA_distrib
	cp -f *.gmt ../SHEBA_distrib
	cp -f *.csh ../SHEBA_distrib
	cp -f split_sheba ../SHEBA_distrib
	cp -f Makefile ../SHEBA_distrib
	cp -f INSTALL ../SHEBA_distrib
	cp -f cleansheba ../SHEBA_distrib	
clean:
	rm -f *.o *.mod

# 
# 	Compile Instuctions
#
%.o: %.F90
	$(FC) $(OPT) $(OPT90) -c $*.F90 
%.o: %.f90
	$(FC) $(OPT) $(OPT90) -c $*.f90 
%.o: %.f
	$(FC) $(OPT) $(OPT77) -c $*.f

#
#	Dependencies
#
$(SUBROUTINES)    :	$(MODULES)		
