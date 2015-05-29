#===============================================================================
#-------------------------------------------------------------------------------
#
#  MAKEFILE FOR SHEBA
#
#-------------------------------------------------------------------------------
#===============================================================================
#
#  (C) James Wookey, December 2003 - February 2011
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
#
#   SHEBA requires a FORTRAN 90 compiler. Compilers known to work are:    
#      Linux : ifc/ifort, g95, gfortran
#      Mac (PowerPC) : g95, gfortran, XLF
#      Mac (i386) : gfortran, ifort (maybe)
#      Sun (sparc) : Solaris f90
#


#===============================================================================
# Path options
#===============================================================================

# Delivery path
MACRODIR=/usr/local/sac/macros
EXECDIR=/usr/local/sac/macros

#===============================================================================
# Compiler and options:
#===============================================================================

## -- flags which work for gfortran on Mac (and probably Linux)
FC = gfortran 
# Production
FFLAGS = -O2
# Debug
# FFLAGS = -fcheck=bounds -C -g
# additional flags for .f/.F files 
F77FLAGS = -w
# additional flags for .f90/.F90 files 
F90FLAGS =

#===============================================================================
# SAC compatibility options
#===============================================================================

## Uncomment for sac/sac2000 (binaries sac files are native endian)
#F90SAC_FLAGS = 

## Uncomment for MacSAC (binary sac files are always big-endian)
F90SAC_FLAGS = -DFORCE_BIGENDIAN_SACFILES 

#===============================================================================
# GMT Compatibility options
#===============================================================================

## gmt prefix (required by some gmt installations)
GMT_PREFIX = gmt4 

## For traditional gmt installations
#GMT_PREFIX = 

#===============================================================================
#===============================================================================
# No editing *should* be required below here ...
#===============================================================================
#===============================================================================
#
#	Code Objects
#
MODULES = sheba_config.o event_info.o array_sizes.o
SUBROUTINES = sheba_core.o misc.o input.o desplit.o output.o teanby.o \
                  rumpker.o cluster.o split.o\
                  traceops.o particle.o ndf.o\
                  calcsnr.o reorient.o crosscorr.o

F90SAC = f90sac_distrib.o 
#
#	Executable info
#
# 
#

all:$(EXECDIR)/sheba_exec \
      $(EXECDIR)/sheba_plot_errclu.gmt \
      $(EXECDIR)/sheba_plot_stackerr.gmt \
      $(EXECDIR)/sheba_combine_plots.csh \
      $(EXECDIR)/cleansheba \
      $(MACRODIR)/split_sheba\
      $(MACRODIR)/sheba\
      $(EXECDIR)/sheba_stack\
		$(EXECDIR)/stack_wgtcalc
      
#
#     SHEBA EXECUTABLE
#
$(EXECDIR)/sheba_exec:${F90SAC} ${MODULES} sheba_main.o ${SUBROUTINES}
	$(FC) $(FFLAGS) -o $(EXECDIR)/sheba_exec ${F90SAC} ${MODULES} sheba_main.o ${SUBROUTINES}
#
$(EXECDIR)/sheba_stack:${F90SAC} ${MODULES} sheba_stack.o ${SUBROUTINES}
	$(FC) $(FFLAGS) -o $(EXECDIR)/sheba_stack ${F90SAC} ${MODULES} sheba_stack.o ${SUBROUTINES}

$(EXECDIR)/stack_wgtcalc:stack_wgtcalc.o 
	$(FC) $(FFLAGS) -o $(EXECDIR)/stack_wgtcalc stack_wgtcalc.o

#	F90SAC requires special options to compile ...
f90sac_distrib.o: f90sac_distrib.F90
	$(FC) $(FFLAGS) $(F90FLAGS) $(F90SAC_FLAGS) -c $*.F90 

#
#     GMT PLOTTING SCRIPTS + OTHER SHELL SCRIPTS
#
$(EXECDIR)/sheba_plot_stackerr.gmt:sheba_plot_stackerr.gmt
	chmod +x sheba_plot_stackerr.gmt; cat sheba_plot_stackerr.gmt | sed 's/X1X/$(GMT_PREFIX)/g' > $(EXECDIR)/sheba_plot_stackerr.gmt

$(EXECDIR)/sheba_plot_errclu.gmt:sheba_plot_errclu.gmt
	chmod +x sheba_plot_errclu.gmt; cat sheba_plot_errclu.gmt | sed 's/X1X/$(GMT_PREFIX)/g' > $(EXECDIR)/sheba_plot_errclu.gmt

$(EXECDIR)/sheba_combine_plots.csh:sheba_combine_plots.csh
	chmod +x sheba_combine_plots.csh; \cp sheba_combine_plots.csh $(EXECDIR)

$(EXECDIR)/cleansheba:cleansheba
	chmod +x cleansheba; \cp cleansheba $(EXECDIR)

#
#     SAC Macro
#
$(MACRODIR)/split_sheba:sheba
	cp -f sheba $(MACRODIR)/split_sheba; cp -f sheba $(MACRODIR)

$(MACRODIR)/sheba:sheba
	cp -f sheba $(MACRODIR)

distrib:
	rm -rf ../SHEBA_distrib
	mkdir ../SHEBA_distrib
	cp -f *.f90 ../SHEBA_distrib
	cp -f *.F90 ../SHEBA_distrib
	cp -f *.f ../SHEBA_distrib
	cp -f *.gmt ../SHEBA_distrib
	cp -f *.csh ../SHEBA_distrib
	cp -f sheba ../SHEBA_distrib
	cp -f Makefile ../SHEBA_distrib
	cp -f INSTALL ../SHEBA_distrib
	cp -f cleansheba ../SHEBA_distrib	

clean:
	rm -f *.o *.mod

# 
# 	Compile Instuctions
#
%.o: %.F90
	$(FC) $(FFLAGS) $(F90FLAGS) -c $*.F90 
%.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -c $*.f90 
%.o: %.f
	$(FC) $(FFLAGS) $(F77FLAGS) -c $*.f

#
#	Dependencies
#
$(SUBROUTINES)    :	$(MODULES)		
