SHELL = /bin/bash

include ../make.inc

PROGS=prep_obs\
prep_obs_waves

SFMAKEDEPEND := ../Config/sfmakedepend
DIR_LIB=lib
DIR_FILES=src
DIR_PPSRC=ppsrc
DIR_DOC=doc
DIR_OBJS=obj

FILENAMES =$(notdir $(basename $(wildcard src/*.F90))) $(notdir $(basename $(wildcard src/*.F))) $(notdir $(basename $(wildcard src/*.c)))
FILES=src/*
OBJS=$(addsuffix .o,$(addprefix $(DIR_OBJS)/,$(FILENAMES)))

all: depend prepprep $(DIR_LIB)/libprep_obs.a $(PROGS)
		@echo
		@echo ===================================================
		@echo   Prep_Routines is ok
		@echo ===================================================
		@echo

prepprep:
		@echo
		@echo ===================================================
		@echo   Create working directories for Prep_Routines ...
		@echo ===================================================
		mkdir -p $(DIR_PPSRC)
		mkdir -p $(DIR_OBJS)
		mkdir -p $(DIR_LIB)
		@echo

$(DIR_LIB)/libprep_obs.a: $(OBJS)
		@echo
		@echo ===================================================
		@echo   Compile Prep_Routines library...
		@echo ===================================================
		$(AR) -r $@ $(OBJS)
		ranlib $@
		@echo

prep_obs: $(DIR_LIB)/libprep_obs.a $(DIR_PPSRC)/p_prep_obs.f90
		$(FC) $(FFLAGS) -o prep_obs  $(DIR_PPSRC)/p_prep_obs.f90 $(DIR_LIB)/libprep_obs.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

prep_obs_OOPS: $(DIR_LIB)/libprep_obs.a $(DIR_PPSRC)/p_prep_obs_OOPS.f90
		$(FC) $(FFLAGS) -o prep_obs_OOPS  $(DIR_PPSRC)/p_prep_obs_OOPS.f90 $(DIR_LIB)/libprep_obs.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

prep_obs_waves: $(DIR_LIB)/libprep_obs.a $(DIR_PPSRC)/p_prep_obs_waves.f90
		$(FC) $(FFLAGS) -o prep_obs_waves  $(DIR_PPSRC)/p_prep_obs_waves.f90 $(DIR_LIB)/libprep_obs.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.F90 Makefile MODEL.CPP
		$(CPP) -include MODEL.CPP $(CPPFLAGS) $(DIR_FILES)/$(*F).F90 > $(DIR_PPSRC)/$(*F).f90
		$(FC) $(FFLAGS) -c $(DIR_PPSRC)/$(*F).f90 -o $(DIR_OBJS)/$(*F).o -I../EnKF-MPI-Waves/obj

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.F Makefile MODEL.CPP
		$(CPP) -include MODEL.CPP $(CPPFLAGS) $(DIR_FILES)/$(*F).F > $(DIR_PPSRC)/$(*F).f
		$(FC) $(FFLAGS) $(F77FLG)-c $(DIR_PPSRC)/$(*F).f -o $(DIR_OBJS)/$(*F).o

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.c Makefile MODEL.CPP
	$(CC) -c $(CFLAGS) -c $(DIR_FILES)/$(*F).c -o $(DIR_OBJS)/$(*F).o -I../EnKF-MPI-Waves/src

#--------------------------------------------------------------------------
#  Target to create dependeanies.
#--------------------------------------------------------------------------

MakeDepend: $(DIR_FILES)
		$(SFMAKEDEPEND) --silent --fext=F90 --fext=F --fext=c --file=- --moddir=$(DIR_OBJS) --objdir=$(DIR_OBJS)  $(FILES) > MakeDepend

depend: $(DIR_FILES)
		@echo
		@echo ===================================================
		@echo   Create dependancies for Prep_Routines...
		@echo ===================================================
		$(SFMAKEDEPEND) --silent --fext=F90 --fext=F --fext=c --file=- --moddir=$(DIR_OBJS) --objdir=$(DIR_OBJS)  $(FILES) > MakeDepend

#ifneq "$(MAKECMDGOALS)" "clean"
-include ./MakeDepend
#endif

#--------------------------------------------------------------------------
#  Target to clean
#--------------------------------------------------------------------------

clean:
		@echo
		@echo ===================================================
		@echo   Clean Prep_Routines ...
		@echo ===================================================
		$(RM) -r $(DIR_PPSRC)/* $(DIR_LIB)/libprep_obs.a $(DIR_OBJS)/*.o $(DIR_OBJS)/*.mod *.mod $(PROGS) $(DIR_OBJS) $(DIR_LIB) $(DIR_PPSRC) $(DIR_DOC)
	
# DO NOT DELETE THIS LINE - used by make depend	
