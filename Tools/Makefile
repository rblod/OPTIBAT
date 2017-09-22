SHELL = /bin/bash

include ../make.inc

PROGS = EnKF_assemble\
checkice\
consistency\
fixmicom\
anamorphosis_state\
build_dp\
pert_param\
fixenkf\
ensemble_stats

PROGS = EnKF_assemble\
checkice\
consistency\
fixmicom\
anamorphosis_state\
build_dp\
pert_param\
fixenkf\
ensemble_stats

SFMAKEDEPEND := ../Config/sfmakedepend
DIR_LIB=lib
DIR_FILES=src
DIR_PPSRC=ppsrc
DIR_DOC=doc
DIR_OBJS=obj

FILENAMES =$(notdir $(basename $(wildcard src/*.F90))) $(notdir $(basename $(wildcard src/*.F))) $(notdir $(basename $(wildcard src/*.c)))
FILES=src/*
OBJS=$(addsuffix .o,$(addprefix $(DIR_OBJS)/,$(FILENAMES)))

all: depend prepprep $(DIR_LIB)/libtools.a $(PROGS)
		@echo
		@echo ===================================================
		@echo   Tools is ok
		@echo ===================================================
		@echo

prepprep:
		@echo
		@echo ===================================================
		@echo   Create working directories for Tools ...
		@echo ===================================================
		mkdir -p $(DIR_PPSRC)
		mkdir -p $(DIR_OBJS)
		mkdir -p $(DIR_LIB)
		@echo

$(DIR_LIB)/libtools.a: $(OBJS)
		@echo
		@echo ===================================================
		@echo   Compile Tools library...
		@echo ===================================================
		$(AR) -r $@ $(OBJS)
		ranlib $@
		@echo

EnKF_assemble: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_EnKF_assemble.f90
		$(FC) $(FFLAGS) -o EnKF_assemble  $(DIR_PPSRC)/p_EnKF_assemble.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)


consistency: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_consistency.f90
		$(FC) $(FFLAGS) -o consistency  $(DIR_PPSRC)/p_consistency.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

checkice: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_check_ice.f90
		$(FC) $(FFLAGS) -o checkice  $(DIR_PPSRC)/p_check_ice.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

fixmicom: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_fixmicom.f90
		$(FC) $(FFLAGS) -o fixmicom  $(DIR_PPSRC)/p_fixmicom.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

anamorphosis_state: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_anamorphosis_state.f90
		$(FC) $(FFLAGS) -o anamorphosis_state  $(DIR_PPSRC)/p_anamorphosis_state.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

build_dp: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_build_dp.f90
		$(FC) $(FFLAGS) -o build_dp  $(DIR_PPSRC)/p_build_dp.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

pert_param: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_pert_param.f90
		$(FC) $(FFLAGS) -o pert_param  $(DIR_PPSRC)/p_pert_param.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

fixenkf: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_fixenkf.f90
		$(FC) $(FFLAGS) -o fixenkf  $(DIR_PPSRC)/p_fixenkf.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

conservation: $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_conservation.f90
		$(FC) $(FFLAGS) -o conservation  $(DIR_PPSRC)/p_conservation.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

ensemble_stats:  $(DIR_LIB)/libtools.a $(DIR_PPSRC)/p_ensemble_stats.f90
		$(FC) $(FFLAGS) -o ensemble_stats  $(DIR_PPSRC)/p_ensemble_stats.f90 $(DIR_LIB)/libtools.a -I../EnKF-MPI-Waves/obj -L../EnKF-MPI-Waves/lib -lEnKF $(LIBS)

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.F90 Makefile MODEL.CPP
		$(CPP) $(INCS) -include MODEL.CPP $(CPPFLAGS) $(DIR_FILES)/$(*F).F90 > $(DIR_PPSRC)/$(*F).f90
		$(FC) $(FFLAGS) -c $(DIR_PPSRC)/$(*F).f90 -o $(DIR_OBJS)/$(*F).o -I../EnKF-MPI-Waves/obj $(INCS) 

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.F Makefile MODEL.CPP
		$(CPP) $(INCS) -include MODEL.CPP $(CPPFLAGS) $(DIR_FILES)/$(*F).F > $(DIR_PPSRC)/$(*F).f
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
		@echo   Create dependancies for Tools...
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
		@echo   Clean Tools ...
		@echo ===================================================
		$(RM) -r $(DIR_PPSRC)/* $(DIR_LIB)/libtools.a $(DIR_OBJS)/*.o $(DIR_OBJS)/*.mod *.mod $(PROGS) $(DIR_OBJS) $(DIR_LIB) $(DIR_PPSRC) $(DIR_DOC)
	
# DO NOT DELETE THIS LINE - used by make depend	
