SHELL = /bin/bash
MPI = YES

include ../make.inc

PROGS = EnKF
SFMAKEDEPEND := ../Config/sfmakedepend
DIR_LIB=lib
DIR_FILES=src
DIR_PPSRC=ppsrc
DIR_DOC=doc
DIR_OBJS=obj

FILENAMES =$(notdir $(basename $(wildcard src/*.F90))) $(notdir $(basename $(wildcard src/*.F))) $(notdir $(basename $(wildcard src/*.c)))
FILES=src/*
OBJS=$(addsuffix .o,$(addprefix $(DIR_OBJS)/,$(FILENAMES)))

all: depend prepEnKF $(DIR_LIB)/libEnKF.a $(PROGS)
		@echo
		@echo ===================================================
		@echo   EnKF is ok
		@echo ===================================================
		@echo

prepEnKF:
		@echo
		@echo ===================================================
		@echo   Create working directories for EnKF ...
		@echo ===================================================
		mkdir -p $(DIR_PPSRC)
		mkdir -p $(DIR_OBJS)
		mkdir -p $(DIR_LIB)
		@echo

$(DIR_LIB)/libEnKF.a: $(OBJS)
		@echo
		@echo ===================================================
		@echo   Compile EnKF library...
		@echo ===================================================
		$(AR) -r $@ $(OBJS)
		ranlib $@
		@echo

EnKF: $(DIR_LIB)/libEnKF.a $(DIR_PPSRC)/EnKF.f90
		$(FC) $(LINKFLAGS) $(FFLAGS) -o EnKF $(DIR_PPSRC)/EnKF.f90 $(DIR_LIB)/libEnKF.a $(LIBS) -I$(DIR_OBJS)

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.F90 Makefile MODEL.CPP
		$(CPP) -include MODEL.CPP $(CPPFLAGS) $(DIR_FILES)/$(*F).F90 > $(DIR_PPSRC)/$(*F).f90
		$(FC) $(FFLAGS) -c $(DIR_PPSRC)/$(*F).f90 -o $(DIR_OBJS)/$(*F).o

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.F Makefile MODEL.CPP
		$(CPP) -include MODEL.CPP $(CPPFLAGS) $(DIR_FILES)/$(*F).F > $(DIR_PPSRC)/$(*F).f
		$(FC) $(FFLAGS) $(F77FLG) -c $(DIR_PPSRC)/$(*F).f -o $(DIR_OBJS)/$(*F).o

$(DIR_OBJS)/%.o: $(DIR_FILES)/%.c Makefile MODEL.CPP
	$(CC) -c $(CFLAGS) -c $(DIR_FILES)/$(*F).c -o $(DIR_OBJS)/$(*F).o

#--------------------------------------------------------------------------
#  Target to create dependeanies.
#--------------------------------------------------------------------------

MakeDepend: $(DIR_FILES)
		$(SFMAKEDEPEND) --silent --fext=F90 --fext=F --fext=c --file=- --moddir=$(DIR_OBJS) --objdir=$(DIR_OBJS)  $(FILES) > MakeDepend

depend: $(DIR_FILES)
		@echo
		@echo ===================================================
		@echo   Create dependancies for EnKF ...
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
		@echo   Clean ...
		@echo ===================================================
		$(RM) -r $(DIR_PPSRC)/* $(DIR_LIB)/libprep_obs.a $(DIR_OBJS)/*.o $(DIR_OBJS)/*.mod *.mod $(PROGS) $(DIR_OBJS) $(DIR_LIB) $(DIR_PPSRC) $(DIR_DOC)
	
# DO NOT DELETE THIS LINE - used by make depend	
