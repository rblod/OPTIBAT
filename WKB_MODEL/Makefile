SHELL = /bin/bash

include ../make.inc

SFMAKEDEPEND := ../Config/sfmakedepend

DIR_LIB=lib
DIR_FILES=src
DIR_PPSRC=ppsrc
DIR_DOC=doc
DIR_OBJS=obj

FILENAMES =$(notdir $(basename $(wildcard src/*.F90)))
FILES=$(addsuffix .F90,$(addprefix $(DIR_FILES)/,$(FILENAMES)))
OBJS=$(addsuffix .o,$(addprefix $(DIR_OBJS)/,$(FILENAMES)))

all: depend  prepwkb $(DIR_LIB)/libwkb.a  wkb.exe
	@echo
	@echo ===================================================
	@echo   WKB is ok
	@echo ===================================================
	@echo

prepwkb:
	@echo
	@echo ===================================================
	@echo   Create working directories for WKB ...
	@echo ===================================================
	mkdir -p $(DIR_PPSRC)
	mkdir -p $(DIR_OBJS)
	mkdir -p $(DIR_LIB)
	@echo

$(DIR_LIB)/libwkb.a:  $(OBJS)
	@echo
	@echo ===================================================
	@echo   Compile WKB library...
	@echo ===================================================
	$(AR) -r $@ $(OBJS)
	ranlib $@
	@echo

wkb.exe: $(DIR_LIB)/libwkb.a $(DIR_FILES)/wkb.F90 
	$(FC) -o wkb.exe $(DIR_PPSRC)/wkb.f90 -I$(DIR_OBJS) -I. $(DIR_LIB)/libwkb.a -I$(NCDF_ROOT)/include -L$(NCDF_ROOT)/lib -lnetcdf -lnetcdff


$(DIR_OBJS)/%.o : $(DIR_FILES)/%.F90
	$(CPP) $(CPPFLAGS) $(DIR_FILES)/$(*F).F90 > $(DIR_PPSRC)/$(*F).f90
	$(FC) $(FFLAGS) -c $(DIR_PPSRC)/$(*F).f90 -o $(DIR_OBJS)/$(*F).o        			

#--------------------------------------------------------------------------
#  Target to create dependeanies.
#--------------------------------------------------------------------------

MakeDepend: $(DIR_FILES)
	$(SFMAKEDEPEND) --silent --fext=F90 --file=- --moddir=$(DIR_OBJS) --objdir=$(DIR_OBJS)  $(FILES) > MakeDepend

depend: $(DIR_FILES) 
	@echo
	@echo ===================================================
	@echo   Create dependancies ...
	@echo ===================================================
	$(SFMAKEDEPEND) --silent --fext=F90 --file=- --moddir=$(DIR_OBJS) --objdir=$(DIR_OBJS)  $(FILES) > MakeDepend

#ifneq "$(MAKECMDGOALS)" "clean"
-include ./MakeDepend
#endif

#--------------------------------------------------------------------------
#  Target for doxygen
#--------------------------------------------------------------------------

doc: $(FILES) 
	@echo
	@echo ===================================================
	@echo   Make documentation ...
	@echo ===================================================
	doxygen Doxyfile
	
#--------------------------------------------------------------------------
#  Target to clean
#--------------------------------------------------------------------------

clean:
	@echo
	@echo ===================================================
	@echo   Clean ...
	@echo ===================================================
	$(RM) -r $(DIR_PPSRC)/* $(DIR_LIB)/libwkb.a $(DIR_OBJS)/*.o $(DIR_OBJS)/*.mod *.mod wkb.exe $(DIR_OBJS) $(DIR_LIB) $(DIR_PPSRC) $(DIR_DOC)
# DO NOT DELETE THIS LINE - used by make depend
