MPI = NO
include ../make.inc

VPATH = ../:./:TMP

PROGS = EnKF_assemble\
checkice\
consistency\
fixmicom\
anamorphosis_state\
build_dp\
pert_param\
fixenkf

all: $(PROGS)

ENKF_ASSEMBLE_OBJ = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_put_mod_fld.o\
p_EnKF_assemble.o

CONSISTENCY_OBJ = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_put_mod_fld.o\
m_get_mod_grid.o\
m_get_mod_fld.o\
nfw.o\
p_consistency.o

CHECK_ICE_OBJ = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_get_mod_grid.o\
p_check_ice.o


FIXMICOM_OBJ = qmpi.o\
nfw.o\
p_fixmicom.o



TESTX4_SRC_F90 = qmpi.F90\
m_parse_blkdat.F90\
m_get_mod_nrens.F90\
testX4.F90


TESTX4_OBJ = $(TESTX4_SRC_F90:.F90=.o)


ANAMORPHOSIS_STATE_SRC_F90=qmpi.F90\
mod_raw_io.F90\
m_parse_blkdat.F90\
m_get_micom_nrens.F90\
m_get_micom_grid.F90\
m_get_micom_dim.F90\
m_get_micom_fld.F90\
mod_measurement.F90\
m_sort2.F90\
m_random.F90\
m_ana_exp_enkf.F90\
nfw.F90\
p_anamorphosis_state.F90

ANAMORPHOSIS_STATE_OBJ=$(ANAMORPHOSIS_STATE_SRC_F90:.F90=.o)

BUILD_DP_OBJ = qmpi.o\
nfw.o\
mod_raw_io.o\
m_get_micom_fld.o\
m_sort2.o\
m_random.o\
m_ana_exp_enkf.o\
p_build_dp.o


MKENSEMBLE_F90= \
qmpi.F90\
nfw.F90\
func1D.F90\
func2D.F90\
m_zeroin.F90\
m_set_random_seed2.F90\
m_pseudo2D.F90\
m_fixsample.F90\
m_sample2D.F90\
m_param_ensemble.F90\
p_pert_param.F90

MKENSEMBLE_F77= eosdat.F \
mod_eosfun.F

MKENSEMBLE_OBJ = $(MKENSEMBLE_F77:.F=.o) $(MKENSEMBLE_F90:.F90=.o)

# MKENSEMBLE_OBJ = qmpi.o\
# nfw.o\
# eosdat.o\
# mod_eosfun.o\
# func1D.o\
# func2D.o\
# m_zeroin.o\
# m_pseudo2D.o\
# p_mkensemble.o


FIXENKF_OBJ = qmpi.o\
nfw.o\
eosdat.o\
mod_eosfun.o\
m_fixhycom_eco_metno.o\
p_fixenkf.o


CONSERV_OBJ = qmpi.o\
nfw.o\
eosdat.o\
mod_eosfun.o\
m_get_micom_dim.o\
m_get_micom_grid.o\
p_conservation.o

m_Generate_element_Si.o: m_parse_blkdat.o mod_measurement.o m_get_mod_fld.o m_insitu.o m_obs.o
m_insitu.o: nfw.o
m_get_mod_nrens.o: qmpi.o
m_obs.o : m_parameters.o m_uobs.o

$(CONSISTENCY_OBJ) $(FIXMICOM_OBJ) $(MKENSEMBLE_OBJ) : MODEL.CPP makefile

EnKF_assemble: $(ENKF_ASSEMBLE_OBJ)
	@echo "->EnKF_assemble"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../EnKF_assemble $(ENKF_ASSEMBLE_OBJ) 

consistency: $(CONSISTENCY_OBJ)
	@echo "->consistency"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../consistency $(CONSISTENCY_OBJ) $(LIBS)

checkice: $(CHECK_ICE_OBJ)
	   @echo "->checkice"
	   @cd ./TMP ; $(LD) $(LINKFLAGS) -o ../checkice $(CHECK_ICE_OBJ) $(LIBS)

fixmicom: $(FIXMICOM_OBJ)
	@echo "->fixmicom"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../fixmicom $(FIXMICOM_OBJ)

anamorphosis_state: $(ANAMORPHOSIS_STATE_OBJ)
	@echo "->anamorphosis_state"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../anamorphosis_state $(ANAMORPHOSIS_STATE_OBJ) $(LIBS)

build_dp: $(BUILD_DP_OBJ)
	@echo "->build_dp"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../build_dp $(BUILD_DP_OBJ)

pert_param: $(MKENSEMBLE_OBJ)
	@echo "->pert_param"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../pert_param $(MKENSEMBLE_OBJ) $(LIBS) 

fixenkf: $(FIXENKF_OBJ)
	@echo "->fixenkf"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../fixenkf $(FIXENKF_OBJ)

conservation: $(CONSERV_OBJ)
	@echo "-> conservation"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../conservation $(CONSERV_OBJ)

clean:
	@rm -f *.o *.mod TMP/*.* $(PROGS)

%.o: %.F90
	@echo "  $*".F90
	@rm -f ./TMP/$*.f90
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90
	@cd ./TMP; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90

%.o: %.F
	@echo "  $*".F
	@rm -f ./TMP/$*.f
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f
	@cd ./TMP; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f 2> /dev/null

%.o: %.c
	@echo "  $*".c
	@cd ./TMP ; $(CC) -c $(CFLAGS) -o $*.o -I.. -I../.. ../../$*.c
