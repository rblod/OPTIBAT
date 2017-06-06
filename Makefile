SHELL = /bin/bash

all: WKB_MODEL EnKF-MPI-Waves Prep_Routines Tools

WKB_MODEL: 
	gmake -C WKB_MODEL
	
EnKF-MPI-Waves:
	gmake -C EnKF-MPI-Waves

Prep_Routines: 
	gmake -C Prep_Routines

Tools:

.PHONY: WKB_MODEL EnKF-MPI-Waves Prep_Routines Tools
