SHELL = /bin/bash

all: WKB_MODEL EnKF-MPI-Waves Prep_Routines Tools

WKB_MODEL: 
	gmake -C WKB_MODEL
	
EnKF-MPI-Waves:
	gmake -C EnKF-MPI-Waves

Prep_Routines: 
	gmake -C Prep_Routines

Tools:
	gmake -C Tools

.PHONY: WKB_MODEL EnKF-MPI-Waves Prep_Routines Tools clean

clean:
		gmake -C WKB_MODEL clean
		gmake -C EnKF-MPI-Waves clean
		gmake -C Prep_Routines clean
		gmake -C Tools clean
 
