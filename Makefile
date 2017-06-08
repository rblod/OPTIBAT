SHELL = /bin/bash

all: WKB_MODEL EnKF-MPI-Waves Prep_Routines Tools

WKB_MODEL: 
	make -C WKB_MODEL
	
EnKF-MPI-Waves:
	make -C EnKF-MPI-Waves

Prep_Routines: 
	make -C Prep_Routines

Tools:
	make -C Tools

.PHONY: WKB_MODEL EnKF-MPI-Waves Prep_Routines Tools clean

clean:
		make -C WKB_MODEL clean
		make -C EnKF-MPI-Waves clean
		make -C Prep_Routines clean
		make -C Tools clean
 
