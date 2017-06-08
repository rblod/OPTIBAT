Répertoires:
- TIDE : élevation (deltaH) en Netcdf et en txt
- BATHY : les bathymétries 2D et 1D:
		- originale
		- en changeant la pente (up/down) 
		- en rajoutant un banc de sable (bank) 
- OBS : les obs parfaites sur un cycle de marée (+/- 1m) avec 40 échéances
- OBS_1D : les mêmes en 1D		

Les paramètres utilisés pour produire les obs:
BREAK_TG86
wkb_wwave:  amp [m], ang [deg], prd [s] , B_tg, gamma_tg
            0.7      -10.0      8.3       0.7    0.3

wkb_roller:  roller_sinb  roller_fraction
                  0.1          1.0
