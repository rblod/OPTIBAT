#!/bin/bash


\cp ../WKB_MODEL/wkb.exe .
for i in 0{0..9} {10..30}
do
ln -sf bryfile_${i}.nc shoreface_bry.nc
ln -sf depth_file_${i}.nc shoreface_depth.nc
./wkb.exe
\rm shoreface_bry.nc shoreface_depth.nc
\mv shoreface_out.nc shoreface_out_${i}.nc

done