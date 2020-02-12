#!/bin/bash

name=$1
dx=$2
dt=$3
dtmax=$dt
dtmin=$4
nspool=$5
stack=$6

mkdir -p $name
(cd $name; python ${SCRIPTS_DIR}/create_grid.py $dx)
cp param.nml $name
sed -i "/dt =/c\  dt = ${dt}" $name/param.nml
sed -i "/dtb_max =/c\  dtb_max = ${dtmax}" $name/param.nml
sed -i "/dtb_min =/c\  dtb_min = ${dtmin}" $name/param.nml
sed -i "/ nspool =/c\  nspool = ${nspool}" $name/param.nml
sed -i "/ ihfskip =/c\  ihfskip = ${stack}" $name/param.nml
cp vgrid.in $name
mkdir -p $name/outputs

echo "  run model"
(cd $name; time ${SCHISM_BUILD_DIR}/bin/pschism_TVD-VL &> log.txt)
#; ${SCHISM_BUILD_DIR}/bin/combine_output11 1 1)
