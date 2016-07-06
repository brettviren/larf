#!/bin/bash

rm -f larf.db
log="speed.log"
date > $log 2>&1

do_mesh () {
    cap=$1; shift
    echo "mesh capacitor$cap" >> $log 2>&1
    (time larf mesh capacitor$cap) >> $log 2>&1
    echo "mesh capacitor$cap" >> $log 2>&1
}

do_boundary () {
    res=$1; shift
    cap=$1; shift
    echo "boundary capacitor$cap" >> $log
    (time larf boundary -m $res capdrift) >> $log 2>&1
    echo "boundary capacitor$cap" >> $log
}

do_raster () {
    res=$1; shift
    rez=$1; shift
    echo "raster capacitor$rez" >> $log
    (time larf raster -b $res capacitor$rez) >> $log 2>&1
    echo "raster capacitor$rez" >> $log
}

do_mesh 1
do_mesh 5
do_mesh 10

do_boundary 1 1
do_boundary 2 5
do_boundary 3 10

do_raster 4 10
do_raster 5 10
do_raster 6 10

do_raster 4 100
do_raster 5 100
do_raster 6 100


