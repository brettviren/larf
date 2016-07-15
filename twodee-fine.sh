#!/bin/bash

mydir="$(dirname $(readlink -f $BASH_SOURCE))"

thisname="twodee-fine"

thislarf () {
    larf -c "${mydir}/${thisname}.cfg" -s "${mydir}/${thisname}.db" $@
}

ident=0

doifneeded () {
    ident=$(( $ident + 1 ))
    name=$1 ; shift

    if [ -n "$(thislarf list | grep ^${ident}'\b')" ] ; then
	return 
    fi

    fullname="${thisname}-${name}"

    time thislarf $@ "$fullname"
    thislarf export -r $ident "${fullname}.vtk"
}

set -x

doifneeded mesh		mesh -m twodee
meshid=$ident

doifneeded boundary-u7	-P domain=7	boundary -b weighting -m $meshid
doifneeded raster-u7			raster -r fine -b $ident
u7id=$ident

doifneeded boundary-v7	-P domain=20	boundary -b weighting -m $meshid
doifneeded raster-v7			raster -r fine -b $ident
v7id=$ident

doifneeded boundary-w7	-P domain=33	boundary -b weighting -m $meshid
doifneeded raster-w7			raster -r fine -b $ident
w7id=$ident

doifneeded boundary-drift		boundary -b drift -m $meshid
doifneeded raster-drift			raster -r fine -b $ident
driftid=$ident


