#!/bin/bash

mydir="$(dirname $(readlink -f $BASH_SOURCE))"

thisname="uboone"

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

doifneeded mesh		mesh -m uboone
meshid=$ident

doifneeded boundary-drift		boundary -b drift -m $meshid
driftbounds_id=$ident
doifneeded raster-drift			raster -r fine -b $ident
driftid=$ident

doifneeded boundary-u7	-P domain=7	boundary -b weighting -m $meshid
u7bounds_id=$ident
doifneeded raster-u7			raster -r fine -b $ident
u7id=$ident

doifneeded boundary-v7	-P domain=20	boundary -b weighting -m $meshid
v7bounds_id=$ident
doifneeded raster-v7			raster -r fine -b $ident
v7id=$ident

doifneeded boundary-w7	-P domain=33	boundary -b weighting -m $meshid
w7bounds_id=$ident
doifneeded raster-w7			raster -r fine -b $ident
w7id=$ident

# zoom
doifneeded raster-drift			raster -r zoom -b $driftbounds_id
doifneeded raster-u7			raster -r zoom -b $u7bounds_id
doifneeded raster-v7			raster -r zoom -b $v7bounds_id
doifneeded raster-w7			raster -r zoom -b $w7bounds_id


