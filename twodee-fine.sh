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
doifneeded boundary-u7	-P domain=7	boundary -b weighting -m $ident
doifneeded raster-u7			raster -r fine -b $ident

# # boundaries
# oifneeded 2              boundary -b drift drift
# oifneeded 3 -P domain=18 boundary -b weighting v6
# oifneeded 4              raster   -r coarse -b 3 v6
# oifneeded 5 -P domain=6  boundary -b weighting u6
# oifneeded 6              raster   -r coarse -b 5 u6
# oifneeded 7 -P domain=30 boundary -b weighting w6
# oifneeded 8              raster   -r coarse -b 7 w6
# oifneeded 9              raster   -r coarse -b 2 drift

# oifneeded 10             raster   -r fine -b 3 v6
# oifneeded 11             raster   -r fine -b 5 u6
# oifneeded 12             raster   -r fine -b 7 w6
# oifneeded 13             raster   -r fine -b 2 drift

