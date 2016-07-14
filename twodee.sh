#!/bin/bash

mydir="$(dirname $(readlink -f $BASH_SOURCE))"

twodee () {
    larf -c "$mydir/twodee.cfg" -s "$mydir/twodee.db" $@
}


doifneeded () {
    ident=$1 ; shift
    if [ -n "$(twodee list | grep ^${ident}'\b')" ] ; then
	return
    fi

    twodee $@
}

doifneeded 1 mesh twodee

set -x

# boundaries
time doifneeded 2              boundary -b drift drift
time doifneeded 3 -P domain=18 boundary -b weighting v6
time doifneeded 4              raster   -r coarse -b 3 v6
time doifneeded 5 -P domain=6  boundary -b weighting u6
time doifneeded 6              raster   -r coarse -b 5 u6
time doifneeded 7 -P domain=30 boundary -b weighting w6
time doifneeded 8              raster   -r coarse -b 7 w6
time doifneeded 9              raster   -r coarse -b 2 drift

time doifneeded 10             raster   -r fine -b 3 v6
time doifneeded 11             raster   -r fine -b 5 u6
time doifneeded 12             raster   -r fine -b 7 w6
time doifneeded 13             raster   -r fine -b 2 drift

