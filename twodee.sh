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

