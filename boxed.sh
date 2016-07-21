#!/bin/bash

# get rid of annoying error:
# Error retrieving accessibility bus address: org.freedesktop.DBus.Error.ServiceUnknown: The name org.a11y.Bus was not provided by any .service files
export NO_AT_BRIDGE=1

mydir="$(dirname $(readlink -f $BASH_SOURCE))"

thisname="boxed"

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
set -e

doifneeded box			mesh -m box
doifneeded invbox		mesh -m invbox
doifneeded boxed		mesh -m boxed
doifneeded invboxed		mesh -m invboxed

doifneeded drift-boundary-box		boundary -b drift -m 1
doifneeded drift-boundary-invbox	boundary -b drift -m 2
doifneeded drift-boundary-boxed		boundary -b drift -m 3
doifneeded drift-boundary-invboxed	boundary -b drift -m 4

doifneeded bare-wires			mesh -m wires
doifneeded bare-wires-boundary		boundary -b drift -m 9

exit 0


doifneeded drift-field			raster -r coarse -b $driftbounds_id
driftfield_id=$ident



doifneeded uweight-boundary		boundary -b uweighting -m $meshid
ubounds_id=$ident

doifneeded uweight-field		raster -r fine -b $ubounds_id
ufield_id=$ident



