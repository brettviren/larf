#!/bin/bash

# Quick and dirty idempotent workflow system.  As long as "doifneeded"
# commands are always *added* to the end, everything will remain
# idempotent.  The "ident" variable is incremented to hold the result
# ID made by each doifneeded call.  It can be saved to a variable for
# later reference.

# Jump to "end boilerplate" for servicible parts


######### begin boilerplate
# get rid of annoying error:
##
#   Error retrieving accessibility bus address:
#   org.freedesktop.DBus.Error.ServiceUnknown: The name org.a11y.Bus
#   was not provided by any .service files
##
export NO_AT_BRIDGE=1

mydir="$(dirname $(readlink -f $BASH_SOURCE))"
thisname="$(basename $BASH_SOURCE .sh)"
logfile="$mydir/${thisname}.log"
thislarf () {
    date
    t1=$(date +%s)
    larf -c "${mydir}/${thisname}.cfg" -s "${mydir}/${thisname}.db" $@
    date
    t2=$(date +%s)
    dt=$(( $t2 - $t1 ))
    echo "elapsed: $dt"
}
ident=0
doifneeded () {
    ident=$(( $ident + 1 ))
    name=$1 ; shift
    if [ -n "$(thislarf list | grep ^${ident}'\b')" ] ; then
	return 
    fi
    fullname="${thisname}-${name}"

    thislarf $@ "$fullname" >> $logfile 2>&1
    thislarf export -r $ident "${fullname}.vtk" >> $logfile 2>&1
}
set -x
set -e
######### end boilerplate


doifneeded mesh				mesh -m $thisname
mesh_nominal=$ident

doifneeded drift-boundary		boundary -b drift -m $mesh_nominal
bound_drift=$ident

doifneeded drift-field			raster -r coarse -b $bound_drift
field_drift=$ident

doifneeded mesh				mesh -m capped
mesh_wires=$ident

doifneeded drift-boundary-wires		boundary -b drift -m $mesh_wires
bound_drift_wires=$ident

doifneeded drift-field-wires		raster -r coarse -b $bound_drift_wires
field_drift_wires=$ident

# changed wire lcar to 0.5mm

doifneeded mesh				mesh -m capped
mesh_wires=$ident

doifneeded drift-boundary-wires		boundary -b drift -m $mesh_wires
bound_drift_wires=$ident

doifneeded drift-field-wires		raster -r coarse -b $bound_drift_wires
field_drift_wires=$ident

# remove super-anode plane

doifneeded mesh				mesh -m cappedone
mesh_wires=$ident

doifneeded drift-boundary-wires		boundary -b drift -m $mesh_wires
bound_drift_wires=$ident

doifneeded drift-field-wires		raster -r coarse -b $bound_drift_wires
field_drift_wires=$ident

doifneeded drift-field-wires		raster -r fine -b $bound_drift_wires
field_drift_wires_fine=$ident

# U

doifneeded uweight-boundary-wires	boundary -b uweight -m $mesh_wires
bound_uweight_wires=$ident

doifneeded uweight-field-wires		raster -r coarse -b $bound_uweight_wires
field_uweight_wires=$ident

doifneeded uweight-field-wires-fine	raster -r fine -b $bound_uweight_wires
field_uweight_wires_fine=$ident

# V

doifneeded vweight-boundary-wires	boundary -b vweight -m $mesh_wires
bound_vweight_wires=$ident

doifneeded vweight-field-wires		raster -r coarse -b $bound_vweight_wires
field_vweight_wires=$ident

doifneeded vweight-field-wires-fine	raster -r fine -b $bound_vweight_wires
field_vweight_wires_fine=$ident

# W

doifneeded wweight-boundary-wires	boundary -b wweight -m $mesh_wires
bound_wweight_wires=$ident

doifneeded wweight-field-wires		raster -r coarse -b $bound_wweight_wires
field_wweight_wires=$ident

doifneeded wweight-field-wires-fine	raster -r fine -b $bound_wweight_wires
field_wweight_wires_fine=$ident


# velocity

doifneeded velocity			velocity -r $field_drift_wires
drift_velo=$ident
doifneeded velocity-fine		velocity -r $field_drift_wires_fine
drift_velo_fine=$ident

# instantaneous current

doifneeded ucurrent			current -v $drift_velo -w $field_uweight_wires 
ucurrent=$ident
doifneeded vcurrent			current -v $drift_velo -w $field_vweight_wires 
vcurrent=$ident
doifneeded wcurrent			current -v $drift_velo -w $field_wweight_wires 
wcurrent=$ident

doifneeded ucurrent-fine		current -v $drift_velo_fine -w $field_uweight_wires_fine
doifneeded vcurrent-fine		current -v $drift_velo_fine -w $field_vweight_wires_fine 
doifneeded wcurrent-fine		current -v $drift_velo_fine -w $field_wweight_wires_fine

doifneeded uwaveform			waveforms -w patch -v $drift_velo -c $ucurrent
doifneeded vwaveform			waveforms -w patch -v $drift_velo -c $vcurrent
doifneeded wwaveform			waveforms -w patch -v $drift_velo -c $wcurrent

exit 0

