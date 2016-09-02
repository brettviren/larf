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
outdir=/data1/users/bv/larf-run/${thisname}
mkdir -p $outdir
dbfile=$outdir/larf.db

thislarf () {
    larf -c "${mydir}/${thisname}.cfg" -s "${dbfile}" $@
}
doifneeded () {
    name=$1 ; shift
    if thislarf has -n $name ; then
	echo "Already have result $name"
	return
    fi
    t1=$(date +%s)
    echo $@ "$name"
    thislarf $@ "$name" >> $logfile 2>&1
    t2=$(date +%s)
    dt=$(( $t2 - $t1 ))
    echo "  elapsed: $dt seconds"
    echo "  elapsed: $dt seconds" >> $logfile
    vtkfile="${outdir}/${name}.vtk"
    echo "export to $vtkfile"
    thislarf export -n $name $vtkfile >> $logfile 2>&1
}
#set -x
set -e
######### end boilerplate


# wire geometry and surface
doifneeded wires			wires -w wireplanes
doifneeded surface			surface -s wirescreen	-w wires

exit 

# calculate boundary conditions
doifneeded drift-boundary		boundary -b drift	-s surface 
doifneeded uweight-boundary		boundary -b uweight	-s surface 
doifneeded vweight-boundary		boundary -b vweight	-s surface 
doifneeded wweight-boundary		boundary -b wweight	-s surface 


doifneeded drift-potential		raster -r coarse -b drift-boundary
doifneeded uweight-potential		raster -r coarse -b uweight-boundary
doifneeded vweight-potential		raster -r coarse -b vweight-boundary
doifneeded wweight-potential		raster -r coarse -b wweight-boundary


# make initial starting points
doifneeded vpoints			points -p uwires w sixty-wires
doifneeded upoints			points -p vwires w sixty-wires
doifneeded wpoints			points -p wwires w sixty-wires


# doifneeded vpaths			drift -d vline -b sixty-drift-boundary -w sixty-wires
# doifneeded wpaths			drift -d wline -b sixty-drift-boundary -w sixty-wires
# doifneeded upaths			drift -d uline -b sixty-drift-boundary -w sixty-wires

#doifneeded ucurrent			current -c dqdt -b sixty-uweight-boundary -d sixty-upaths
#doifneeded vcurrent			current -c dqdt -b sixty-vweight-boundary -d sixty-vpaths
#doifneeded wcurrent			current -c dqdt -b sixty-wweight-boundary -d sixty-wpaths



