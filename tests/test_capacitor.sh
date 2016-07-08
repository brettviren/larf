#!/bin/bash



dbfile=test_capacitor.db
rm -f $dbfile

mylarf () {
    larf -s $dbfile $@
}

set -e
set -x

mylarf mesh capacitor
mylarf boundary capdrift
mylarf raster -r capacitor capdrift
mylarf velocity capacitor
mylarf step capacitor

mylarf -P domain=1 boundary -b weighting capbound1
mylarf -P domain=2 boundary -b weighting capbound2

mylarf raster -b capbound1 -r capacitor capweight1
mylarf raster -b capbound2 -r capacitor capweight2

mylarf current -c capacitor -s capacitor -w capweight1 capcur1
mylarf current -c capacitor -s capacitor -w capweight2 capcur2
