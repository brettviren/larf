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
mylarf velocity capvelo
mylarf step capacitor
