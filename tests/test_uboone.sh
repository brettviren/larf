#!/bin/bash

# Run larf on uboone wire patterns.

dbfile=test_uboone.db
rm -f $dbfile

mylarf () {
    larf -s $dbfile $@
}

set -e
set -x

mylarf mesh uboone
mylarf boundary capdrift
mylarf raster -r uboone -b drift drift
mylarf velocity -r drift drift
mylarf step -v drift uboone

mylarf -P domain=6 boundary -b weighting uboone-u6
mylarf -P domain=18 boundary -b weighting uboone-v6
mylarf -P domain=30 boundary -b weighting uboone-w6

mylarf raster -b uboone-u6 -r uboone uboone-u6
mylarf raster -b uboone-v6 -r uboone uboone-v6
mylarf raster -b uboone-w6 -r uboone uboone-w6

mylarf current -c uboone -s uboone -w uboone-u6 uboone-u6
mylarf current -c uboone -s uboone -w uboone-v6 uboone-v6
mylarf current -c uboone -s uboone -w uboone-w6 uboone-w6


