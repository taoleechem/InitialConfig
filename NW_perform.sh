#!/bin/bash

cd './DATA/'
nwchem NW.nw >& NW.out
wait
grep -i 'DTF energy' <NW.out |awk '{print $5}' >temp.txt
