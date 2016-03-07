#!/bin/bash

cd '/home/learner/InitialConfig/DATA/'
nwchem NW.nw >& NW.out
wait
grep -i 'scf energy' <NW.out |awk '{print $5}' >temp.txt
