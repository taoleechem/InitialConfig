#!/bin/bash

cd /home/learner/InitialConfig/DATA/
nwchem NW.nw >& NW.out
wait
grep -i 'Total DFT energy' <NW.out |awk '{print $5}' >temp.txt
