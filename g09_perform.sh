#!/bin/bash

cd '/home/learner/InitialConfig/DATA/'
g09 system.gjf
wait
grep -i 'SCF Done:  E(RB3LYP) = ' <system.log |awk '{print $5}' >temp.txt
