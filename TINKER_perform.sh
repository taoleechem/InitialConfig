#!/bin/bash
cd './DATA/'
if
 ./analyze TinkerMolecules.xyz E >& TinkerMolecules.txt
then 
grep -i 'Intermolecular Energy :  ' <TinkerMolecules.txt |awk '{print $4}' >temp.txt
fi

