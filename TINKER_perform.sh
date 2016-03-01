#!/bin/bash
cd '/home/litao/molecule-recognition/DATA/'
if
 ./analyze TinkerMolecules.xyz E >& TinkerMolecules.txt
then 
grep -i 'Intermolecular Energy :  ' <TinkerMolecules.txt |awk '{print $4}' >temp.txt
fi

