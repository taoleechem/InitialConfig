#!/bin/bash

count=1
for x in *.gjf
do
cp $x $count.xyz
count=$(($count +1))
done

