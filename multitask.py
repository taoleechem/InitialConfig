#!/usr/bin/python
import subprocess

#This script could perform multi-similar task automatically

#1-Set variables and pathes

taskFile = 'task_mymodel.txt' ##the template task file
SaveDir = 'SaveConfigs/mymodel/'
EachOutPutNum = 50
Temperature = 500
ClosestDistance = 2.30
a_num = 19
b_num = 1
a_fragments = 5
b_fragments = 2
initialConfigPath = 'InitiConfig/'
a_xyzfilename = 'nc2'
b_xyzfilename = 'choch3'

## Name rule a_xyzfilename+'0'+'.xyz'

def GenerateTaskFile(tsk, a, b, na, nb, num_output=50, T=500):
    fhandle = open(tsk, 'w')
    fhandle.write('//explanation\n')
    fhandle.write(a+'.xyz\n')
    fhandle.write(b+'.xyz\n')
    fhandle.write(str(na)+"\n")
    fhandle.write(str(nb)+"\n")
    fhandle.write(initialConfigPath+ a_xyzfilename+"\n")
    fhandle.write(initialConfigPath+ b_xyzfilename+"\n")
    fhandle.write(str(num_output)+"\n")
    fhandle.write(str(T)+"\n")
    fhandle.write(str(2.30)+'\n')   #B1_default_value
    fhandle.write('0 0'+"\n")
    fhandle.flush()
    fhandle.close()

subprocess.call(["mkdir "+SaveDir], shell=True)

count = 0
for i in range(a_num):
    for j in range(b_num):
        count += 1
        print '#!Perform No.',count, "Initial Configs Calculation!"
        filename_a = initialConfigPath + a_xyzfilename+"_" + str(i+1)
        filename_b = initialConfigPath + b_xyzfilename +"_"+ str(j+1)
        GenerateTaskFile(taskFile, filename_a, filename_b,a_fragments, b_fragments, EachOutPutNum, Temperature)
	s1 = "./NWCHEM_Configs.exe "+taskFile+"\n wait"
        subprocess.call([s1], shell=True)
	s2 = "mv SaveConfigs/final.mol2 "+SaveDir+str(count)+".mol2"
	subprocess.call([s2],shell=True)
	s3 = "i=1\n for x in SaveConfigs/*.xyz\n do\n"+"	mv $x "+SaveDir+str(count)+"_$i.xyz\n i=$(($i+1))\n  done"
	subprocess.call([s3],shell=True)


