#!/usr/bin/python
import subprocess
import time
from xyz import get_b3lyp_energy
import math

#  This script could perform multi-similar task automatically

#  1-Set variables and paths

taskFile = 'task_mymodel.txt'  # the template task file
SaveDir = 'SaveConfigs/mymodel_nc2_cc2h6cho/'
# 1000K -- 2.30 Angstrom, 500K -- 2.80 Angstrom
Temperature = 1000
ClosestDistance = 2.30
TotalOutputNum = 300
a_num = 19
b_num = 1
EachOutPutNum = []
a_fragments = 5
b_fragments = 2
initialConfigPath = 'InitiConfig/'
a_xyzfilename = 'nc2'
b_xyzfilename = 'cc2h6cho'
functional = "b3lyp"
basis = "6-31g"
othercommand = ''

HARTREE = 4.359744e-18
K_B_BOLTZMAN = 1.38037e-23


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
    fhandle.write(str(ClosestDistance)+'\n')    # B1_default_value
    fhandle.write('0 0'+"\n")
    fhandle.write('%s %s %s' %(functional, basis, othercommand))
    fhandle.flush()
    fhandle.close()


def initial_partition():
    energy_a = []
    for i in range(a_num):
        filename_a = initialConfigPath + a_xyzfilename + "_" + str(i + 1) + '.xyz'
        temp = get_b3lyp_energy(filename_a)
	print 'No.%s molecule A enenrgy: %s' %(i+1, temp)
        energy_a.append(temp)
    energy_b = []
    for j in range(b_num):
        filename_b = initialConfigPath + b_xyzfilename + "_" + str(j + 1) + '.xyz'
        temp = get_b3lyp_energy(filename_b)
        print 'No.%s molecule B energy: %s' %(j+1, temp)
        energy_b.append(temp)
    e1_min = min(energy_a)
    e2_min = min(energy_b)
    total_partition = 0
    Energy = []
    print 'Energy(unitless) is:'
    counts = 0
    for e1 in energy_a:
        counts += 1
        temp = []
        for e2 in energy_b:
            e = (e1 + e2 - e1_min - e2_min)*HARTREE / K_B_BOLTZMAN/Temperature
            temp.append(e)
            total_partition += math.exp(-1*e)
        print counts, temp
        Energy.append(temp)

    print 'partition function is:'
    Partition = []
    for line in Energy:
        temp = []
        for single in line:
            p = math.exp(-1*single)/total_partition
            temp.append(p)
            print p,
        Partition.append(temp)
        print ' '

    print 'Each pair output number is:'
    for line in Partition:
        temp = []
        for single in line:
            if single > 0.5:
                p = int(0.5*TotalOutputNum)
            elif single*TotalOutputNum <2:
                p = 2
            else:
                p = int(single*TotalOutputNum)
            print p,
            temp.append(p)
        print ' '
        EachOutPutNum.append(temp)






subprocess.call(["mkdir "+SaveDir], shell=True)


initial_partition()

count = 0
for i in range(a_num):
    for j in range(b_num):
        count += 1
        print '\n\n\n#!Perform No.',count, "Initial Configs Calculation!\n"
        filename_a = initialConfigPath + a_xyzfilename+"_" + str(i+1)
        filename_b = initialConfigPath + b_xyzfilename + "_" + str(j+1)
        GenerateTaskFile(taskFile, filename_a, filename_b, a_fragments, b_fragments, EachOutPutNum[i][j], Temperature)
        start = time.time()
        s1 = "./Gaussian_Config.exe "+taskFile+"\n wait"
        subprocess.call([s1], shell=True)
        print "\n\n This run time is:", time.time()-start,"\n\n"
        s2 = "rm SaveConfigs/final.xyz\n mv SaveConfigs/final.mol2 "+SaveDir+str(count)+".mol2"
        subprocess.call([s2],shell=True)
        s3 = "i=1\n for x in SaveConfigs/*.xyz\n do\n"+"	mv $x "+SaveDir+str(count)+"_$i.xyz\n i=$(($i+1))\n  done"
        subprocess.call([s3],shell=True)


