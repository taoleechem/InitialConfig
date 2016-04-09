from molecule import Molecule
import subprocess
import re


def to_gaussian(xyz_file, gaussian_filename, method = 'b3lyp', basis='6-31g', **args):
    handle = open(xyz_file)
    data = handle.readlines()
    handle.close()
    num = int(data[0].strip())
    handle2 = open(gaussian_filename, 'w')
    command = '%nprocshared=7 \n%mem=8GB \nchk=system.chk \n#p '
    command += '%s %s ' % (method, basis)
    for keyword in args.values():
        command += '%s ' % keyword
    command += ' \n\nTitle Card Required \n\n0 1\n'
    handle2.write(command)
    info = data[2:]
    for line in info:
        handle2.write(line)
    handle2.write('\n\n\n\n\n')
    handle2.flush()
    handle2.close()


def get_b3lyp_energy(xyz_file):
    to_gaussian(xyz_file, 'DATA/system.gjf')
    subprocess.call(['./g09_perform.sh'], shell=True)
    subprocess.call(['wait'], shell=True)
    handle = open('DATA/temp.txt')
    data = handle.readline().strip()
    energy = float(data)
    handle.close()
    return energy

if __name__ == '__main__':
    data = 'RMSDP=3.93D-09 MaxDP=1.06D-07 DE= 1.36D-12 OVMax= 1.90D-07\n SCF Done:  E(RB3LYP) =  -671.494323486     A.U. after   16 cycles    \n    Convg  =    0.3935D-08             -V/T =  2.0057\nKE= 6.676754218499D+02 PE=-3.316904475541D+03 EE= 1.102226639932D+03'
    s = (re.findall('E\(RB3LYP\) =\s+([-0-9.]+)', data))
    print float(s[0])

