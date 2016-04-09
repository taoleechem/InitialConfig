# Created by Tao Li

import numpy as np
import re


class Molecule:
    number = 0
    name = []
    corr = []

    def __init__(self, file = '', xyz_format = True, gaussian_format = False):
        if file:
            if gaussian_format:
                print 'generate g09 format'
            elif xyz_format:
                self.from_xyz(file)
        else:
            self.number = 0
            self.name = []
            self.corr = []

    def __eq__(self, other):
        self.number = other.number
        self.name = other.name
        self.corr = other.corr

    def from_xyz(self, filename):
        handle = open(filename)
        n = int(handle.readline().strip().split()[0])
        handle.readline()
        count = 0
        while count < n:
            count += 1
            words = handle.readline().strip().split()
            self.number += 1
            self.name.append(words[0])
            self.corr.append([float(words[1]), float(words[2]), float(words[3])])
        handle.close()

    def to_xyz(self, filename):
        handle = open(filename, 'w')
        handle.write(str(self.number) + '\n\n')
        for atom, position in zip(self.name, self.corr):
            line = '%s %.6f %.6f %.6f\n' % (atom, position[0], position[1], position[2])
            handle.write(line)
        handle.flush()
        handle.close()

    def to_gaussian(self, filename, method = 'b3lyp', basis = '6-31g', **args):
        handle = open(filename, 'w')
        command = '%nprocshared=7 \n%mem=8GB \nchk=system.chk \n#p '
        command += '%s %s ' % (method, basis)
        for keyword in args.values():
            command += '%s ' % keyword
        command += ' \n\nTitle Card Required \n\n0 1\n'
        handle.write(command)
	print command,
        for atom, position in zip(self.name, self.corr):
            line = ' %s %.7f %.7f %.7f\n' %(atom, position[0], position[1], position[2])
            handle.write(line)
            print line,
        handle.write('\n\n\n\n\n')
        print '\n\n\n\n\n'
        handle.flush()
        handle.close()

    def __str__(self):
        x = '%d \n\n' % self.number
        for atom, position in zip(self.name, self.corr):
            line = '%s %.6f %.6f %.6f\n' % (atom, position[0], position[1], position[2])
            x += line
        return x


def test():
    x = Molecule('1.xyz')
    x.to_gaussian('1.gjf', method='b3lyp', basis='6-31g')


if __name__ == '__main__':
    test()
