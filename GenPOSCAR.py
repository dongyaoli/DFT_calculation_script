from __future__ import division
import sys
import os
import re
import numpy as np

def read_POSCAR(file_name):
    f1 = open(file_name, 'r')
    lines = f1.readlines()
    scale_factor = float(lines[1].strip())
    cell = np.zeros((3,3))
    
    cell[0,:] = lines[2].strip().split()
    cell[1,:] = lines[3].strip().split()
    cell[2,:] = lines[4].strip().split()
    f1.close()
    
    return scale_factor, cell

def Vsweep_POSCAR(enlarge, c_a_ratio, element='Graphite'):
    a_length = 2.4604
    c1 = [1, 0, 0]
    c2 = [-1/2, np.sqrt(3)/2, 0]
    c3 = [0, 0, c_a_ratio]
    a_change = a_length * enlarge
        
    fout = open("POSCAR","w")
    fout.write(element + " Hex\n")
    fout.write(str(a_change) + "\n")
    fout.write(" %22.16f  %22.16f  %22.16f\n"%(c1[0],c1[1],c1[2]))
    fout.write(" %22.16f  %22.16f  %22.16f\n"%(c2[0],c2[1],c2[2]))
    fout.write(" %22.16f  %22.16f  %22.16f\n"%(c3[0],c3[1],c3[2]))
    fout.write(str(4) + "\n")
    fout.write("Selective Dynamics \n")
    fout.write("Direct\n")
    fout.write(" 0.00000 0.00000 0.0000 F F F \n")
    fout.write(" 0.66667 0.33333 0.0000 T T T \n")
    fout.write(" 0.33333 0.66667 0.5000 T T T \n")
    fout.write(" 0.66667 0.33333 0.5000 T T T \n")
    fout.close


def Calcu_volume(file_name):
    scale_factor, cell = read_POSCAR(file_name)
    axis1 = cell[0,:]
    axis2 = cell[1,:]
    axis3 = cell[2,:]
    volume = np.abs(np.dot(axis1, np.cross(axis2, axis3))) * scale_factor**3
    return volume
    
if (len(sys.argv) == 3):
   enlarge = float(sys.argv[1])
   c_a_ratio = float(sys.argv[2])
   Vsweep_POSCAR(enlarge, c_a_ratio)

