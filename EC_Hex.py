from __future__ import division
import sys
import numpy as np

def read_CONTCAR(file_name):
    f1 = open(file_name, 'r')
    lines = f1.readlines()
    file_length = len(lines)
    scale_factor = lines[1].strip()
    cell = np.zeros((3,3))
    atom_detail = np.array([])
    
    cell[0,:] = lines[2].strip().split()
    cell[1,:] = lines[3].strip().split()
    cell[2,:] = lines[4].strip().split()
    
    atom_detail = lines[5:file_length]
    f1.close()
    
    return scale_factor, cell, atom_detail


def gen_POSCAR(file_name, deform, element='MoS2'):
    '''
    
    strain: contain strain tesor
    '''
    scale_factor, cell, atom_detail = read_CONTCAR(file_name)    
    
    c1 = cell[0,:]
    c2 = cell[1,:]
    c3 = cell[2,:]
    
    c1_s = np.dot(deform, c1)
    c2_s = np.dot(deform, c2)
    c3_s = np.dot(deform, c3)
    
    fout = open("POSCAR","w")
    fout.write(element + " Hex\n")
    fout.write(str(scale_factor) + "\n")
    fout.write(" %22.16f  %22.16f  %22.16f\n"%(c1_s[0],c1_s[1],c1_s[2]))
    fout.write(" %22.16f  %22.16f  %22.16f\n"%(c2_s[0],c2_s[1],c2_s[2]))
    fout.write(" %22.16f  %22.16f  %22.16f\n"%(c3_s[0],c3_s[1],c3_s[2]))
    for line in atom_detail:
        fout.write(line)
    fout.close
    
def gen_deform(epsilon, mode='c11-c12'):
    '''
    There are five kinds of strain tensor for Hexagonal crystal:
    mode = 
    c11-c12
    c11+c12
    c33
    c13
    c44
    epsilon is the strain
    '''
    strain = np.zeros((3,3))
    if mode == 'c11-c12':
        strain[0,0] = epsilon
        strain[1,1] = -epsilon
    elif mode == 'c11+c12':
        strain[0,0] = epsilon
        strain[1,1] = epsilon        
    elif mode == 'c44':
        strain[0,2] = epsilon
        strain[2,0] = epsilon
    elif mode == 'c13':
        strain[0,0] = epsilon
        strain[1,1] = epsilon
        strain[2,2] = -(2 * epsilon + epsilon * epsilon)/(1 + epsilon)**2
    elif mode == 'c33':
        strain[2,2] = epsilon       
    else:
        print 'Wrong mode'
        return
    deform = np.identity(3) + strain
    return deform

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

def Calcu_volume(file_name):
    scale_factor, cell = read_POSCAR(file_name)
    axis1 = cell[0,:]
    axis2 = cell[1,:]
    axis3 = cell[2,:]
    volume = np.abs(np.dot(axis1, np.cross(axis2, axis3))) * scale_factor**3
    return volume

if (len(sys.argv) == 4):
    file_name = (sys.argv[1])
    strain = float(sys.argv[2])
    mode = sys.argv[3]
    deform = gen_deform(strain, mode=mode)
    gen_POSCAR(file_name, deform)