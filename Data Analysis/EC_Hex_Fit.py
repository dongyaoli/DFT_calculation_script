from __future__ import division
import numpy as np
import scipy.optimize as spopti
import matplotlib.pyplot as plt
import numpy.linalg as npla 

def read_CONTCAR(file_name):
    f1 = open(file_name, 'r')
    lines = f1.readlines()
    file_length = len(lines)
    scale_factor = float(lines[1].strip())
    cell = np.zeros((3,3))
    atom_detail = np.array([])
    atom_number = np.sum(float(x) for x in lines[6].split())
    
    cell[0,:] = lines[2].strip().split()
    cell[1,:] = lines[3].strip().split()
    cell[2,:] = lines[4].strip().split()
    
    atom_detail = lines[5:file_length]
    f1.close()
    
    return scale_factor, cell, atom_number, atom_detail


def gen_POSCAR(file_name, deform, element='Graphite'):
    '''
    
    strain: contain strain tesor
    '''
    scale_factor, cell, atome_number, atom_detail = read_CONTCAR(file_name)    
    
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


def read_strain_summary(summary_name):
    strain = np.array([])
    E = np.array([])
    Volume = np.array([])
    ionicStep = np.array([])
    f = open(summary_name, 'r')
    for line in f:
        line = line.strip()
        column = line.split()
        strain = np.append(strain, float(column[0]))
        E = np.append(E, float(column[7]))
        Volume = np.append(Volume, float(column[1]))
        ionicStep = np.append(ionicStep, float(column[3]))
    f.close()
    Volume = Volume * 10**-30  # Using m as unit for length
    E_final = np.array([])
    V_final = np.array([])
    strain_final = np.array([])
    index = 0
    for step in ionicStep:
        if step == 1:
            strain_final = np.append(strain_final, strain[index])
            E_final = np.append(E_final, E[index])
            V_final = np.append(V_final, Volume[index])
        index += 1
        
    return strain, E, Volume, strain_final, E_final, V_final


def Calcu_volume(file_name):
    scale_factor, cell, atome_number, atom_detail = read_CONTCAR(file_name)   
    axis1 = cell[0,:]
    axis2 = cell[1,:]
    axis3 = cell[2,:]
    volume = np.abs(np.dot(axis1, np.cross(axis2, axis3))) * scale_factor**3
    return volume


def parabola(x, a, b, c):
    return a * x**2 + b * x + c

## ********* Fitting *************
summary_name = 'SUMMARY_MoS2-1500-8x8x8-EC-PBE-TSvdW-c11+c12'
mode = 'c11+c12'
#mode = 'c11-c12'
#mode = 'c33'
#mode = 'c44'
#mode = 'c13'
if mode == 'c13':
    mode = 'c11+c12-4*c13+2*c33'

CONTCAR = 'CONTCAR_MoS2_TSvdW_final'
scale_factor, cell, atom_number, atom_detail = read_CONTCAR(CONTCAR)

V0 = Calcu_volume(CONTCAR)

lattice_a = npla.norm(cell[0,:] * float(scale_factor)) 
lattice_b = npla.norm(cell[1,:] * float(scale_factor)) 
lattice_c = npla.norm(cell[2,:] * float(scale_factor)) 
lattice_constant = "Lattice constants a=%.3f A and c=%.3f A \n" %(lattice_a, lattice_c)

# 1 eV/A^3 = 1.60217646E-19 / 1E-30 Pa = 160.217646 GPa
unit_conversion = 160.217646 



#strain = np.array([])
#E = np.array([])
#ionicStep = np.array([])
#E_final = np.array([])
#strain_final = np.array([])
#f = open(summary_name, 'r')
#for line in f:
#    line = line.strip()
#    column = line.split()
#    strain = np.append(strain, float(column[0]))
#    E = np.append(E, float(column[6]))
#    ionicStep = np.append(ionicStep, float(column[2]))
#f.close()    
#
#index = 0
#for step in ionicStep:
#    if step == 1:
#       strain_final = np.append(strain_final, strain[index])
#       E_final = np.append(E_final, E[index])
#    index += 1

strain, E, Volume, strain_final, E_final, V_final = read_strain_summary(summary_name)

p0 = np.array([40, 0, min(E_final)]) # Initial Guess
popt, _ = spopti.curve_fit(parabola, strain_final, E_final, p0=p0)

if mode == 'c44' or mode == 'c33':
    EC_result = (2 * popt[0] * unit_conversion)/V0 
else:
    EC_result = (popt[0] * unit_conversion)/V0

epsilon = np.linspace(min(strain), max(strain), 50)
energy = parabola(epsilon, popt[0], popt[1], popt[2])

p1 = plt.figure(1)
ax = p1.add_subplot(111)
plt.plot(strain_final, E_final, '*', label='data')
plt.plot(epsilon, energy, label='fitting')
plt.legend(loc = 'upper center')
plt.xlabel('Strain')
plt.ylabel('Energy(eV)')
plt.xlim([-0.013, 0.013])
ax.ticklabel_format(style = 'sci', useOffset=False)

notation = 'Volume: %.3f A^3\n' % (V0) \
        + 'Equlibruim Energy/atom: %.4f eV\n' % (popt[2]/float(atom_number))\
        + lattice_constant\
        + 'Plot to fit ' + mode + '\n' + mode + '=%.2f GPa \n'%(EC_result)

plt.text(1.01, 0.5, notation, ha='left', va='center', transform=ax.transAxes)       
plt.savefig('fitting-' + mode + '.pdf', bbox_inches='tight')