from __future__ import division
import numpy as np
import sys
import os
import linecache as linecache
import scipy.optimize as spopti
import matplotlib.pyplot as plt

#path = "D:\Research\Projects\FP calculation\Elastic Constant\MoS2\Vsweep"
#os.chdir(path)

def Birch_Murnaghan(V, E0, V0, B0, dB0):
    '''
    E0: equilibrium energy
    V0: equilibrium volume
    B0: bulk modulus
    V:  deformed volume
    dB0: derivative of the bulk modulus with respect to pressure
    '''
    E = E0 + 9 * V0 * B0 *(((V0/V)**(2/3) - 1)**3 * dB0 + \
    ((V0/V)**(2/3) - 1)**2 * (6 - 4 * (V0/V)**(2/3)))/16
    return E
    
def Murnaghan(V, E0, V0, B0, dB0):
    '''
    E0: equilibrium energy
    V0: equilibrium volume
    B0: bulk modulus
    V:  deformed volume
    dB0: derivative of the bulk modulus with respect to pressure
    '''
    E = E0 + B0 * V0 * (((V0/V)**(dB0 - 1))/(dB0 - 1) + V/V0 - dB0/(dB0 - 1))/dB0
    return E

def read_summary_vsweep(summary):
    scaling_factor = np.array([])
    E = np.array([])
    Volume = np.array([])
    ionicStep = np.array([])
    f1 = open(summary, 'r')
    for line in f1:
        line = line.strip()
        column = line.split()
        scaling_factor = np.append(scaling_factor, float(column[0]))
        Volume = np.append(Volume, float(column[1]))
        ionicStep = np.append(ionicStep, float(column[2]))
        E = np.append(E, float(column[7]))
    f1.close()
#    E = E * eV
#    Volume = Volume * 10**-30  # Using m as unit for length
    E_final = np.array([])
    V_final = np.array([])
    index = 0
    
    # Store the energy with 1 ionic step
    for step in ionicStep:
        if step == 1:
            E_final = np.append(E_final, E[index])
            V_final = np.append(V_final, Volume[index])
        index += 1
    
    return E, Volume, E_final, V_final, scaling_factor
    
    
def fit_EOS(summary, POSCAR, ENCUT='ENCUT=1500\n', KPOINT='k-point 8x8x8\n'):
    eV = 1.60217 * 10**-19
    E, Volume, E_final, V_final, scaling_factor = read_summary_vsweep(summary)
    Volume = Volume * 10**-30  # Using m as unit for length
    V_final = V_final * 10**-30
    E = E * eV
    E_final = E_final * eV
           
    p0 = [min(E), Volume[0], 50 * 10**9, 8] # initial guess
#    popt, _ = spopti.curve_fit(Birch_Murnaghan, V_min, E_min, p0=p0)  
    popt, _ = spopti.curve_fit(Murnaghan, V_final, E_final, p0=p0)  
    v_plot = np.linspace(min(Volume), max(Volume), 50)
    energy = Birch_Murnaghan(v_plot, popt[0], popt[1], popt[2], popt[3])
    
    axis1 = [float(x) for x in linecache.getline(POSCAR, 3).split()]
    axis2 = [float(x) for x in linecache.getline(POSCAR, 4).split()]
    axis3 = [float(x) for x in linecache.getline(POSCAR, 5).split()] 
    atom_number = np.sum([float(x) for x in linecache.getline(POSCAR, 7).split()])
    V0 = np.dot(axis1, np.cross(axis2, axis3)) # Triple product of three axis to get the volume
    optimScale = (popt[1] / (V0 * 10**-30))**(1/3)   
    notation = 'Volume: %.3f A^3\n' % (popt[1]/10**-30) \
            + 'Atom in cell: %.3f\n' % (atom_number) \
            + 'Equlibruim Energy/atom: %.4f eV\n' % (popt[0]/eV/atom_number)\
            + 'Bulk modulus: %.2f GPa\n' % (popt[2]/10**9)\
            + 'Derivitive of modulus: %.3f \n' %(popt[3])\
            + 'Optimum Scale Factor: %.6f \n\n' %(optimScale)\
            + ENCUT + KPOINT + '\n'
    
    return V_final, E_final/eV/atom_number, Volume, E/eV/atom_number, v_plot, energy/eV/atom_number, notation


#Notes = 'PBE + TS vdW approximiation\n'

#file_name1 = 'SUMMARY_Graphite-16x16x8-1000-PBE-TSvdW-2.671ratio'
#POSCAR1 = 'POSCAR_Graphite-16x16x8-1000-PBE-TSvdW-2.671ratio'
#Volume1, E_per1, v_plot1, E_plot1, notation1 = fit_EOS(file_name1, POSCAR1)


#file_name2 = 'SUMMARY_Graphite-16x16x8-1000-PBE-TSvdW-2.726ratio'
#POSCAR2 = 'POSCAR_Graphite-16x16x8-1000-PBE-TSvdW-2.726ratio'
#Volume2, E_per2, v_plot2, E_plot2, notation2 = fit_EOS(file_name2, POSCAR2)

file_name3 = 'SUMMARY_MoS2-1500-8x8x8'
POSCAR3 = 'POSCAR_MoS2-1500-8x8x8-1-point-2'
V_min3, E_minper3, Volume3, E_per3, v_plot3, E_plot3, notation3 = fit_EOS(file_name3, POSCAR3)


p1 = plt.figure(1)
ax = p1.add_subplot(111)
#plt.plot(Volume1, E_per1, '+', label='Relaxiation Data')
#plt.plot(v_plot1, E_plot1, label='Fit the points with 1 ionic step')
#plt.text(1.01, 0.8, notation1, ha='left', va='center', transform=ax.transAxes)

#plt.plot(Volume2, E_per2, '+', label='Relaxiation Data')
#plt.plot(v_plot2, E_plot2, label='Fit the points with 1 ionic step')
#plt.text(1.01, 0.3, notation2, ha='left', va='center', transform=ax.transAxes)

plt.plot(Volume3, E_per3, '+', label='Different Steps')
plt.plot(V_min3, E_minper3, '*', label='Energy-1 ionic step')
plt.plot(v_plot3, E_plot3, label='Fitting')
plt.text(1.01, 0.3, notation3, ha='left', va='center', transform=ax.transAxes)

plt.legend(loc = 'upper center')
plt.xlabel('Volume (m^3)')
plt.ylabel('Energy/atom(eV)')
ax.ticklabel_format(style = 'sci', useOffset=False)



#cell_length = '%.3f, %.3f, %.3f' % (optimScale * np.linalg.norm(axis1), \
#    optimScale * np.linalg.norm(axis2), \
#    optimScale * np.linalg.norm(axis3))
#    
#notation = 'Length of cell vectors:\n' + cell_length + ' A\n'\
#        + 'Equlibruim Energy/atom: ' + '%.4f' % (popt[0]/eV/atom_number) + ' eV\n'\
#        + 'Bulk modulus: ' + '%.2f' % (popt[2]/10**9) + ' GPa\n'\
#        + 'Derivitive of modulus: ' + '%.3f' %(popt[3]) + '\n\n'\
#        + ENCUT + KPOINT
        

#plt.savefig('fitting.pdf', bbox_inches='tight')
