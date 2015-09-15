from __future__ import division
import numpy as np
import sys
import os
import glob
import matplotlib.pyplot as plt

current_path = os.path
sys.path.append(current_path)
all_file = [open(f) for f in glob.glob("SUMMARY_MoS2*")]

p1 = plt.figure(1)
ax = p1.add_subplot(111)
for f in all_file:
    ENCUT = np.array([])
    E = np.array([])
    name = f.name
    for line in f:
        line = line.strip()
        column = line.split()
        ENCUT = np.append(ENCUT, float(column[0]))
        E = np.append(E, float(column[5]))
    plt.plot(ENCUT, E, '*', alpha =1, label='k-point: ' + name[12:])
    f.close()

#plt.legend(loc = 'upper right')

plt.xlim([200, 1700])
plt.xlabel('Energy Cut(eV)')
plt.ylabel('Equilibrium Energy(eV)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
ax.ticklabel_format(style = 'sci', useOffset=False)
plt.savefig('Convergence.pdf', bbox_inches='tight')
