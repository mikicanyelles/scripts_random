
# coding: utf-8


import MDAnalysis as mda
import sys

qm_name  = sys.argv[1]
pdb_name = sys.argv[2]

qm_list = open(qm_name, 'r').readlines()

qm_num = []

for i in range(len(qm_list)):
    qm_ar = qm_list[i].split()
    
    for j in range(len(qm_ar)):
        try :
            qm_num.append(int(qm_ar[j]))
        except ValueError:
            pass

u = mda.Universe(pdb_name)

qm_sel_list = ''
print('The QM zone has %i atoms (link atoms are not counted)' % len(qm_num))

for i in range(len(qm_num)):
    if i == len(qm_num) - 1:
        qm_sel_list = qm_sel_list + 'bynum ' + str(qm_num[i])
    else :
        qm_sel_list = qm_sel_list + 'bynum ' + str(qm_num[i]) + ' or '

qm_sel = u.select_atoms(qm_sel_list)

qm_sel.write(pdb_name[:-4] + '_QM.pdb')

