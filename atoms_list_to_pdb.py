
# coding: utf-8


import MDAnalysis as mda
import sys

list_name  = sys.argv[1]
pdb_name = sys.argv[2]

atoms_list = open(list_name, 'r').readlines()

atoms_num = []

for i in range(len(atoms_list)):
    atoms_ar = atoms_list[i].split()
    
    for j in range(len(atoms_ar)):
        try :
            atoms_num.append(int(atoms_ar[j]))
        except ValueError:
            pass

u = mda.Universe(pdb_name)

atoms_sel_list = ''
print('The QM zone has %i atoms (link atoms are not counted)' % len(atoms_num))

for i in range(len(atoms_num)):
    if i == len(atoms_num) - 1:
        atoms_sel_list = atoms_sel_list + 'bynum ' + str(atoms_num[i])
    else :
        atoms_sel_list = atoms_sel_list + 'bynum ' + str(atoms_num[i]) + ' or '

atoms_sel = u.select_atoms(atoms_sel_list)

atoms_sel.write(list_name + '.pdb')
