# IDEAS

#Collect two arguments: topology and coordinates (either pdb or inpcrd or whatever)


# Import packages
from parmed import load_file
from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError
import sys


# File
if len(sys.argv) < 2:
    print('Input the topology/parameters and the coordinates as arguments. A name for the output file can be specified as the third argument (without an extension).')

elif len(sys.argv) == 2:
    name_pdb = sys.argv[1]
    only_list = True

elif len(sys.argv) >= 3:
    name_prmtop = sys.argv[1]
    #name_inpcrd = sys.argv[2]
    only_list = False

    if str(sys.argv[2]).find('.pdb') == (len(sys.argv[2]) - 4):
        name_pdb = sys.argv[2]
        name_inpcrd = None
    elif str(sys.argv[2]).find('.inpcrd') == (len(sys.argv[2]) - 7):
        name_inpcrd = sys.argv[2]
        name_pdb = None

if only_list == False:
    # Ask for Ligand and Radius
    ligand = input('Type the name or the three letters code (only if it is a non-standard residue) of the central ligand: ')
    while True:
        quest = input('Is %s correct ([y]/n)? ' % ligand)
        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS'):
            break
        elif quest in ('no', 'NO', 'No', 'nO'):
            ligand = input('Type the name or the three letters code (only if it is a non-standard residue) of the central ligand: ')
            continue
        else :
            print('Answer \'yes\' or \'no\'.')


    while True:
        try:
            radius = float(input('Which is the desired radius for the water drop (in Å)? '))
            break
        except ValueError:
            print('Type only a number.')
            continue
    while True:
        quest = input('Is %s correct ([y]/n)? ' % radius)
        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS'):
            break
        elif quest in ('no', 'NO', 'No', 'nO'):
            radius = input('Which is the desired radius for the water drop (in  Å)? ')
            continue
        else :
            print('Answer \'yes\' or \'no\'.')
    
    # Output files names

    if len(sys.argv) < 4:
        name_prmtop_out = name_prmtop[:-7] + '.cropped.prmtop'
        if name_inpcrd == None:
            name_inpcrd_out = name_prmtop[:-7] + '.cropped.inpcrd'
        else :
            name_inpcrd_out = name_inpcrd[:-7] + '.cropped.inpcrd'

    elif len(sys.argv) == 4:
        name_prmtop_out = sys.argv[3] + '.prmtop'
        name_inpcrd_out = sys.argv[3] + '.inpcrd'
        

    # Crop the system

    if name_pdb != None:
        topology = load_file(name_prmtop, xyz=name_pdb)
    elif name_inpcrd != None:
        topology = load_file(name_prmtop, xyz=name_inpcrd)
    topology.box = None
    topology.strip(':WAT&:%s<@%s' % (ligand, radius))
    topology.strip(':Na+&:%s<@%s' % (ligand, radius))
    topology.strip(':Cl-&:%s<@%s' % (ligand, radius))

    topology.write_parm(name_prmtop_out)
    topology.save(name_inpcrd_out)

    print('Cropped topology and coordinates have been saved as \'%s\' and \'%s\'' % (name_prmtop_out, name_inpcrd_out))
    ##################


# Generate active atoms list
if only_list == False:
    u_set_act = Universe(name_prmtop_out, name_inpcrd_out)
    name = name_prmtop_out[:-7]
elif only_list == True:
    u_set_act = Universe(name_pdb)
    name = name_pdb[:-4]
    radius = 100000000

while True:
    try :
        carbon = input("Type the number of the central atom: ")
        sel_carbon = str(u_set_act.select_atoms("bynum %s" % carbon))
        locA = sel_carbon.find('[<') + 2
        locB = sel_carbon.find(' and segid')
        break
    except SelectionError:
        print('The selection does not exists. Please, type an atom that exists.')
        continue


while True:
    try :
        quest = input('You selected this atom: %s. Is it correct ([y]/n)? ' % sel_carbon[locA:locB])

        if quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
            while True:
                try:
                    carbon = input("Type the number of the central carbon: ")
                    sel_carbon = str(u_set_act.select_atoms("bynum %s" % carbon))
                    locA = sel_carbon.find('[<') + 2
                    locB = sel_carbon.find(' and segid')
                    break
                except SelectionError:
                    continue
            continue

        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            break

    except ValueError:
        print("Type just \'yes\' or \'no\'.")

while True:
    try :
        radius_set_act = float(input("Which is the desired radius (in Å)? "))
        if radius_set_act < radius:
            break
        elif radius_set_act >= radius:
            print('The selected radius for the set_act list is smaller than the radius of solvent. Please, choose a new radius for the active atoms.')
            continue
    except ValueError:
        print("Type just the number, please.")
        continue

selection = u_set_act.select_atoms(str('byres around %s bynum %s' % (radius_set_act, carbon)))

txt = open('set_act_%s_%s' % (name, radius_set_act), 'w')
txt.write('set act { ')
for i in range(0, len(selection)):
    locA = str(selection[i]).find('<Atom ') + 6
    locB = str(selection[i]).find(': ')
    txt.write(str(selection[i])[locA:locB] + " ")
txt.write("} ")
txt.close()

print('%s atoms have been set as active for the QM/MM calculation using ChemShell.' % (len(selection)))


