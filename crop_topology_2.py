# Import packages
from parmed import load_file
from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError
import sys
from argparse import ArgumentParser

parser = ArgumentParser(description="Simple script based on parmed and MDAnalysis for cropping AMBER topologies for ChemShell QM/MM calculations.")

parser.add_argument(
    '-p', '--parameters',
    help='File containing the parameters and topology of the system in the AMBER format. Accepted file formats: prmtop.\nIt needs the coordinates file.',
    required=False
    )
parser.add_argument(
    '-c', '--coordinates',
    help='File containing the coordinates of the system in AMBER format. Accepted file formats: inpcrd, rst7.\nIt needs the paramters file.',
    required=False
    )
parser.add_argument(
    '-pdb',
    help='File containing the coordinates of the system in the PDB format. If parameters are not given, ionly the generation of the active atoms list will be available.',
    required=False
    )
parser.add_argument(
    '-o', '--output',
    help='Name for the output file (the cropped parameters and coordinates).',
    required=False,
    )

argsdict = vars(parser.parse_args())

#if argsdict['parameters'] != None:
#    parameters = argsdict['parameters']

if argsdict['coordinates'] != None and argsdict['pdb'] != None:
    print('Both coordinates and pdb have been inputed. Only coordinates from the pdb will be used.')
    argsdict['coordinates'] = None
#    coordinates = argsdict['pdb']

#if argsdict['coordinates'] != None and argsdict['pdb'] == None:
#    coordinates = argsdict['coordinates']

#if argsdict['coordinates'] == None and argsdict['pdb'] != None:
#    coordinates = argsdict['pdb']



def options_parser(argsdict=argsdict):
    '''
    This function takes a the dictionary created using the argparse module
    and returns a dictionary with the available options. The argsdict is 
    also output because it may has suffered some modifications.
    The input dictionary has to include the following keys: 'parameters', 
    'coordinates', 'pdb' and 'output'. 
    The output dictionary has to include the following keys: 'crop topology', 
    'active list', 'parameters adaption'.
    '''

    options      = {'crop parameters' : None, 'active list' : None, 'parameters adaption' : None}
    options_done = {'crop parameters' : None, 'active list' : None, 'parameters adaption' : None}

    if argsdict['parameters'] != None and (argsdict['coordinates'] != None or argsdict['pdb'] != None):
        options['crop parameters']     = True
        options['active list']         = True
        options['parameters adaption'] = True
    
    elif argsdict['parameters'] != None and argsdict['coordinates'] == None and argsdict['pdb'] == None:
        options['crop parameters']     = False
        options['active list']         = False
        options['parameters adaption'] = True

    elif (argsdict['parameters'] == None and argsdict['coordinates'] == None) and argsdict['pdb'] != None:
        options['crop parameters']     = False
        options['active list']         = True
        options['parameters adaption'] = False

    elif argsdict['parameters'] == None and argsdict['coordinates'] == None and argsdict['pdb'] == None:
        print('No filename has been input, so there is nothing to do.')
        sys.exit(0)


    return options, options_done
    

def crop_top(argsdict=argsdict):

    if argsdict['parameters'] != None:
        parameters = argsdict['parameters']

    if argsdict['coordinates'] != None:
        coordinates = argsdict['coordinates']
    if argsdict['pdb'] != None:
        coordinates = argsdict['pdb']

    if argsdict['output'] != None:
        output = argsdict['output']
    elif argsdict['output'] == None:
        output = parameters[:-7]

    u_top = Universe(parameters)

    ligand = input('Type the number or the three letters code (only if it is a non-standard residue) of the central residue: ')
    try:
        if type(ligand) is int:
            sel = str(u_top.select_atoms('resid %s' % ligand).residues)
        if type(ligand) is str:
            sel = str(u_top.select_atoms('resname %s' % ligand).residues)
            print('The selected residue is %s' % sel[24:-3])
            err = False

    except SelectionError:
        print('The selected residue does not exist. Please, reselect it.')
        err = True
    
    while True:
        if err == True:
            quest = 'no'
        else :
            quest = input('Is the residue correct ([y]/n)? ')
        
        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
            break
        elif quest in ('no', 'NO', 'No', 'nO', 0):
            if err != True:
                print('Please, reselect the ligand then.')
            
            err = False
            ligand = input('Type the number or the three letters code (only if it is a non-standard residue) of the central residue: ')
            try:
                if type(ligand) is int:
                    sel = str(u_top.select_atoms('resid %s' % ligand).residues)
                if type(ligand) is str:
                    sel = str(u_top.select_atoms('resname %s' % ligand).residues)
                    print('The selected residue is %s' % sel[24:-3])

            except SelectionError:
                print('The selected residue does not exist. Please, reselect it.')
                err = True
            
            continue
        else :
            print('Type only yes or no')
            continue


    err = False
    radius = input('Type the desired radius around the selected ligand for the water drop (in Å): ')
    try:
        float(radius)
        print('The selected radius is %s Å' % radius)

    except ValueError:
        print('Type a number')
        err = True

    while True:
        if err == True:
            quest = 'no'
        else :
            quest = input('Is this correct ([y]/n)? ')
        
        if quest in ('', 'y', 'Y', 'yes', 'YES', 'Yes', 'yES', 'YEs', 'yeS', 'yEs', 'YeS', 1):
            break
        elif ('no', 'NO', 'No', 'nO', 0):
            print('Type the number again.')

            radius = input('Type the desired radius around the selected ligand for the water drop (in Å): ')
            try:
                float(radius)
                print('The selected radius is %s Å' % radius)
                err = False

            except ValueError:
                print('Type a number')
                err = True
            continue
        else :
            print('I did not understand the answer. Please, answer again.')
            continue

    del u_top, sel, quest, err

    topology = load_file(parameters, coordinates)

    topology.box = None
    topology.strip(':WAT&:%s<@%s' % (ligand, radius))
    topology.strip(':Na+&:%s<@%s' % (ligand, radius))
    topology.strip(':Cl-&:%s<@%s' % (ligand, radius))  

    topology.write_parm(output + '.cropped.prmtop')
    topology.save(output + '.cropped.inpcrd')

    print('Cropped topology and coordinates have been saved as \'%s\' and \'%s\'' % (output + '.cropped.prmtop', output + '.cropped.prmtop'))

    return radius

radius = 10000
def active_atoms_list(argsdict=argsdict,radius=radius):
    if argsdict['pdb'] != None:
        u_set_act = Universe(argsdict['pdb'])
    elif argsdict['parameters'] != None and argsdict['coordinates'] != None:
        u_set_act = Universe(argsdict['parameters'], argsdict['coordinates'])

    if argsdict['output'] != None:
        output = str(argsdict['output']) + '.act_list'
    elif argsdict['output'] == None:
        output = str(argsdict['parameters'])[:-7] + '.act_list'

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
            if radius_set_act < float(radius):
                break
            elif radius_set_act >= float(radius):
                print('The selected radius for the set_act list is smaller than the radius of solvent. Please, choose a new radius for the active atoms.')
                continue
        except ValueError:
            print("Type just the number, please.")
            continue
    
    selection = u_set_act.select_atoms(str('byres around %s bynum %s' % (radius_set_act, carbon)))

    txt = open('set_act_%s_%s' % (output, radius_set_act), 'w')
    txt.write('set act { ')
    for i in range(0, len(selection)):
        locA = str(selection[i]).find('<Atom ') + 6
        locB = str(selection[i]).find(': ')
        txt.write(str(selection[i])[locA:locB] + " ")
    txt.write("} ")
    txt.close()

    print('%s atoms have been set as active for the QM/MM calculation using ChemShell.' % (len(selection)))


def topology_adapter(argsdict=argsdict):

    top    = open(argsdict['parameters'], 'r').readlines()
    top_out = open(str(argsdict['parameters'])[:-7] + '.mod.prmtop', 'w')

    heterotypes_in = []
    metals_in      = []
    for i in range(0, len(top)):
        if '%FLAG AMBER_ATOM_TYPE' in top[i]:
            initial = i
        if '%FLAG TREE_CHAIN_CLASSIFICATION' in top[i]:
            final   = i-1
        if 'Y' in top[i] and '%' not in top[i]:
            loc = top[i].find('Y')
            heterotypes_in.append(str(top[i])[loc:loc+3])
        if 'M' in top[i] and '%' not in top[i]:
            loc = top[i].find('M')
            metals_in.append(str(top[i])[loc:loc+3])

    heterotypes_out = []
    for i in range(len(heterotypes_in)):
        ht_ = input('Which atom is %' % heterotypes_in[i])
        heterotypes_out.append('{:<3}'.format(ht_))
    
    metals_out = []
    for i in range(len(metals_in)):
        m_ = input('Which atom is %' % metals_in[i])
        metals_out.append('{:<3}'.format(m_))


    for line in top:
        l_ = line.replace('CO', 'C ')
        l_ = l_.replace('CX', 'C ')
        l_ = l_.replace('c ', 'C ')
        l_ = l_.replace('c2', 'C ')
        l_ = l_.replace('c3', 'C ')
        l_ = l_.replace('cx', 'C ')
        l_ = l_.replace('ce', 'C ')
        l_ = l_.replace('cf', 'C ')

        l_ = l_.replace('o ', 'O ')
        l_ = l_.replace('op', 'O ')
        l_ = l_.replace('os', 'O ')

        l_ = l_.replace('hc', 'H ')
        l_ = l_.replace('ha', 'H ')
        l_ = l_.replace('h1', 'H ')

        for j in range(len(heterotypes_out)):
            l_ = l_.replace(heterotypes_in[j], heterotypes_out[j])
        
        for j in range(len(metals_out)):
            l_ = l_.replace(metals_in[j], metals_out[j])

        top_out.write(l_)

    top_out.close()
    top.close()



options, options_done = options_parser()

if options['crop parameters'] == True:
    
    while True:    
        quest = input('Do you want to crop the system ([y]/n)? ')
        
        if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            radius = crop_top()
            options_done['crop parameters'] = True
            break
        elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
            break
        else :
            print('Type just \'yes\' or \'no\'.')
            continue
        

if options['active list'] == True:
    while True:    
        quest = input('Do you want to create the tcl list of the active atoms ([y]/n)? ')
        
        if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            active_atoms_list(argsdict, radius)
            options_done['active list'] = True
            break
        elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
            break
        else :
            print('Type just \'yes\' or \'no\'.')
            continue


if options['parameters adaption'] == True:
    while True:    
        quest = input('Do you want to adapt the atom types of the parameters file to ChemShell ([y]/n)? ')
        
        if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            while options_done['crop parameters'] == True:
                quest2 = input('Do you want to use the new topology (1) or the source one (2) ([1]/2)? ')

                if quest2 in ('', '1'):
                    if argsdict['output'] == None:
                        argsdict['parameters'] = str(argsdict['parameters'])[:-7] + '.cropped.prmtop'
                    elif argsdict['output'] != None:
                        argsdict['parameters'] = argsdict['output'] + '.cropped.prmtop'

                    topology_adapter()
                    options_done['crop parameters'] == False
                elif quest2 == '2':
                    topology_adapter()
                    options_done['crop parameters'] == False

                else :
                    print('Type just \'1\' or \'2\'.')
                    options_done['crop parameters'] == True
            break
        elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
            break
        else :
            print('Type just \'yes\' or \'no\'.')
            continue

