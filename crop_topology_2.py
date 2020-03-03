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

if argsdict['coordinates'] != None and argsdict['pdb'] != None
    print('Both coordinates and pdb have been inputed. Only coordinates from the pdb will be used.')
    argsdict['coordinates'] = None

def options_parser(argsdict):
    '''
    This function takes a the dictionary created using the argparse module
    and returns a dictionary with the available options. The argsdict is 
    also output because it may has suffered some modifications.
    The input dictionary has to include the following keys: 'parameters', 
    'coordinates', 'pdb' and 'output'. 
    The output dictionary has to include the following keys: 'crop topology', 
    'active list', 'parameters adaption'.
    '''

    options = {'crop parameters' : None, 'active list' : None, 'parameters adaption' : None}
 
    if argsdict['parameters'] != None and (argsdict['coordinates'] != None or argsdict['pdb'] != None]:
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

    elif argsdict['parameters'] == None and argsdict['coordinates'] == None) and argsdict['pdb'] == None:    
        print('No filename has been input, so there is nothing to do.')
        sys.exit(0)


    return options
    




    # if len(argsdict) < 2:
    #     print('Input the topology/parameters and the coordinates as arguments, respectively. A name for the output file can be specified as the third argument (without an extension).')

    # elif len(argsdict) == 2:
    #     name_pdb = argsdict[1]
    #     only_list = True

    # elif len(argsdict) >= 3:
    #     name_prmtop = argsdict[1]
    #     only_list = False

    #     if str(argsdict[2]).find('.pdb') == (len(argsdict[2]) - 4):
    #         name_pdb = argsdict[2]
    #         name_inpcrd = None
    #     elif str(argsdict[2]).find('.inpcrd') == (len(argsdict[2]) - 7):
    #         name_inpcrd = argsdict[2]
    #         name_pdb = None