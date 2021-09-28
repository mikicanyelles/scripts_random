import numpy as np
import argparse 



def pdb_to_pdbdict(pdbfilename):
    """
    DESCRIPTION:
        Function for converting pdb file into dict, parsing the data of each atom 
        to a list of dictionaries for easier access.
    """
    
    
    
    pdbfile = open(pdbfilename, 'r').readlines()
    pdbdict = []
    for line in pdbfile:
        if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
            pdbdict.append({
                'record'    : line[0:6].rstrip(),
                'at_num'    : int(line[6:11]),
                'at_name'   : line[11:16],
                'alt_loc'   : line[16],
                'res_name'  : line[17:20],
                'chain'     : line[21],
                'res_num'   : int(line[22:26]),
                'ins_res'   : line[26],
                'x'         : float(line[30:38]),
                'y'         : float(line[38:46]),
                'z'         : float(line[46:54]),
                'occupancy' : line[54:60],
                'temp_fact' : line[60:66],
                'seg_id'    : line[72:76],
                'element'   : line[76:78],
            })
            
        elif line[0:3] == 'TER':
            pdbdict.append({
                'record'    : line[0:6],
                'at_num'    : int(line[6:11]),
                'res_name'  : line[17:20],
                'chain'     : line[21],
                'res_num'   : int(line[22:26]),
                'ins_res'   : line[26],
            })
            
            
    return pdbdict


def calculate_cog(pdbdict, selection):
    """
    DESCRIPTION:
        Function for calculating the center of geometry (cog) of a given selection.
    """
    
    points_x, points_y, points_z   = [], [], []
    
    for atom in pdbdict:
        if (atom['record'] == 'ATOM' or atom['record'] == 'HETATM') and atom['res_num'] in selection:
            points_x.append(atom['x'])
            points_y.append(atom['y'])
            points_z.append(atom['z'])
            
    
    x = round(sum(points_x)/len(points_x), 3)
    y = round(sum(points_y)/len(points_y), 3)
    z = round(sum(points_z)/len(points_z), 3)
    
    return [x,y,z]
    
def pdbdict_move_to_coc(pdbdict, cog, selection):
    """
    DESCRIPTION
        Function for moving the given selection to the center of coordinates (point 0,0,0) from its center of geometry.
    """
    
    for atom in pdbdict:
        if atom['res_num'] in selection:
            try : 
                atom['x'], atom['y'], atom['z'] = round(atom['x'] - cog[0], 3), round(atom['y'] - cog[1], 3), round(atom['z'] - cog[2], 3)
            except KeyError:
                pass
        
    return pdbdict

    
def pdbdict_move_to_cog(pdbdict, cog, selection):
    """
    DESCRIPTION
        Function for moving the given selection to the center of geometry from its center of coordinates (point 0,0,0).
    """
    
    for atom in pdbdict:
        if atom['res_num'] in selection:
            try : 
                atom['x'], atom['y'], atom['z'] = round(atom['x'] + cog[0], 3), round(atom['y'] + cog[1], 3), round(atom['z'] + cog[2], 3)
            except KeyError:
                pass
        
    return pdbdict


def pdbdict_to_pdb(pdbfilename, pdbdict):#, selection):
    """
    DESCRIPTION
        Function for generating a pdb file from a pdbdict.
    """
    
    #pdbfilein  = open(pdbfilename).readlines()
    pdbfileout = open(''.join(pdbfilename.split('.')[:-1]) + '.out.' + str(pdbfilename.split('.')[-1]), 'w')
    
    for atom in pdbdict:
        
        if atom['record'] == 'ATOM' or atom['record'] == 'HETATM':
            if len(atom['at_name']) == 3:
                at_name = atom['at_name']
            
            line =  '{:6}'.format(atom['record'])     + \
                    '{:>5}'.format(atom['at_num'])    + \
                    '{:4}'.format(atom['at_name'])    + \
                    '{:1}'.format(atom['alt_loc'])    + \
                    '{:>3}'.format(atom['res_name'])  + \
                    ' '                               + \
                    '{:1}'.format(atom['chain'])      + \
                    '{:4}'.format(atom['res_num'])    + \
                    '{:1}'.format(atom['ins_res'])    + \
                    '   '                             + \
                    '{:>8}'.format(atom['x'])         + \
                    '{:>8}'.format(atom['y'])         + \
                    '{:>8}'.format(atom['z'])         + \
                    '{:6}'.format(atom['occupancy'])  + \
                    '{:6}'.format(atom['temp_fact'])  + \
                    '      '                          + \
                    '{:4}'.format(atom['seg_id'])     + \
                    '{:2}'.format(atom['element']) + '\n'
                    
        
        if atom['record'] == 'TER':
            line =  '{:6}'.format(atom['record'])     + \
                    '{:>5}'.format(atom['at_num'])    + \
                    '{:4}'.format(atom['at_name'])    + \
                    '{:1}'.format(atom['alt_loc'])    + \
                    '{:>4}'.format(atom['res_name'])  + \
                    '{:1}'.format(atom['chain'])      + \
                    ' '                               + \
                    '{:>4}'.format(atom['res_num'])   + \
                    '{:1}'.format(atom['ins_res']) + '\n'
    
        
        pdbfileout.write(line)

    
    pdbfileout.close()
        


def rotate_over_axis(pdbdict, selection, angle, axis='z'):
    """
    DESCRIPTION
        Function for rotating a given selection with a given angle over 
        the selected axis (x, y or z (z is the default)).
    """

    angle = np.deg2rad(angle)
    
    if axis == 'x':
    
        rotation_matrix = np.array(
            [
                [1, 0,             0],
                [0, np.cos(angle), -np.sin(angle)],
                [0, np.sin(angle), np.cos(angle)]
            ]
        )
    
    elif axis == 'y':
        rotation_matrix = np.array(
            [
                [np.cos(angle),  0, np.sin(angle)],
                [0,              1, 0],
                [-np.sin(angle), 0, np.cos(angle)]
            ]
        )
    
    elif axis == 'z':
        rotation_matrix = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle),  0],
                [0,             0,              1]
            ]
        )

    
    for atom in pdbdict:
        if atom['res_num'] in selection:
            try : 
                v = np.dot(rotation_matrix, [atom['x'], atom['y'], atom['z']])
                
                atom['x'] = round(float(v[0]), 3)
                atom['y'] = round(float(v[1]), 3)
                atom['z'] = round(float(v[2]), 3)
                
            except KeyError:
                pass
        
    return pdbdict


def parse_selection(selection):

    selection_list = []
    for sel in selection:
        if sel.find('-') != -1:
            selection_list = selection_list + list(range(int(sel.split('-')[0]), int(sel.split('-')[-1])+1))

        else :
            selection_list.append(int(sel))

    return selection_list


def main():

    parser = argparse.ArgumentParser(description='Program for moving and rotating a structure from a PDB')

    parser.add_argument('-a', '--angle', type=float, nargs=1, help='Angle of rotation')
    parser.add_argument('-ax', '--axis', type=str, nargs=1, help='Axis of rotation', choices=['x', 'y', 'z'], default='z')


    parser.add_argument('-o', '--output', type=str, nargs=1, help='Name of the generated PDB file')

    parser.add_argument('-p', '--pdb_file', type=str, nargs=1, required=True, help='Input a pdb structure.')

    parser.add_argument('-s', '--selection', type=str, nargs='+', help='Selection of residues (residue number) for applying the rotation. They can be listed one by one and/or using a dash for selection of several consecutive residues.')

    parser.add_argument('-m', '--movement', type=list, nargs=1, help='Vector (x,y,z) for translating the structure.')


    args = vars(parser.parse_args())
    

    pdbdict = pdb_to_pdbdict(args['pdb_file'][0])

    if args['output'] == None:
        args['output'] = ''.join(args['pdb_file'][0].split('.')[:-1]) + '.out.' + str(args['pdb_file'][0].split('.')[-1])

    if args['selection'] != None:
        selection = parse_selection(args['selection'])
    
    elif args['selection'] == None:
        selection = list(range(pdbdict[0]['res_num'], pdbdict[-1]['res_num'] + 1))

    cog = calculate_cog(pdbdict, selection)

    # only rotation
    if args['movement'] == None and args['angle'] != None:

        pdbdict = pdbdict_move_to_cog(
            rotate_over_axis(
                pdbdict_move_to_coc(
                    pdbdict,
                    cog,
                    selection
                ),
            selection,
            args['angle'],
            args['axis']
            ),
            cog,
            selection
        )

    # only movement (atom_coords - movement_vector)
    if args['movement'] != None and args['angle'] == None:
        pdbdict = pdbdict_move_to_coc(
            pdbdict,
            args['movement'],
            selection
        )


    # rotation and movement
    if args['movement'] != None and args['angle'] != None:
        pdbdict = pdbdict_move_to_coc(
            pdbdict_move_to_cog(
                rotate_over_axis(
                    pdbdict_move_to_coc(
                        pdbdict,
                        cog,
                        selection
                    ),
                selection,
                args['angle'],
                args['axis']
                ),
                cog,
                selection
            ),
            args['movement'],
            selection
        )
    

    pdbdict_to_pdb(args['output'], pdbdict)


if __name__ == '__main__':
    main()