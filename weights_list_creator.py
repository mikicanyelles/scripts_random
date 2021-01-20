from os import remove, path
from argparse import ArgumentParser


parser = ArgumentParser(description="Simple script for creating lists of weights for some ChemShell calculations")

parser.add_argument(
    '-n', '--number_atoms',
    help='Integer number corresponding to the total number of atoms of the system',
    required=True
    )
parser.add_argument(
    '-act',
    help='tcl list of the active atoms. This will be the reference for building the list of weights.',
    required=True
    )
parser.add_argument(
    '-qm',
    help='tcl list of the QM atoms. This is the second list, not necessarily the QM list',
    required=False
    )
parser.add_argument(
    '-core',
    help='tcl list of the core atoms. This is thought to be a little list of greater weight than the QM atoms',
    required=False,
    )
parser.add_argument(
    '-non_act_w', '--non_act_weigth',
    help='Float value for the non active atoms weigth (the non QM nor core). The default value is 0.0',
    required=False,
    default=0.0
    )
parser.add_argument(
    '-act_w', '--act_weigth',
    help='Float value for the active atoms weigth (the non QM nor core). The default value is 0.0',
    required=False,
    default=0.0
    )
parser.add_argument(
    '-qm_w', '--qm_weigth',
    help='Float value for the QM atoms weigth. The default value is 1.0',
    required=False,
    default=1.0
    )
parser.add_argument(
    '-core_w', '--core_weigth',
    help='Float value for the core atoms weigth. The default value is 2.0',
    required=False,
    default=2.0
    )

parser.add_argument(
    '-of', '--output_file_name',
    help='Name of the output file name',
    required=False,
    )

argsdict = vars(parser.parse_args())

try :
    float(argsdict['act_weigth'])
    float(argsdict['qm_weigth'])
    float(argsdict['core_weigth'])
except ValueError:
    print('The specified weights are not float numbers. Rerun the scripts specifying the weights in the proper way.')
    exit()



def list_from_file(filename):
    file_ = open(filename).readlines()#[0].split()
    list_ = []

    for i in range(len(file_)):
        for j in range(len(file_[i].split())):
            try :
                list_.append(int(file_[i].split()[j]))
            except ValueError:
                pass

    return list_

list_act = list_from_file(argsdict['act'])

if argsdict['qm'] != None:
    list_qm = list_from_file(argsdict['qm'])

if argsdict['core'] != None:
    list_ts = list_from_file(argsdict['core'])


if argsdict['output_file_name'] == None:
    filename = ''
    if len(argsdict['act'].split('/')) ==  1:#and argsdict['act'].split('/')[0] == '':
        filename = 'weigths_' + str(argsdict['act'])
        
    elif len(argsdict['act'].split('/')) > 1:
        for i in range(0,len(argsdict['act'].split('/'))-1):
            filename = filename + str(argsdict['act'].split('/')[i]) + '/'
            

        filename = filename + 'weigths_' + str(argsdict['act'].split('/')[-1])

else:
    filename = argsdict['output_file_name']

try:
    f = open(filename)
    exists = True
    f.close()
except IOError:
    exists = False
    
if exists == True:
    while True:
        quest = input('\'%s\' existis. Do you want to remove it ([y]/no)? ' % filename)
        if quest.lower() in ('', 'y', 'yes', '1'):
            remove(filename)
            break
        elif quest.lower() in ('n', 'no', '0'):
            filename = input('Specify the new name for the output file: ')
            print('The new file will be saved as %s' % filename)
            break
        else :
            print('Please, answer \'yes\' or \'no\'.')

weights = open(filename, 'w')

weights.write('set weights {')
for i in range(1,int(argsdict['number_atoms'])+1):
    if argsdict['qm'] != None:
        if argsdict['core'] != None:
            if i in list_act and i not in list_qm and i not in list_ts:
                weights.write(' ' + str(argsdict['act_weigth']))
            elif i in list_qm and i not in list_ts:
                weights.write(' ' + str(argsdict['qm_weigth']))
            elif i in list_ts:
                weights.write(' ' + str(argsdict['core_weigth']))
            else:
                weights.write(' ' + str(argsdict['non_act_weigth']))
        elif argsdict['core'] == None:
            if i in list_act and i not in list_qm:
                weights.write(' ' + str(argsdict['act_weigth']))
            if i in list_qm:
                weights.write(' ' + str(argsdict['qm_weigth']))
            else:
                weights.write(' ' + str(argsdict['non_act_weigth']))

    if argsdict['qm'] == None:
        if argsdict['core'] != None:
            if i in list_act and i not in list_ts:
                weights.write(' ' + str(argsdict['act_weigth']))
            elif i in list_ts:
                weights.write(' ' + str(argsdict['core_weigth']))
            else:
                weights.write(' ' + str(argsdict['non_act_weigth']))
        elif argsdict['core'] == None:
            if i in list_act:
                weights.write(' ' + str(argsdict['act_weigth']))
            else:
                weights.write(' ' + str(argsdict['non_act_weigth']))

weights.write(' }')
weights.close()
