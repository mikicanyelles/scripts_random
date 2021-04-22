"""
DESCRIPTION
    Script for fixing the numbering of the models in multimodel pdbs. If several pdbs are introduced, they will be concatenated.
"""


from sys import argv

out_file = open(str(argv[1]).split('.')[0] + '.out.pdb', 'w')

model = 1
for a in argv[1:]:
    file = open(a).readlines()
    for line in file:
        if line[:5] == 'MODEL':
            out_file.write(f'MODEL {model}\n')
            model += 1
        else :
            out_file.write(line)

out_file.close()




