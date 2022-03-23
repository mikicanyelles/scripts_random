import re
from MDAnalysis import Universe
import MDAnalysis.analysis.rms as rms
import argparse
import  matplotlib.pyplot as plt
import re

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

parser = argparse.ArgumentParser(description='Script for performing RMSD of the CA for a given trajectory.')

parser.add_argument(
    '-p', '--parameters',
    help='Flag for specify the parameters in any format compatible with MDAnalysis (prmtop, for example).',
    required=True
)

parser.add_argument(
    '-t', '--trajectory',
    help='Flag for specify the trajectory or trajectories files.',
    required=True,
    nargs='+'
)

argsdict = vars(parser.parse_args())

if len(argsdict['trajectory']) > 1:
    argsdict['trajectory'] = natural_sort(argsdict['trajectory'])

u = Universe(argsdict['parameters'], argsdict['trajectory'])

sel = u.select_atoms('name CA')

R = rms.RMSD(
    sel,
    )
R.run()


rmsd = R.rmsd.T
time = rmsd[1]

plt.plot(time, rmsd[2])
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (Å)')
plt.show()
plt.close()
