import MDAnalysis as mda
import numpy as np
import warnings
from math import sqrt


pdbfile = 'opt_pd_0940_DHA-12LOX.pdb'

def create_o2_universe():
    n_residues = 1
    n_atoms = n_residues * 2

    resindices = np.repeat(range(n_residues), 2)
    #assert len(resindices) == n_atoms
    #print("resindices:", resindices[:10])

    segindices = [0] * n_residues
    #print("segindices:", segindices[:10])

    oxy = mda.Universe.empty(n_atoms,
                            n_residues=n_residues,
                            atom_resindex=resindices,
                            residue_segindex=segindices,
                            trajectory=True)
    #print(oxy)

    oxy.add_TopologyAttr('name', ['O1', 'O2']*n_residues)
    #print(oxy.atoms.names)

    oxy.add_TopologyAttr('type', ['O', 'O']*n_residues)
    #print(oxy.atoms.types)

    oxy.add_TopologyAttr('resname', ['OXY']*n_residues)
    oxy.add_TopologyAttr('resid', [1])
    #print(oxy.atoms.resids)

    #print(oxy.atoms.positions)
    print(oxy)

    return oxy

def load_protein(pdbfilename):
    
    u = mda.Universe(pdbfilename)

    #solvent = u.select_atoms('resname WAT or resname Na+ or resname NA+ or resname Cl- or resname CL-')

    #if len(solvent.residues) != 0:
    #    first_solvent = int(str(solvent.residues[0])[(str(solvent.residues[0]).find(', ')+2):str(solvent.residues[0]).find('>')]) - 1
    #    protein, solvent = u.residues[:first_solvent], u.residues[first_solvent:]

    #else :
    #    protein, solvent = u.residues[:], None


    #return protein.atoms, solvent.atoms
    return u


positions = np.array([
    [1,0,0],
    [-1,0,0],    
    [0,1,0],
    [0,-1,0],
    [0,0,1],
    [0,0,-1],
    
    [(1/sqrt(2)),(1/sqrt(2)),0],
    [(1/sqrt(2)),0,(1/sqrt(2))],
    [0,(1/sqrt(2)),(1/sqrt(2))],
    [-(1/sqrt(2)),(1/sqrt(2)),0],
    [-(1/sqrt(2)),0,(1/sqrt(2))],
    [0,-(1/sqrt(2)),(1/sqrt(2))],
    [(1/sqrt(2)),-(1/sqrt(2)),0],
    [(1/sqrt(2)),0,-(1/sqrt(2))],
    [0,(1/sqrt(2)),-(1/sqrt(2))],
    [-(1/sqrt(2)),-(1/sqrt(2)),0],
    [-(1/sqrt(2)),0,-(1/sqrt(2))],
    [0,-(1/sqrt(2)),-(1/sqrt(2))],
    
    [(1/sqrt(3)),(1/sqrt(3)),(1/sqrt(3))],
    [-(1/sqrt(3)),(1/sqrt(3)),(1/sqrt(3))],
    [(1/sqrt(3)),-(1/sqrt(3)),(1/sqrt(3))],
    [(1/sqrt(3)),(1/sqrt(3)),-(1/sqrt(3))],
    [-(1/sqrt(3)),-(1/sqrt(3)),(1/sqrt(3))],
    [-(1/sqrt(3)),(1/sqrt(3)),-(1/sqrt(3))],
    [(1/sqrt(3)),-(1/sqrt(3)),-(1/sqrt(3))],
    [-(1/sqrt(3)),-(1/sqrt(3)),-(1/sqrt(3))]
])


def add_o2(protein, oxy):

#    print(protein)
#    print(oxy)
#    print(solvent)

    print(len(protein.residues))
    solvent = protein.select_atoms('resname WAT or resname Na+ or resname NA+ or resname Cl- or resname CL-')

    if len(solvent.residues) != 0:
        first_solvent = int(str(solvent.residues[0])[(str(solvent.residues[0]).find(', ')+2):str(solvent.residues[0]).find('>')])
        last_solvent  = int(str(solvent.residues[-1])[(str(solvent.residues[-1]).find(', ')+2):str(solvent.residues[-1]).find('>')])

        print(first_solvent, last_solvent)
    
        protein.residues[first_solvent:last_solvent].resids = list(range(first_solvent +1, last_solvent+1))
        oxy.residues.resids     = first_solvent

        u = mda.Merge(protein.residues[:first_solvent-1].atoms, oxy.residues[:].atoms, protein.residues[first_solvent-1:].atoms)
    
    elif solvent == None:
        oxy.residues.resids     = list(range(len(protein.residues + 1)))

        u = mda.Merge(protein, oxy)


    u.atoms.write('test.pdb')


    return u






#sol = create_o2_universe()

#sol.atoms.positions = np.array([[41.6030, 34.3650, 32.4820], [40.4570, 34.8090, 32.4400]])
#print(sol.atoms.positions)

#warnings.filterwarnings("ignore")
#u.atoms.write('test.pdb')
#warnings.filterwarnings("default")

u = load_protein(pdbfile)
oxy = create_o2_universe()

