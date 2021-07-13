#from _typeshed import Self
#from o2_addition import create_o2_universe, merge_protein_o2
import MDAnalysis as mda
from MDAnalysis.lib import distances as mdadist
#from MDAnalysis.core.groups import Residue
import numpy as np
import warnings
from math import dist, sqrt
import os


pdbfile = 'opt_pd_0940_DHA-12LOX.pdb'


class AddOxygen:

    def __init__(self, pdbfilename, at_num, distance=3.0, resolution='high', savepdb=True, prefix='', hide_saving_warnings=True):
        """
        DESCRIPTION:
            __init__ method. Loads protein, generates coc and saves distance and resolutions parameters.
        """

        self.protein       = mda.Universe(pdbfilename)
        self.coc           = self.protein.select_atoms('bynum %s' % at_num).positions[0]
        self.distance      = distance
        self.resolution    = resolution
        self.savepdb      = savepdb
        self.prefix        = prefix
        self.hide_warnings = hide_saving_warnings

        self.positions = {
            'high' : np.array([
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
            ]),

            'medium' : np.array([
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
            ]),

            'low' : np.array([
                [1,0,0],
                [-1,0,0],    
                [0,1,0],
                [0,-1,0],
                [0,0,1],
                [0,0,-1],
            ])
        }


    def create_o2_universe(self):
        n_residues = 1
        n_atoms = n_residues * 2

        resindices = np.repeat(range(n_residues), 2)

        segindices = [0] * n_residues

        oxy = mda.Universe.empty(n_atoms,
                                n_residues=n_residues,
                                atom_resindex=resindices,
                                residue_segindex=segindices,
                                trajectory=True)

        oxy.add_TopologyAttr('name', ['O1', 'O2']*n_residues)
        oxy.add_TopologyAttr('type', ['O', 'O']*n_residues)
        oxy.add_TopologyAttr('resname', ['OXY']*n_residues)
        oxy.add_TopologyAttr('resid', [1])


        return oxy


    def build_positions(self, oxy, distance=3.0, resolution='high'):
        """
        DESCRIPTION:

        """

        if distance != 3.0:
            self.distance = distance
        if resolution != 'high':
            self.resolution = resolution

        o2_interatomic = 1.209152596

        if resolution.lower() in ('high', 'h'):
            positions = self.positions['high']
        elif resolution.lower() in ('medium', 'm'):
            positions = self.positions['medium']
        elif resolution.lower() in ('low', 'l'):
            positions = self.positions['low']

        oxys = []

        for p in range(len(positions)):
            #o = self.create_o2_universe()
            oxys.append(oxy.copy())
            

            oxys[p].atoms.positions = np.array(
                [positions[p]*(distance)+self.coc,
                positions[p]*(distance+o2_interatomic)+self.coc])

        return oxys



    def merge_protein_oxy(self, oxys, savepdb=True, prefix=''):

        if prefix != '':
            self.prefix = prefix
        
        if self.prefix == '':

            if 'protein_oxy' not in os.listdir():
                os.mkdir('protein_oxy')

            self.prefix = 'protein_oxy/' + pdbfile.split('.')[0]


        if savepdb != True:
            self.savepdb = savepdb


        
        # Finding first solvent molec (if any)
        solvent = self.protein.select_atoms('resname WAT or resname Na+ or resname NA+ or resname Cl- or resname CL-')

        if len(solvent.residues) != 0:
            first_solvent = int(str(solvent.residues[0])[(str(solvent.residues[0]).find(', ')+2):str(solvent.residues[0]).find('>')])
            last_solvent  = int(str(solvent.residues[-1])[(str(solvent.residues[-1]).find(', ')+2):str(solvent.residues[-1]).find('>')])

        else :
            first_solvent, last_solvent = None



        # Generating each universe
        #if merge_all == True:
        #    protein_all = self.protein.copy()
        #    protein_all.residues[first_solvent:last_solvent].resids = list(range(first_solvent + len(oxys), last_solvent+len(oxys)))

        # Modifying resids indexes
        if first_solvent != None and last_solvent != None:
            self.protein.residues[first_solvent:last_solvent].resids = list(range(first_solvent +1, last_solvent+1))

        proteins = []

       
        
        if isinstance(oxys, list):
            if self.hide_warnings == True:
                warnings.filterwarnings('ignore')
            for oxy in range(len(oxys)):
                if first_solvent != None and last_solvent != None:
                    oxys[oxy].residues.resids = first_solvent

                    proteins.append(
                        mda.Merge(self.protein.residues[:first_solvent-1].atoms, oxys[oxy].residues[:].atoms, self.protein.residues[first_solvent-1:].atoms)
                    )

                    if self.savepdb == True:
                        proteins[oxy].atoms.write('%s_OXY_%s.pdb' % (self.prefix,oxy+1))
        
            if self.hide_warnings == True:
                warnings.filterwarnings('default')

            return proteins

        else :
            if self.hide_warnings == True:
                warnings.filterwarnings('ignore')

            if first_solvent != None and last_solvent != None:
                oxys.residues.resids = first_solvent
            
            protein = mda.Merge(self.protein.residues[:first_solvent-1].atoms, oxys.residues[:].atoms, self.protein.residues[first_solvent-1:].atoms)


            if self.savepdb == True:
                protein.atoms.write('%s_OXY.pdb' % (self.prefix))
            
            if self.hide_warnings == True:
                warnings.filterwarnings('default')

            return protein
    
    
    def merge_protein_all_oxys(self, oxys, savepdb=True, prefix=''):

        if prefix != '':
            self.prefix = prefix
        
        if self.prefix == '':

            if 'protein_oxy' not in os.listdir():
                os.mkdir('protein_oxy')

            self.prefix = 'protein_oxy/' + pdbfile.split('.')[0]


        if savepdb != True:
            self.savepdb = savepdb


        
        # Finding first solvent molec (if any)
        solvent = self.protein.select_atoms('resname WAT or resname Na+ or resname NA+ or resname Cl- or resname CL-')

        if len(solvent.residues) != 0:
            first_solvent = int(str(solvent.residues[0])[(str(solvent.residues[0]).find(', ')+2):str(solvent.residues[0]).find('>')])
            last_solvent  = int(str(solvent.residues[-1])[(str(solvent.residues[-1]).find(', ')+2):str(solvent.residues[-1]).find('>')])

        else :
            first_solvent, last_solvent = None


        protein_all = self.protein.copy()
        #protein_all.residues[first_solvent:last_solvent].resids = list(range(first_solvent + len(oxys), last_solvent + len(oxys)))

        # Modifying resids indexes
        if first_solvent != None and last_solvent != None:
            protein_all.residues[first_solvent-1:last_solvent-1].resids = list(range(first_solvent + len(oxys), last_solvent + len(oxys)))
        
        if self.hide_warnings == True:
            warnings.filterwarnings('ignore')

        for oxy in range(len(oxys)):
            if first_solvent != None and last_solvent != None:
                oxys[oxy].residues.resids = first_solvent + oxy

                if oxy == 0:
                    protein = mda.Merge(protein_all.residues[:first_solvent-1].atoms, oxys[oxy].residues.atoms)
                elif oxy == (len(oxys) - 1):
                    protein = mda.Merge(protein.residues[:].atoms, oxys[oxy].residues.atoms, protein_all.residues[first_solvent-1:].atoms)
                else :
                    protein = mda.Merge(protein.residues[:].atoms, oxys[oxy].residues.atoms)

        if self.savepdb == True:
            protein.atoms.write('%s_ALL_OXY.pdb' % (self.prefix))
    
        if self.hide_warnings == True:
            warnings.filterwarnings('default')

        return protein


    def check_clashes(self, proteins, prefix=''):
        """
        DESCRIPTION:

        """

        if prefix != '':
            self.prefix = prefix

        f = open(self.prefix + '_clashes.txt', 'w')

        f.write('The following positions have a clash with an atom of the protein or the substrate.\n\n')
        print('The following positions have a clash with an atom of the protein or the substrate.\n\n')

        for protein in range(len(proteins)):

            sel_oxyenv  = proteins[protein].select_atoms('(around 4 resname OXY) and not resname OXY').positions
            sel_oxy     = proteins[protein].select_atoms('resname OXY').positions

            oxyenv_dist = np.min(mdadist.distance_array(sel_oxy, sel_oxyenv))

            if oxyenv_dist < 1.5:
                f.write('\tPosition number %s has a clash.' % protein)
                print('\tPosition number %s has a clash.' % protein)

        f.close()



    def run(self, saveall=True):

        #u             = AddOxygen()
        oxys          = self.build_positions(self.create_o2_universe())
        proteins      = self.merge_protein_oxy(oxys)
        self.check_clashes(proteins)
        if saveall == True:
            self.merge_protein_all_oxys(oxys)


                
        
