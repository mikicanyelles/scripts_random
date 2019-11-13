# README #

Sort of different scripts written mainly in Python or Shell which are useful for QM/MM calculations using ChemShell, Turbomole and Amber.

# atoms_list_to_pdb.py

Takes a list of atoms (like the ones from ChemShell scripts) and the complete pdb and generates a new pdb which includes only the atoms in the list.


# crop_topology.py

Takes an AMBER topology file and the coordinates (either inpcrd or pdb or mol2) and crops the water box removing all atoms outside a drop of a desired radius drawn around the selected residue. It also generates a list of atoms at a selected radius from an specified atom.


# energy_plots.py

Takes the output of a QM/MM calculation from ChemShell (interfaced to either Turbomole or Gaussian) and plots the energies (QM/MM, QM and MM) along the optimisation of the structure or the optimisation of each point of a scan calculation


# gantt_phd.py

Creates a Gantt Diagram using the python-gantt package.


# taglatex.py

Parser for tags in LaTeX source files. It is able to read either a tex source file or a dict file (which is a plain text with a python dictionary in one line) and write an itemized line with the 'Table of Tags' chapter in a LaTeX source file or a dict file.

---

# To-Do list
- [ ] Add saving the pngs as an option (which can be modified modifing the script)
