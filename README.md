# README #

# Description of the repo

Sort of different scripts written mainly in Python or Shell which are useful for QM/MM calculations using ChemShell, Turbomole and Amber.


# List of scripts

## atoms_list_to_pdb

Takes a list of atoms (like the ones from ChemShell scripts) and the complete pdb and generates a new pdb which includes only the atoms in the list.


## calcs_to_DB*

Reads the status of calculations running in HPC servers and dumps it into a YAML Database.


## charges_nbo_turbomole

Script for simplifying the analysis of NBO (by now) charges calculated using the dscf module from Turbomole, which can be obtained by adding the `$pop nba' keyword into the control file.


## crop_topology

Takes an AMBER topology file and the coordinates (either inpcrd or pdb or mol2) and crops the water box removing all atoms outside a drop of a desired radius drawn around the selected residue. It also generates a list of atoms at a selected radius from an specified atom.


## energy_plots

Takes the output of a QM/MM calculation from ChemShell (interfaced to either Turbomole or Gaussian) and plots the energies (QM/MM, QM and MM) along the optimisation of the structure or the optimisation of each point of a scan calculation


## hessian_parser (working on it...)

Takes the output of a hessian calculation of ChemShell (using the `force` module), parsers it and lets the user extract the most useful information such as the scructure of the QM zone (including the link atoms), a list of the calculated frequencies, the most contributing atom to a selected frequency and the contribution to a selected frequency of a selected atom. 


## multiframe

Script from creating a multimodel/multiframe pdb file from a list of input files.


## printPES

Script for printing a PES.plt file (generated with ChemShell for a PES exploration job) of a running job by inputting the job ID.


## taglatex

Parser for tags in LaTeX source files. It is able to read either a tex source file or a dict file (which is a plain text with a python dictionary in one line) and write an itemized line with the 'Table of Tags' chapter in a LaTeX source file or a dict file.


## weights_list_creator

Script for creating a list of weights from a list of atoms. Four types of atoms can be considered. The output is a file containing a tcl list. Useful for some kind of jobs using ChemShell (like TS optimisations using the dimer method).


---

# To-Do list

- Add example pics to description of each script which generates a plot.

### energy_plots

- Add saving the _pngs_ as an option (which can be modified modifing the script)
