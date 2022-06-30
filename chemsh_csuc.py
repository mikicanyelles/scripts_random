from os import getcwd


template_sh = """#!/bin/bash
#SBATCH -J opt_rc_{frame}
#SBATCH -e {stdouterr}.err
#SBATCH -o {stdouterr}.log
#SBATCH --mail-user=miquel.canyelles@uab.cat
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH -p std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -t {timelimit}

. /home/uabqut17/soft/chemshell_turbomole.sh

#Enviornment variables
export PARNODES=`echo $SLURM_JOB_NODELIST | wc -w`
export PARA_ARCH=SMP

#Location input/output data files
DIR=$SLURM_SUBMIT_DIR
SOURCEDIR={sourcedir}
file={filename}

#Copy input data files to the scratch folder
#mkdir -p $SCRATCH/
cd $SCRATCH/
cp -f $DIR/* $SCRATCH/
cp -f $SOURCEDIR/* $SCRATCH/
cp -f $SOURCEDIR/{frame}/* $SCRATCH/

#Run chemshell script
chemsh ${filename}.chm > $DIR/${filename}.out

#Copy output files from scratch to user folder
cp $SCRATCH/*.pdb $DIR/structures/
cp $SCRATCH/*.c   $DIR/structures/
"""

template_chm_opt_min = """# Creating vars
global sys_name_id
set frame {frame}
set structure_type {structure_type}
set calculation_type {kind}
set sys_name_id ${{calculation_type}}_${{structure_type}}_${{frame}}

# Load source tcls
## if Picard
    #source /QFsoft/applic/CHEMSHELL/3.7/intel2017-serial/src/interface_turbomole/turbomole.tcl
    #source /QFsoft/applic/CHEMSHELL/3.7/intel2017-serial/src/interface_amber/parse_amber.tcl
## if CSUC
    source /home/uabqut17/soft/chemshell-3.7/src/interface_amber/parse_amber.tcl
    source /home/uabqut17/soft/chemshell-3.7/src/interface_turbomole/turbomole.tcl

# Load coordinates, parameters and pdb
{no_c}
set prmtop_file ${{frame}}.cropped.mod.prmtop
set inpcrd_file ${{frame}}.cropped.inpcrd
set pdb_file    ${{frame}}.cropped.pdb

# Loading lists of atoms
source ./${{frame}}_15.0.act_list
source ./qm_list

# Defining residues for HDLC
set pdbresidues [ pdb_to_res "$pdb_file" ]
set fil [open "pdbresidues" {{ RDWR CREAT TRUNC }}]
puts $fil "set pdbresidues [ list $pdbresidues ]"
close $fil

# Performing QM/MM calculation

dl-find coords=$c_file \\
        result=${{sys_name_id}}.c \\
        residues= $pdbresidues \\
        coordinates= hdlc \\
        optimizer= lbfgs \\
        active_atoms= $act \\
        maxene= 100000 \\
        tolerance= 0.00045 \\
        maxstep= 0.1 \\
        lbfgs_mem= 1000 \\
        list_option= none \\{microiterative}
        theory= hybrid : [ list \\
                coupling= shift \\
                debug= no \\
                qm_region= $qm_list \\
                qm_theory= turbomole : [ list \\
                        hamiltonian= b3lyp \\
                        scftype= uhf \\
                        basisspec= {{ {{ lanl2dz fe }} {{ 6-31g* {{ h c n o }} }} }} \\
                        maxcyc= 1000 \\
                        charge= 1 \\
                        mult= 6 ] \\
                mm_theory= dl_poly : [ list \\
                        conn= $c_file \\
                        list_option= none \\
                        debug= no \\
                        scale14= [ list [ expr 1 /1.2 ] 0.5 ] \\
                        amber_prmtop_file=$prmtop_file ] \\s
        ] \\

# Converting .c into .pdb
read_pdb  file=$pdb_file coords=dummy.coords
write_pdb file=./${{sys_name_id}}.pdb coords=./${{sys_name_id}}.c

# End script
times

exit


"""
# no_c
##if yes:
############## Creating .c from AMBER files
##############load_amber_coords inpcrd=$inpcrd_file prmtop=$prmtop_file coords=$c_file
##### set c_file ${{frame}}.cropped.c
##if no:
###### set c_file {starting_structure}


template_chm_scan_r2 = """# Creating vars
    global sys_name_id
    set frame {frame}
    set structure_type   {structure_type}
    set calculation_type {kind}
    set Hname            {Hname}

    set localdir {localdir}
    set sys_name_id ${{calculation_type}}_${{structure_type}}_${{frame}}

# Load source tcls
## if Picard
    #source /QFsoft/applic/CHEMSHELL/3.7/intel2017-serial/src/interface_turbomole/turbomole.tcl
    #source /QFsoft/applic/CHEMSHELL/3.7/intel2017-serial/src/interface_amber/parse_amber.tcl
## if CSUC
    source /home/uabqut17/soft/chemshell-3.7/src/interface_amber/parse_amber.tcl
    source /home/uabqut17/soft/chemshell-3.7/src/interface_turbomole/turbomole.tcl

# Load files
    set c_file      opt_${{structure_type}}_${{frame}}.c
    set pdb_file    opt_${{structure_type}}_${{frame}}.pdb
    set prmtop_file ${{frame}}.cropped.mod.prmtop
    #set inpcrd_file ${{frame}}.cropped.inpcrd

# Loading lists of atoms
    source ./${{frame}}_15.0.act_list
    source ./qm_list

# Defining residues for HDLC
    set pdbresidues [ pdb_to_res "$pdb_file" ]
    set fil [open "pdbresidues" {{ RDWR CREAT TRUNC }}]
    puts $fil "set pdbresidues [ list $pdbresidues ]"
    close $fil

# Create PES.plt file for storing the PES
    set PESdata "j, rc, E, deltaE, d_CH, d_HO, rc_m"
    set PESnode "PES.plt"
    set PEShome "${{localdir}}/PES.plt"
    ## PES in node
    set PES [open $PESnode {{ RDWR CREAT TRUNC }} ]
    puts $PES $PESdata
    close $PES
    ## PES in home
    set PES [ open $PEShome {{ RDWR CREAT TRUNC }} ]
    puts $PES $PESdata
    close $PES

# Set atom numbers
    {atoms}


# Create point 0 of the scan from .c file
    exec cp $c_file ${{sys_name_id}}_0.c

# Set initial coordinates for scan
    set rc_short [interatomic_distance coords=./${{sys_name_id}}_0.c i=${atoms_rx_A}  j=${atoms_rx_B}]
    set rc_long [interatomic_distance coords=./${{sys_name_id}}_0.c i=${atoms_rx_B} j=${atoms_rx_C}]
    set rc_i [ expr {{ $rc_short - $rc_long }} ]
    set rc_f [ expr {{ {final_rc} * 1.889726 }} ]
    set step [ expr {{ {step} * 1.889726 }} ]

# Set initial counter
    set k 0

# Start performing QM/MM calculation

for {{ set rc $rc_i }} {{ $rc < $rc_f }} {{ set rc [ expr $rc + $step ] }} {{
    # Update counter
    set j [ expr {{ $k + 1 }} ]

    # Reinitalise internal energy variable
    matrix dl-find.energy new volatile

    # Start QM/MM calculation for step
    dl-find coords=${{sys_name_id}}_${{k}}.c \\
            result=${{sys_name_id}}_${{j}}.c \\
            residues= $pdbresidues \\
            coordinates= hdlc \\
            optimizer= lbfgs \\
            active_atoms= $act \\
            maxene= 100000 \\
            tolerance= 0.001 \\
            maxstep= 0.1 \\
            lbfgs_mem= 1000 \\
            list_option= none \\{microiterative}
            restraints= [ list [ list bonddiff2 ${atoms_rx_A} ${atoms_rx_B} ${atoms_rx_B} ${atoms_rx_C} $rc 3.0 ] ] \\
            theory= hybrid : [ list \\
                    coupling= shift \\
                    debug= no \\
                    qm_region= $qm_list \\
                    qm_theory= turbomole : [ list \\
                            hamiltonian= b3lyp \\
                            scftype= uhf \\
                            basisspec= {{ {{ lanl2dz fe }} {{ 6-31g* {{ h c n o }} }} }} \\
                            maxcyc= 1000 \\
                            charge= 1 \\
                            mult= 6 ] \\
                    mm_theory= dl_poly : [ list \\
                            conn= $c_file \\
                            list_option= none \\
                            debug= no \\
                            scale14= [ list [ expr 1 /1.2 ] 0.5 ] \\
                            amber_prmtop_file=$prmtop_file ] \\
        ] \

    # Convert .c into .pdb
    read_pdb  file=${{pdb_file}} coords=dummy.coords
    write_pdb file=./${{sys_name_id}}_${{j}}.pdb coords=./${{sys_name_id}}_${{j}}.c

    puts stdout "\nThe new coordinates file is: ${{sys_name_id}}_${{j}}.c\n"
    # Calculate relative energy
    set final_energy [expr {{[ get_matrix_element matrix=dl-find.energy indices= {{ 0 0 }} format=%12.12f ] * 627.509 }}]

        # Save first energy
        if {{ $k == 0 }} {{
            set energy_0 $final_energy
        }}


    set delta [ expr {{ $final_energy - $energy_0 }} ]

    # Measure distances
    set rc_  [ expr {{ $rc * 0.5291772 }} ]
    set d_CH [ interatomic_distance coords=./${{sys_name_id}}_${{j}}.c i=${atoms_rx_A} j=${atoms_rx_B} ]
    set d_HO [ interatomic_distance coords=./${{sys_name_id}}_${{j}}.c i=${atoms_rx_B} j=${atoms_rx_C} ]
    set rc_m [ expr {{ $d_CH - $d_HO }} ]

    # Save data to PESs files
    set PES [open $PESnode {{ RDWR APPEND }} ]
    puts $PES [ format \"%s %-5.5f %-5.5f %-5.5f %-5.5f %-5.5f %-5.5f\"  $j $rc_ $final_energy $delta $d_CH $d_HO $rc_m ]
    close $PES
    set PES [open $PEShome {{ RDWR APPEND }} ]
    puts $PES [ format \"%s %-5.5f %-5.5f %-5.5f %-5.5f %-5.5f %-5.5f\"  $j $rc_ $final_energy $delta $d_CH $d_HO $rc_m ]
    close $PES

    # Delete energy object
    delete_object dl-find.energy

    # Copy files to localdir
#    exec cp -f PES.plt $localdir
    exec cp -f ${{sys_name_id}}_${{j}}.c ${{localdir}}/structures
    exec cp -f ${{sys_name_id}}_${{j}}.pdb ${{localdir}}/structures

    # Increase counter
    incr k

}}

# End script
times

exit


"""


values = {
    "frame" : '10106',
    "timelimit" : '5-00:00',
    "sourcedir" : '/here/',
    "filename" : 'testing',
    "stdouterr" : '%x.%n',
    "structure_type" : 'rc',
    "kind" : 'scan_for1',
    "no_c" : 'set c_file',
    "localdir" : getcwd(),
    "Hname" : 'h12a',
    "atoms" : 'set OH 10104\nset C12 10105\nset H12A 10106\n',
    "atoms_rx_A" : 'C12',
    "atoms_rx_B" : 'H12A',
    "atoms_rx_C" : 'OH',
    "final_rc" : '4.0',
    "step" : '0.1',
    "microiterative" : '',
}



with open('sh.txt', 'w') as f:
    f.write(template_sh.format(**values))
    f.close()

with open('opt.txt', 'w') as f:
    f.write(template_chm_opt_min.format(**values))
    f.close()

with open('scan.txt', 'w') as f:
    f.write(template_chm_scan_r2.format(**values))
    f.close()