#!/home/mcanyelles/miniconda3/envs/py_env/bin/python3
# coding: utf-8

from matplotlib import pyplot as plt
from os import listdir
from sys import argv, exit

while True:
    try :
        if argv[1] in listdir():
            file_name = argv[1]
            break
        elif argv[1] not in listdir():
            print("This file doesn't exist.")
            file_name = input('Which is the name of the .out file?')
            break
    except IndexError: #len(argv[1]) == 0 or len(argv[2]) == 0:
        #file_name = input('Which is the name of the .out file?')
        #break
        print('.out file is not specified. Please, specify it as an argument. (> energy_plotter file_name.out)')
        exit(0)


file = open('%s' % file_name, 'r')

# Checker for Gaussian/Turbomole, HDLCOpt/FL-Find, Opt/Scan
new = []

for line in file:
    if line.find('Contribution to energy from                  turbomole:') != -1:
        qm_prog = 'turbomole'
    elif line.find('Contribution to energy from                   gaussian:') != -1:
        qm_prog = 'gaussian'

    if line.find('Cycle      1,') != -1:
        new.append('new')
    
if len(new) == 1 or len(new) == 0:
    type = 'optimisation'
else :
    type = 'scan'

del file; del line
#####

# 
file = open('%s' % file_name, 'r')

if type == 'optimisation':
    
    if qm_prog == 'turbomole':
        opt_energy = [[], [], [], []]
        
        for line in file:
            if line.find('Cycle ') != -1:
                opt_energy[0].append(int(line[9:12]))
            elif line.find('QM/MM Energy: ') != -1:
                opt_energy[1].append(float(line[13:-7]))
            elif line.find('Contribution to energy from                  turbomole:') != -1:
                opt_energy[2].append(float(line[55:-7]))
            elif line.find('Contribution to energy from                    dl_poly:') != -1:
                opt_energy[3].append(float(line[55:-7]))
            else :
                pass

    elif qm_prog == 'gaussian':
        opt_energy = [[], [], [], []]
        
        for line in file:
            if line.find('Cycle ') != -1:
                opt_energy[0].append(int(line[9:12]))
            elif line.find('QM/MM Energy: ') != -1:
                opt_energy[1].append(int(line[13:-7]))
            elif line.find('Contribution to energy from                   gaussian:') != -1:
                opt_energy[2].append(int(line[55:-7]))
            elif line.find('Contribution to energy from                    dl_poly:') != -1:
                opt_energy[3].append(int(line[55:-7]))
            else :
                pass

    
    while True:
        try :
            opt_energy[1][len(opt_energy[0])-1]
            break
        except IndexError:
            del opt_energy[0][-1]
            break

    while True:
        try:
            if argv[2] == 'csv':
                csv = open('%s_energies.csv' % file_name, 'w+')
                csv.write("Cycle, QM/MM energy, QM energy, MM energy\n")
                for n in range(0,len(opt_energy[0])):
                    csv.write(str(opt_energy[0][n])+','+str(opt_energy[1][n])+','+str(opt_energy[2][n])+','+str(opt_energy[3][n])+'\n')
                csv.close()
                print("CSV file saved.")
                break
            else :
                break
        except IndexError :
            break

    plt.plot(opt_energy[0], opt_energy[1])
    plt.xlabel('Number of cycles')
    plt.ylabel('QM/MM Potential Energy (Hartree)')
    plt.legend(["QM/MM\nPotential\nEnergy"], bbox_to_anchor=(1.02,1), loc = 2)
    plt.grid(True)

    plt.title('QM/MM energies along the optimisation', y=1.05, loc='center')
    plt.savefig('QMMM_energies.png', transparent=False, dpi=300, bbox_inches='tight')
    plt.show()

    print("Plot for QM/MM energies saved")

    fig, ax1 = plt.subplots()
    ax1.plot(opt_energy[0], opt_energy[2], color='blue')#, marker='o')
    ax1.set_xlabel('Number of cycles')
    ax1.set_ylabel('QM Potential energy (Hartree)', color='blue')
    ax1.tick_params(axis='y', colors='blue')
    #ax1.yaxis.grid(linestyle='-', linewidth='0.5', color='blue')
    ax1.legend(["QM\nenergy"], bbox_to_anchor=(1.20,1), loc = 2)

    ax2=ax1.twinx()
    #ax2.plot(csv_ar[:,2], csv_ar[:,-1], color='red')
    ax2.plot(opt_energy[0], opt_energy[3], color='red')#, marker='^')
    #ax2.set_xlabel('O-H distance (Å)')
    #for k in range(1,i+1):
    #    ax2.plot(range(globals()['len_%s' % (k-1)]+1,globals()['len_%s' % k]+1), globals()['table_%s' % k][3], linestyle='--')
    ax2.set_ylabel('MM Potential energy (Hartree)', color='red')
    #ax2.set_yticks(mm_ticks)
    ax2.tick_params(axis='y', colors='red')
    ax2.legend(["MM\nenergy"], bbox_to_anchor=(1.20,0.90), loc = 2)
    #ax2.yaxis.grid(linestyle='-.', linewidth='0.5', color='red')

    plt.xticks()
    plt.title('QM and MM energies along the optimisation', y=1.05, loc='center')
    plt.savefig('QM_vs_MM_energies.png', transparent=False, dpi=300, bbox_inches='tight')
    plt.show()

    print("Plot for QM and MM energies saved")


elif type == 'scan':
    i=0; j=0; lines= []

    if qm_prog == 'turbomole':
        for line in file:
            j +=1
            if line.find('Cycle      1') != -1:
                i+= 1
                globals()['table_%s' % i] = [[],[],[],[]]
                lines.append(j)
            if line.find('Cycle ') != -1:
                globals()['table_%s' % i][0].append(int(line[9:12]))
            if line.find("QM/MM Energy:") != -1:
                globals()['table_%s' % i][1].append(float(line[13:-7]))
            if line.find('Contribution to energy from                  turbomole:') != -1:
                globals()['table_%s' % i][2].append(float(line[55:-7]))
            if line.find('Contribution to energy from                    dl_poly:') != -1:
                globals()['table_%s' % i][3].append(float(line[55:-7]))  

    if qm_prog == 'gaussian':
        for line in file:
            j +=1
            if line.find('Cycle      1') != -1:
                i+= 1
                globals()['table_%s' % i] = [[],[],[],[]]
                lines.append(j)
            if line.find('Cycle ') != -1:
                globals()['table_%s' % i][0].append(int(line[9:12]))
            if line.find("QM/MM Energy:") != -1:
                globals()['table_%s' % i][1].append(float(line[13:-7]))
            if line.find('Contribution to energy from                   gaussian:') != -1:
                globals()['table_%s' % i][2].append(float(line[55:-7]))
            if line.find('Contribution to energy from                    dl_poly:') != -1:
                globals()['table_%s' % i][3].append(float(line[55:-7]))  


    while True:
        try :
            globals()['table_%s' % i][1][len(globals()['table_%s' % i][0])-1]
            break
        except IndexError:
            del globals()['table_%s' % i][0][-1]
            break

    while True:
        try:
            if argv[2] == 'csv':
                csv = open('%s_energies.csv' % file_name, 'w+')
                csv.write("Point, Cycle, QM/MM energy, QM energy, MM energy\n")
                for k in range(1, i+1):
                    #globals()['table_%s' % i][0]
                    for n in range(0,len(globals()['table_%s' % k][0])):
                        csv.write(str(k)+','+str(globals()['table_%s' % k][0][n])+','+str(globals()['table_%s' % k][1][n])+','+str(globals()['table_%s' % k][2][n])+','+str(globals()['table_%s' % k][3][n])+'\n')
                csv.close()
                print("CSV file saved.")
                break
            else :
                break
        except IndexError :
            break

    table = [[], [], [], []]
    for k in range(1,i):
        table[0].extend(globals()['table_%s' % k][0])
        table[1].extend(globals()['table_%s' % k][1])
        table[2].extend(globals()['table_%s' % k][2])
        table[3].extend(globals()['table_%s' % k][3])




    len_0 = 0
    for k in range(1, i+1):
        globals()['len_%s' % k] = len(globals()['table_%s' % k][1]) + globals()['len_%s' % (k-1)]


    for k in range(1,i+1):
        plt.plot(range(globals()['len_%s' % (k-1)]+1,globals()['len_%s' % k]+1), globals()['table_%s' % k][1])
    plt.xlabel('Number of cycles')
    plt.ylabel('QM/MM Potential energy (Hartree)')
    plt.legend(["QM/MM\nenergy"], bbox_to_anchor=(1.02,1), loc = 2)

    plt.title('QM/MM energies along the scan (%s points)' % i, y=1.05, loc='center')
    plt.savefig('QMMM_energies.png', transparent=False, dpi=300, bbox_inches='tight')
    plt.show()

    print("Plot for QM/MM energies saved")


    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(table[1])), table[2][:], color='blue')#, marker='o')
    ax1.set_xlabel('Number of cycles')
    ax1.set_ylabel('QM Potential energy (Hartree)', color='blue')
    ax1.tick_params(axis='y', colors='blue')
    ax1.legend(["QM\nenergy"], bbox_to_anchor=(1.20,1), loc = 2)

    ax2=ax1.twinx()
    ax2.plot(range(0,len(table[1])), table[3][:], color='red')#, marker='^')
    ax2.set_ylabel('Potential energy (Hartree)', color='red')
    ax2.tick_params(axis='y', colors='red')
    ax2.legend(["MM\nenergy"], bbox_to_anchor=(1.20,0.90), loc = 2)

    plt.xticks()
    plt.title('QM and MM energies along the scan (%s points)' % i, y=1.05, loc='center')
    plt.savefig('QM_vs_MM_energies.png', transparent=False, dpi=300, bbox_inches='tight')
    plt.show()

    print("Plot for QM and MM energies saved")
