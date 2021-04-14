#! /users/mcanyelles/.conda/envs/py_env/bin/python

import os
import subprocess as sp
from sys import argv, exit

if len(argv) == 2:
    ID = argv[1]

else :
    print('Type the job\'s ID as the first argument')
    exit()

def main(ID):
    command = "scontrol show job " + str(ID) + ' | grep WorkDir'

    work_dir = str(sp.getoutput(command))[11:]

    print(work_dir)

main(ID)
