# -*- coding: utf-8 -*-
import argparse
import os
import sys
from yaml import load, dump, Loader, Dumper
from platform import node

mail = 'miquel.canyelles@uab.cat'

class MissingInformationInCLIError(Exception):

    #"Some information is missing. Job name (-j, --jobname), file name (-f, --filename), and source dir (-s, --sourcedir) are required.\nThey can also be input using a YAML file containing the long keys and their values and input as the -i flag."

    def __init__(self, parser, message):
        self.message = message
        super().__init__(self.message)
        parser.print_help()

    pass

class MissingInformationInYAMLError(Exception):
    #"Some information is missing in the YAML file. jobname, filename and/or sourcedir are missing in the YAML file. \nThey can also be input as flags (-j, -f and -sd, respectively)."
    #parser.print_help()
    def __init__(self, parser, message):
        self.message = message
        super().__init__(self.message)
        parser.print_help()

    pass


def arguments_parser():
    """
    DESCRIPTION:
        Function for parsing input arguments and storing them as a dictionary.
    """

    parser = argparse.ArgumentParser(prog='chemsh_jobsender',
                description='Arguments for using chemsh_jobsender as CLI.')

    parser.add_argument(
        '-i', '--input_file',
        type=str, help='Input a YAML file containing the information for faster sending'
    )

    parser.add_argument(
        '-j', '--jobname',
        type=str, action='store', help='Input the name of the job',
    )

    parser.add_argument(
        '-f', '--filename',
        type=str, action='store', help='Specify the filename of the ChemShell job',
    )

    parser.add_argument(
        '-m', '--mail',
        type=str, action='store', default=mail, help='Mail of the user to send SLURM notifications',
    )

    parser.add_argument(
        '-n', '--nodes',
        type=int, action='store', default=1, help='Select number of nodes',
    )

    parser.add_argument(
        '-p', '--processors',
        type=int, action='store', default=8, help='Select the number of processors for the calculation',
    )

    parser.add_argument(
        '-q', '--queue',
        type=str, action='store', help='Select the queue where the job has to be send. For mirfak, LocalQ is the default; for CSUC, std is the default.',
    )

    parser.add_argument(
        '-t', '--timelimit',
        type=int, action='store', default='30', help='Time limit in days',
    )

    parser.add_argument(
        '-d', '--jobdir',
        type=str, action='store', default="$SLURM_SUBMIT_DIR", help='Specify the home path of the job. By default, it is the sending dir specified through the SLURM_SUMBIT_DIR var',
    )

    parser.add_argument(
        '-s', '--sourcedir',
        nargs='+', action='store', help='Specify the path to the source files for the ChemShell execution',
    )

    parser.add_argument(
        '-y', '--save_yaml',
        action='store_true', default=False, help="Activate the saving of the configuration as a YAML file. The name will be 'jobconfig.yaml'."
    )

    parser.add_argument(
        '-ns', '--not_send_job',
        action='store_true', default=False, help="Deactivates sending the job."
    )

    args = vars(parser.parse_args())

    if args['input_file'] == None:
        if args['jobname'] == None:
            raise MissingInformationInCLIError(parser, 'jobname (-j) is missing')
        elif args['sourcedir'] == None:
            raise MissingInformationInCLIError(parser, 'sourcedir (-s) is missing')

        elif args['filename'] == None:
            raise MissingInformationInCLIError(parser, 'filename (-f) is missing')

    elif args['input_file'] != None:
        f = open(args['input_file'])
        yamlargs = load(f, Loader=Loader)
        f.close()

        for arg in list(yamlargs.keys()):
            if yamlargs[arg] == None:
                raise MissingInformationInYAMLError(parser, "Check the input YAML, some required information is missing (ojbname, sourcedir and/or filename.")

        for arg in list(yamlargs.keys()):
            args[arg] = yamlargs[arg]

    if args['queue'] == None:
        if node() == 'mirfak':
            args['queue'] = 'LocalQ'
        elif node() in ('login1', 'login2'):
            args['queue'] = 'std'
        else :
            args['queue'] = input('Type the name of the queue where the job has to be send: ')

    if '.' in args['filename']:
        args['filename'] = str(args['filename'].split('.')[0])

    if args['not_send_job'] == True:
        send = False
    elif args['not_send_job'] == False:
        send = True

    if args['save_yaml'] == True:
        del args['save_yaml']
        del args['input_file']
        del args['not_send_job']

        with open('job_config.yaml', 'w') as f:
            f.write(dump(args, Dumper=Dumper))

    args['module']  = ''
    args['scratch'] = ''

    return args, send

def file_builder(args):
    """
    DESCRIPTION:
        Function for building the SLURM script

    NEEDED ARGS:
        - filename
        - jobdir
        - jobname
        - mail
        - nodes
        - processors
        - queue
        - sourcedir
        - timelimit
    """

    if node() == 'mirfak':
        args['module'] = 'module load chemshell'
        args['scratch'] = "$SCRATCH={}/SCRATCH\nmkdir $SCRATCH".format(args['jobdir'])
    elif node() in ('login1', 'login2'): # CSUC
        args['module'] = '. /home/uabqut17/soft/chemshell_turbomole.sh'

    else:
        args['module'] = '# No module or script loaded'

    sourcedirs = ''
    for sd in args['sourcedir']:
        if sd[-1] == '/':
            sd = sd[:-1]

        sourcedirs += "cp {}/* $SCRATCH\n".format(sd)

    args['sourcedir'] = sourcedirs[:-1]

    script = """#! /bin/bash
#SBATCH -J {jobname}
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.log
#SBATCH --mail-user={mail}
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH -p {queue}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={processors}
#SBATCH -t {timelimit}-00:00

{module}

# Environment variables
export OMP_NUM_THREADS={processors}
export PARNODES={processors}
export PARA_ARCH=SMP

# Copy of files to SCRATCH and create subfolder
{scratch}
cp -f {jobdir}/* $SCRATCH
{sourcedir}
mkdir {jobdir}/structures

# Execute ChemShell
chemsh {filename}.chm > {jobdir}/{filename}.out

# Copy files
cp $SCRATCH/*.pdb {jobdir}/structures
cp $SCRATCH/*.c   {jobdir}/structures
cp $SCRATCH/*.plt {jobdir}/structures
    """.format(**args)

    f = open('script.sh', 'w')
    f.write(script)
    f.close()


if __name__ == '__main__':
    args, send = arguments_parser()

    file_builder(args)

    if send == True:
        os.system('sbatch script.sh')
        print("Job titled '{}' has been sent.".format(args['jobname']))
    else :
        print("Job has not been sent, but SLURM script has been sent as 'script.sh'.")


