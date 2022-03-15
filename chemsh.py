#! /usr/bin/env python3.5

# DESCRIPTION: Script for sending ChemShell jobs to Kirk cluster
# AUTHOR:      Miquel Canyelles (github.com/mikicanyelles)

from os import getcwd, mkdir, path, popen
from sys import exit
from argparse import ArgumentParser
import subprocess

# environment vars
filename       = 'script_chemsh.pbs'
walltime_borg2 = 1800000
walltime_borg3 = None
walltime_borg4 = 21600000
#user           = getlogin()
#hostname       = ''

# args
parser = ArgumentParser(
    description= 'Script for sending ChemShell jobs to Kirk cluster'
)

parser.add_argument('-n', '--nproc', type=int, required=True, help='Number of processors per node. Multiples of 8 work in borg3 (up to 64) and borg4 (up to 32), multiples of 12 work in borg2 (up to 48), 1, 2 or 4 procs can be used in all three queues.')
parser.add_argument('-J', '--jobname', help='Name of the submited job. If none give, it will take the folder\'s name.')
parser.add_argument('-N', '--nosub', action='store_true', help='Do not submit, only create the pbs script.')
parser.add_argument('-q', '--queue', choices=['borg2', 'borg3', 'borg4'], required=True, help='Queue to submit to.')
parser.add_argument('-s', '--source_dir', help='Path to the source dir if files are not in folder.')
parser.add_argument('-sf', '--source_dir_frame', help='Arg to specify a subfolder inside source_dir containing file for the job (useful in case of having source file organised by frames into folders.')
parser.add_argument('-S', '--save_scratch', action='store_true', help='Save all file in scratch in local (in a subfolder).')
parser.add_argument('-v', '--version', choices=['3.7', '3.4'], default='3.7', type=str)
parser.add_argument('-w', '--walltime', help='Custom walltime (in s). Max walltimes are: 1800000 for borg2 and 3600000 for borg4, borg3 has no walltime limit.')
parser.add_argument('input_file', help='Input file')
parser.add_argument('output_file', help='Output file')


group = parser.add_mutually_exclusive_group()
group.add_argument('-G', '--gaussian', action='store_true', help='Use Gaussian09 as the QM backend.')
group.add_argument('-T', '--turbomole', action='store_true', help='Use Turbomole7 as the QM backend.')


args = parser.parse_args()

# Queue selection
def configure_head():

    head = {
        'jobname'  : None,
        'nproc'    : None,
        'queue'    : None,
        'walltime' : None
        }

    head['queue'] = args.queue

    if args.queue == 'borg2':
        if args.nproc not in (1, 2, 4, 12, 24, 36, 48):
            print('The number of processors selected cannot be used in borg2.')
            exit(1)

        else :
            head['nproc'] = args.nproc

        if args.walltime != None:
            if args.walltime <= walltime_borg2:
                head['walltime'] = '#PBS -l walltime=' + str(args.walltime)
            else :
                print('Max walltime has been used instead of the given walltime')
                head['walltime'] = '#PBS -l walltime=1800000'
        else :
            head['walltime'] = '#PBS -l walltime=1800000'

    if args.queue == 'borg3':
        if args.nproc not in (1, 2, 4, 16, 32, 64):
            print('The number of processors selected cannot be used in borg3.')
            exit(1)
        else :
            head['nproc'] = args.nproc

        head['walltime'] = ''


    if args.queue == 'borg4':
        if args.nproc not in (1, 2, 4, 16, 32):
            print('The number of processors selected cannot be used in borg4.')
            exit(1)
        else :
            head['nproc'] = args.nproc

        if args.walltime != None:
            if args.walltime <= walltime_borg4:
                head['walltime'] = '#PBS -l walltime=' + str(args.walltime)
            else :
                print('Max walltime has been used instead of the given walltime')
                head['walltime'] = '#PBS -l walltime=36000000'

        else :
            head['walltime'] = '#PBS -l walltime=%s' % str(walltime_borg4)


    if args.jobname != None:
        head['jobname'] = args.jobname

    else :
        head['jobname'] = args.input_file

    return head


def configure_modules():
    modules = {
        'chemsh'    : None,
        'mpirun'    : None,
        'para_arch' : None,
        'parnodes'  : None,
        'qm'        : None
        }

    # version
    if args.version == '3.7':
        modules['chemsh'] = 'chemshell/3.7_gcc6.3.0_ompi-2.0.1'

    elif args.version == '3.4':
        modules['chemsh'] = 'chemshell-3.4.2-ifort'

    # QM backend
    if args.turbomole == True and args.gaussian == False:
        modules['qm']       = 'turbomole7.0'
        modules['parnodes'] = 'export PARNODES=%s' % str(args.nproc)

        # parallel motor
        if args.queue == 'borg2':
            modules['para_arch'] = 'export PARA_ARCH=MPI'

        else :
            modules['para_arch'] = 'export PARA_ARCH=SMP'

    elif args.turbomole == False and args.gaussian == True:
        modules['qm']        = 'g09_D.01_pgi11.9-ISTANBUL'
        modules['parnodes']  = ''
        modules['para_arch'] = ''

    if args.queue == 'borg4':
        modules['mpirun'] = 'module load openmpi/2.0.1-tm_gcc-6.3.0-borg4'

    else :
        modules['mpirun'] = ''

    return modules


def configure_body():

    body = {
        'input_file'       : None,
        'output_file'      : None,
        'source_dir'       : None,
        'source_dir_frame' : None,
        }

    # look for source_dir and frame and check if they exist
    if args.source_dir != None:
        if path.exists(args.source_dir):
            body['source_dir'] = 'cp -f %s* $TMPDIR' % args.source_dir

            if args.source_dir_frame != None:
                if path.exists(args.source_dir + args.source_dir_frame):
                    body['source_dir_frame'] = 'cp -f %s%s/* $TMPDIR' % (args.source_dir, args.source_dir_frame)

                else :
                        print('Frame in source dir does not exist. Check the frame number.')
                        exit(1)
            else :
                body['source_dir_frame'] = None

        else :
            print('Source dir does not existis. Check the route to source files.')
            exit(1)

    else :
        body['source_dir'] = None


    body['input'] = args.input_file

    if args.output_file == None:
        body['output_file'] = '.'.join(body['input_file'].split('.')[:-1]) + '.out'

    else :
        body['output_file'] = args.output_file


    return body

def configure_tail():

    tail = {
        'save_scratch' : None
        }

    if args.save_scratch == True:
        tail['save_scratch'] = True
        if not path.exists('./node'):
            mkdir('node')

    else :
        tail['save_scratch'] = ''

    if not path.exists('./structures'):
        mkdir('structures')

    return tail


def create_pbs(head, modules, body, tail):


    template = """#PBS -q {queue}
#PBS -N {jobname}
#PBS -l nodes=1:ppn={nproc}
#PBS -k oe
#PBS -r n
{walltime}


. ./QFcomm/modules.profile
module load {chemshell}
{mpirun}
{paraarch}
module load {qm}
{parnodes}

export TMPDIR=/scratch/${{PBS_JOBNAME}}.${{PBS_JOBID::6}}.$USER

mkdir $TMPDIR

cp -f $PBS_O_WORKDIR/* $TMPDIR

{sourcedir}
{sourcedirframe}

INPUT={inputfile}
OUTPUT={outputfile}

cd $TMPDIR

chemsh $INPUT > $PBS_O_WORKDIR/$OUTPUT

cp -f $TMPDIR/*.pdb $PBS_O_WORKDIR/structures
cp -f $TMPDIR/*.c $PBS_O_WORKDIR/structures
cp -f $TMPDIR/PES.plt $PBS_O_WORKDIR
{savescratch}
"""

    filler = {
        'queue'    : head['queue'],
        'jobname'  : head['jobname'],
        'nproc'    : head['nproc'],
        'walltime' : head['walltime'],

        'chemshell' : modules['chemsh'],
        'mpirun'    : modules['mpirun'],
        'paraarch'  : modules['para_arch'],
        'qm'        : modules['qm'],
        'parnodes'  : modules['parnodes'],

        'sourcedir'      : body['source_dir'],
        'sourcedirframe' : body['source_dir_frame'],
        'inputfile'      : args.input_file, #body['input_file'],
        'outputfile'     : body['output_file'],

        'savescratch' : tail['save_scratch'],
    }

    with open('./' + filename, 'w') as f:
        f.write(template.format(**filler))
        f.close()


def submit_job(head, modules):
    if args.nosub == False:
#        jobid = popen("/usr/local/torque/bin/qsub script_chemshell.pbs").read()[:6]
        jobid = subprocess.Popen(["/usr/local/torque/bin/qsub", filename], stdout=subprocess.PIPE).communicate()
        print(jobid[0].decode('utf8'))
        jobid = str(jobid[0].decode('utf8'))[:6]

        print('\n')
        print('---------Job information---------')
        print('Job name:', head['jobname'])
        print('Job ID:', jobid)
        print('Program: Chemshell %s interfacing %s' % (args.version, modules['qm']))
        print('Job sent to:', head['queue'])

    else :
        print(filename, 'created.\n')



def main():
    head, modules, body, tail = configure_head(), configure_modules(), configure_body(), configure_tail()
    create_pbs(head, modules, body, tail)
    submit_job(head, modules)



main()


