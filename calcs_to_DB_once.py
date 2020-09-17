# coding: utf-8

import paramiko
from password import password
from prevs import prev_running_IDs, prev_waiting_IDs
from yaml import safe_load, dump
from time import sleep
import os

host = "picard.uab.es"
port = 22022
username = "mcanyelles"

db_filename = 'calcs.db'


def set_ssh_connection(host=host, port=port, username=username, password=password):
    '''
    Function for setting the SSH connection with the server
    '''

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, port, username, password)


    return ssh

def close_ssh_connection(ssh): ssh.close()

def read_IDs(ssh, username=username):
    command = "qstat -u %s" % username

    stdin, stdout, stderr = ssh.exec_command(command)
    lines = stdout.readlines()
    if lines != []:
        running_IDs = []
        waiting_IDs = []
        for i in lines[5:]:
            if i[86:87] == 'Q':
                waiting_IDs.append(i[:6])
            elif i[86:87] == 'R':
                running_IDs.append(i[:6])

        return running_IDs, waiting_IDs

    else :
        print('Empty prompt!')
        return None, None

def read_info_running(ID, ssh):
    info = {'ID'            : '',
            'job_name'      : '',
            'server'        : '',
            'queue'         : '',
            'node'          : '',
            'status'        : '',
            'max_time'      : '',
            'spent_time'    : '',
            'starting_time' : '',
            'procs'         : '',
            'work_dir'      : ''}

    command = 'qstat -f %s.kirk.uab.es' % ID
    stdin, stdout, stderr = ssh.exec_command(command)
    lines = stdout.readlines()

    for line in range(len(lines)):
        info['ID'] = ID
        info['server'] = 'Picard'
        if 'Job_Name' in lines[line]:
            info['job_name'] = lines[line][15:-1]

        if 'job_state' in lines[line]:
            info['status'] = 'Running'

        if 'queue' in lines[line]:
            info['queue'] = lines[line][12:-1]

        if 'exec_host' in lines[line]:
            loc = lines[line].find('.')
            info['node'] = lines[line][16:loc]

        if 'Resource_List.walltime' in lines[line]:
            info['max_time'] = lines[line][29:-1]

        if 'resources_used.walltime' in lines[line]:
            info['spent_time'] = lines[line][30:-1]

        if 'start_time' in lines[line]:
            info['starting_time'] = lines[line][17:-1]

        if 'Resource_List.nodes' in lines[line]:
            loc = lines[line].find('ppn')
            info['procs'] = lines[line][loc+4:-1]

        if 'init_work_dir' in lines[line]:
            info['work_dir'] = lines[line][20:-1] + ((lines[line + 1].replace(' ', '')).replace('\t', '')).replace('\n', '')

    return info

def read_info_waiting(ID, ssh):
    info = {'ID'            : '',
        'job_name'      : '',
        'server'        : '',
        'queue'         : '',
        'node'          : '',
        'status'        : '',
        'max_time'      : '',
        'spent_time'    : '',
        'starting_time' : '',
        'procs'         : '',
        'work_dir'      : ''}

    command = 'qstat -f %s.kirk.uab.es' % ID
    stdin, stdout, stderr = ssh.exec_command(command)
    lines = stdout.readlines()

    for line in range(len(lines)):
        info['ID'] = ID
        info['server'] = 'Picard'
        if 'Job_Name' in lines[line]:
            info['job_name'] = lines[line][15:-1]

        if 'job_state' in lines[line]:
            info['status'] = 'Waiting'

        if 'queue' in lines[line]:
            info['queue'] = lines[line][12:-1]

        if 'Resource_List.walltime' in lines[line]:
            info['max_time'] = lines[line][29:-1]

        if 'Resource_List.nodes' in lines[line]:
            loc = lines[line].find('ppn')
            info['procs'] = lines[line][loc+4:-1]

        if 'init_work_dir' in lines[line]:
            info['work_dir'] = lines[line][20:-1] + ((lines[line + 1].replace(' ', '')).replace('\t', '')).replace('\n', '')


    return info


def dump_to_db(db_dict, db_filename=db_filename):

    db_yaml = dump(db_dict)
    db_file = open(db_filename + '.tmp', 'w')
    db_file.write(db_yaml)
    db_file.close()

    os.remove(db_filename)
    os.rename(db_filename + '.tmp', db_filename)

def save_prevs(prev_running_IDs, prev_waiting_IDs):

    f = open('prevs.py.tmp', 'w')
    f.write('prev_running_IDs = %s\nprev_waiting_IDs = %s' % (str(prev_running_IDs), str(prev_running_IDs)))
    f.close()

    os.remove('prevs.py')
    os.rename('prevs.py.tmp', 'prevs.py')


def main(prev_running_IDs, prev_waiting_IDs):
    ssh = set_ssh_connection()

    running_IDs, waiting_IDs = read_IDs(ssh)
    print('running: ', running_IDs)
    print('waiting: ', waiting_IDs)

    if running_IDs != None and waiting_IDs != None:
        db = safe_load(open(db_filename, 'r'))

        for ID in waiting_IDs:
            if prev_running_IDs != None:
                if ID not in prev_waiting_IDs:
                    db['P.' + ID] = read_info_waiting(ID, ssh)
            else :
                db['P.' + ID] = read_info_waiting(ID, ssh)


        for ID in running_IDs:
            db['P.' + ID] = read_info_running(ID, ssh)
        if prev_running_IDs != None:
            for ID in prev_running_IDs:
                if ID not in running_IDs:
                    db['P.' + ID]['status'] = 'Completed'


        #dump_to_YAML
        if running_IDs == []:
            prev_running_IDs = None
        else :
            prev_running_IDs = running_IDs

        if waiting_IDs == []:
            prev_running_IDs = None
        else :
            prev_waiting_IDs = waiting_IDs


        dump_to_db(db)
        print('DB updated!')
        save_prevs(prev_running_IDs, prev_waiting_IDs)


        close_ssh_connection(ssh)

if __name__ == '__main__':
    main(prev_running_IDs, prev_waiting_IDs)
