# coding: utf-8

import paramiko
from password import password
from yaml import safe_load, dump
from time import sleep
import os


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
    
    for line in lines:
        info['ID'] = ID
        info['server'] = 'Picard'
        if 'Job_Name' in line:
            info['job_name'] = line[15:-1]

        if 'job_state' in line:
            info['status'] = 'Running'

        if 'queue' in line:
            info['queue'] = line[12:-1]

        if 'exec_host' in line:
            loc = line.find('.')
            info['node'] = line[16:loc]

        if 'Resource_List.walltime' in line:
            info['max_time'] = line[29:-1]

        if 'resources_used.walltime' in line:
            info['spent_time'] = line[30:-1]

        if 'start_time' in line:
            info['starting_time'] = line[17:-1]

        if 'Resource_List.nodes' in line:
            loc = line.find('ppn')
            info['procs'] = line[loc+4:-1]

        if 'init_work_dir' in line:
            info['work_dir'] = line[20:-1]

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

    
    for line in lines:
        info['ID'] = ID
        info['server'] = 'Picard'
        if 'Job_Name' in line:
            info['job_name'] = line[15:-1]

        if 'job_state' in line:
            info['status'] = 'Waiting'

        if 'queue' in line:
            info['queue'] = line[12:-1]

        if 'Resource_List.walltime' in line:
            info['max_time'] = line[29:-1]

        if 'Resource_List.nodes' in line:
            loc = line.find('ppn')
            info['procs'] = line[loc+4:-1]

        if 'init_work_dir' in line:
            info['work_dir'] = line[20:-1]

    return info


def dump_to_db(db_dict, db_filename=db_filename):
    
    db_yaml = dump(db_dict)
    db_file = open(db_filename + '.tmp', 'w')
    db_file.write(db_yaml)
    db_file.close()
    
    os.remove(db_filename)
    os.rename(db_filename + '.tmp', db_filename)
    

try : 
    del prev_running_IDs
    del prev_waiting_IDs
except NameError:
    pass

while True:
    ssh = set_ssh_connection()

    running_IDs, waiting_IDs = read_IDs(ssh)
    print('running: ', running_IDs)
    print('waiting: ', waiting_IDs)
    
    if running_IDs != None and waiting_IDs != None:
        db = safe_load(open(db_filename, 'r'))        
        
        for ID in waiting_IDs:
            try: 
                if ID not in prev_waiting_IDs:
                    db[ID] = read_info_waiting(ID, ssh) 
                else :
                    print('avoided!')
            except NameError:
                db[ID] = read_info_waiting(ID, ssh)


        try :
            for ID in running_IDs:
                db[ID] = read_info_running(ID, ssh)
            
            for ID in prev_running_IDs:
                if ID not in running_IDs:
                    db[ID]['status'] = 'Completed'
                    
        except NameError:
            for ID in running_IDs:
                db[ID] = read_info_running(ID, ssh)


        #dump_to_YAML

        prev_running_IDs = running_IDs
        prev_waiting_IDs = waiting_IDs

        
        dump_to_db(db)
        print('DB updated!')

        
        close_ssh_connection(ssh)
   
        sleep(900) # 15 min = 900 s

    else :
        del prev_running_IDs
        del prev_waiting_IDs
        sleep(3600) # 1 h
