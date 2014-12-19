#!/usr/bin/env python2

import os, shutil, subprocess, json

def clear_directory(directory):
    try: shutil.rmtree(directory)
    except OSError: pass
    os.mkdir(directory)

def running_on_cluster():
    from socket import gethostname
    return gethostname() == 'iq218'

def submit_job(command, params):
    with open('params.json', 'w') as file:
        json.dump(params, file)
    
    prefix = (
            'qsub',
            '-cwd',
            '-o', 'stdout',
            '-e', 'stderr',
            '-l', 'h_rt=6:00:00',
            '-l', 'mem_free=1G',
            '-l', 'arch=linux-x64',
            '-l', 'netapp=1G',
            '-t', '1-{0}'.format(len(params)),
    )
    if command[0].endswith('.py'):
        prefix += '-S', '/netapp/home/kale/.local/bin/python2.7'

    subprocess.call(map(str, prefix + command))

def read_params():
    task_id = int(os.environ['SGE_TASK_ID']) - 1
    with open('params.json') as file:
        return json.load(file)[task_id]


