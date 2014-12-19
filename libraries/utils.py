#!/usr/bin/env python2

import os, shutil, json

params_path = 'inputs/params.json'

def clear_directory(directory):
    try: shutil.rmtree(directory)
    except OSError: pass
    os.mkdir(directory)

def clear_directories(*directories):
    for directory in directories:
        clear_directory(directory)

def running_on_cluster():
    from socket import gethostname
    return gethostname() == 'iq218'

def submit_job(command, params):
    from subprocess import Popen, STDIN, PIPE

    clear_directories('inputs', 'stdout', 'stderr')

    with open(params_path, 'w') as file:
        json.dump(params, file)

    qsub_command = (
            'qsub',
            '-cwd',
            '-o', 'stdout',
            '-e', 'stderr',
            '-l', 'h_rt=6:00:00',
            '-l', 'mem_free=1G',
            '-l', 'arch=linux-x64',
            '-l', 'netapp=1G',
            '-t', '1-{0}'.format(len(params)),
            '-N', command,
    )

    process = Popen(qsub_command, stdin=PIPE)
    process.stdin.write('#!/usr/bin/env sh')
    process.stdin.write('module load imp-fast')
    process.stdin.write('PYTHONPATH=.:$PYTHONPATH')
    process.stdin.write('/netapp/home/kale/.local/bin/python2.7 ' + command)
    process.stdin.close()

def read_params():
    task_id = int(os.environ['SGE_TASK_ID']) - 1
    with open(params_path) as file:
        return json.load(file)[task_id]


