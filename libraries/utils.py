#!/usr/bin/env python2

import os, shutil, subprocess

def clear_directory(directory):
    try: shutil.rmtree(directory)
    except OSError: pass
    os.mkdir(directory)

def running_on_chef():
    from socket import gethostname
    return gethostname() == 'chef.compbio.ucsf.edu'

def load_imp_modules(fast=True):
    command = 'module load imp{}'.format('-fast' if fast else '')
    subprocess.call(command, shell=True)

