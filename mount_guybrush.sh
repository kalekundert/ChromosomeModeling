#!/usr/bin/env sh

if [ $# -ne 1 ]; then
    echo "Usage: mount_guybrush <directory>"
    exit 1
fi

sshfs guybrush:chromo/$1 $1
