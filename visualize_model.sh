#!/usr/bin/env sh

if [[ $# -eq 0 ]]; then
    pdb='toy_model.pdb'
else
    pdb=$1
fi

model=$(basename "$pdb" .pdb)

fork pymol -qx $pdb                                                       \
    -d 'set seq_view, off'                                                  \
    -d 'set connect_mode, 1'                                                \
    -d 'set sphere_scale, 0.2'                                              \
    -d 'split_chains all'                                                   \
    -d "set_name ${model}_A, reference"                                     \
    -d "set_name ${model}_B, restraints"                                    \
    -d "set_name ${model}_C, optimized"                                     \
    -d "delete ${model}"                                                    \
    -d 'as spheres, restraints'                                             \
    -d 'color cyan, optimized'                                              \
