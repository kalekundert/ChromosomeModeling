#!/usr/bin/env sh

mkdir $1
cd $1
ln -s ../libraries .
ln -s ../chromosome_modeler .
ln -s ../visualize_model.sh .
