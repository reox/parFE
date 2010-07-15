#!/bin/bash

if [ $# -ne "1" ]
then
  echo "Usage: $0 filename"
  exit 1
fi

if [  ! -f "$1" ]
then
  echo "File not found"
  exit 2
fi

./convert_to_newformat.py $1
sed -e 's/Node number/Node_number/' -e 's/Nodal freedom/Nodal_freedom/' -i $1 
