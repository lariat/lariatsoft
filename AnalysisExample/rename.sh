#!/bin/sh

# Script to change the AnalysisExample name to something else.
# Usage: ./rename.sh <new-name>

newname=$1

# Rename the files.
mv AnalysisExample_module.cc  ${newname}_module.cc
mv AnalysisExample.fcl        ${newname}.fcl

# Search through all the files, replacing the text "AnalysisExample"
# with the new name. 
sed -i -e "s/AnalysisExample/${newname}/g" ${newname}.fcl ${newname}_module.cc CMakeLists.txt
