#!/bin/bash

python_version="python3.11"
python_version_path=$(which python3.11)
if [ -f $python_version_path ]; then
   echo "Python version $python_version found at $python_version_path"
   echo "You can continue with the INSTALLATION PROTOCOLE."
else
   echo "ERROR: Required python version is not in your system"	
   echo "Python version $python_version not found in your PATH"
   echo "You must install python version $python_version to run to finish setup."
   exit
fi

if [ -d "$PWD/../venv" ]; then
   echo "Virtual environment already exist."
   exit
else
   echo "Creating Virtual Environment"
   cd ..
   virtualenv --python=$python_version_path venv
   cd bin
fi
