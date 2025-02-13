DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#TODO: should not need this see how to add the em to python path
${DIR}/reinstall.sh
nosetests