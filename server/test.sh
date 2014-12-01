#!/bin/bash

if [ ! -d "./p-env/" ]; then
    virtualenv p-env
fi

source ./p-env/bin/activate

TMP_DEPS=/tmp/temp_deps_${RANDOM}
pip freeze -l > ${TMP_DEPS}
if ! cmp ./requirements.txt ${TMP_DEPS} > /dev/null 2>&1
then
  echo "Installing Python dependencies ..."
  cat ${TMP_DEPS}
  pip install -r ./requirements.txt
fi

mkdir -p /tmp/uploads
cp src/sim/example.xlsx /tmp/uploads
NOSE_NOCAPTURE=1 OPTIMA_TEST_CFG="${PWD}/test.cfg" nosetests $@
