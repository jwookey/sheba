#!/bin/bash
#
# Script to bundle a distribution copy of SHEBA to /tmp
#

rm -rf /tmp/sheba
cp -R ../sheba /tmp
cd /tmp/sheba
make clean
rm -rf .git runtests run_unittests tests
cd /tmp
tar zcvf sheba.tgz sheba

