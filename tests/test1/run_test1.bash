#!/bin/bash
#
#  Run a test of the sheba core executable, to make sure it still returns the
#  correct result from a synthetic waveform. 
#
#  Requires: SWAV_test1.BHE, SWAV_test1.BHN, SWAV_test1.BHZ
#  This uses the /tmp directory. It also assumes that sheba_exec is availble
#  on the path. 
#
#  Files were created with: sacsplitwave -op 30 1.3 -noise 0.1
#  then filtered in SAC with: bp bu co 0.01 0.3 n 2 p 2

# get rid of any old test directory
rm -rf /tmp/sheba_autotest1

# create a new one
mkdir /tmp/sheba_autotest1

# copy the data there. 
cp SWAV_test1.BH? /tmp/sheba_autotest1

# save the current directory.
here=`pwd`

# go to the test directory
cd /tmp/sheba_autotest1

# build the input file
cat << EOF > sheba.in
# sheba.in
SWAV_test1
BHE
BHN
BHZ
1
1 1
4.0
0
0
EOF

echo "##################" > sheba_autotest1.log
echo "## SHEBA OUTPUT ##" >> sheba_autotest1.log
echo "##################" >> sheba_autotest1.log
echo "" >> sheba_autotest1.log

# run the SHEBA executable
sheba_exec >> sheba_autotest1.log

echo "##################" > sheba_autotest1.log
echo "## TEST RESULTS ##" >> sheba_autotest1.log
echo "##################" >> sheba_autotest1.log
echo "" >> sheba_autotest1.log


### EXAMINE THE MAIN RESULT.
# Benchmark result is: FAST = 34 +/- 3.25 ; TLAG = 1.20 +/- 0.06

echo "Benchmark result is: FAST = 34 +/- 3.25 ; TLAG = 1.20 +/- 0.06" >> sheba_autotest1.log
tail -1 SWAV_test1_sheba.result | awk '
{\
printf("Checking SHEBA outputs:"); \
printf("   Checking "); if (abs($11-34)>1) print "FAST: FAILED."; else print "FAST: PASSED"; \
printf("   Checking "); if (abs($12-3.50)>1) print "DFAST: FAILED."; else print "DFAST: PASSED"; \
printf("   Checking "); if (abs($13-1.2)>0.05) print "TLAG: FAILED."; else print "TLAG: PASSED"; \
printf("   Checking "); if (abs($14-0.08)>0.01) print "DTLAG: FAILED."; else print "DTLAG: PASSED"; \
}
function abs(value)
{
  return (value<0?-value:value) ;
}
' >> sheba_autotest1.log

# finally, return a value indicating whether tests were passed or failed. 
NFAILS=`grep "FAILED" sheba_autotest1.log | wc -l`
if [ $NFAILS -gt 0 ] 
then
   echo "test1 FAILED"
else
   echo "test1 PASSED"
fi   


# done
