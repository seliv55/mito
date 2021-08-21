#!/bin/sh
#edata="../filesimid/files/mark2" #13C labeling 
edata="a011110rbm.txt" #experimental data
init="i0" #initial values
par="1" #set of parameters
mode="0" #-Fit with Simulated annealing, find -Number of independent parameters,
tst=yes

while getopts ":e:i:p:m:" opt; do
  case $opt in
    e) edata=$OPTARG;;
    i) init=$OPTARG;;
    p) par=$OPTARG;;
    m) mode=$OPTARG;;
    *)
      echo "Valid options: -e edata, -i init value, -p parameters, -m mode" 
      cat help
      tst=no
      ;;
  esac
done
if [ $tst = yes ]
then
./mitodyn.out $edata $init $par $mode
fi
