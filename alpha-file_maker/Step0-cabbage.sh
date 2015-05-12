#!/bin/bash
#Step0-cabbage.sh: This wrapper script contains multiple options for creating/handling cabbage-files:

#Standard error message (general):
if [[ $# -lt 2 || $# -gt 3 ]]
then
echo "usage: Step0-cabbage.sh <Input file/directory> <parameter-1> <parameter-2>"
echo "parameter-1 is non-optional. There are three choices (types):"
echo "  -type1: Converts individual <pocket.pdb files> into individual <cabbage unit files>"
echo "  -type2: Concatenates individual <cabbage unit files> (in IP-directory) into a <cabbage datafile>"
echo "  -type3: Concatenates individual <pocket.pdb files> (in IP-directory) into a <cabbage datafile>"
echo "  -decode: Decodes a given <cabbage unit file>, outputting in a human-readable format"
echo "parameter-2 is optional. There is only one choice:"
echo "  -runPM: Directly run PocketMatch (type1, serial) on outputted <cabbage datafile>"
exit
fi

#pre-execution cleanup:
rm outfile.cabbage

if [[ "$2" == "-type1" ]]
then
  anchor=`pwd`
  ./Step0-cabbage-core $1 > $1.cabbage
  mv $1.cabbage $anchor
fi

if [[ "$2" == "-type2" ]]
then
  anchor=`pwd`
  cd $1
  address=`pwd`
  cd $anchor
  #Concatenate <cabbage unit files> into <cabbage datafile>:
  for i in `ls $address`
  do
    cat $address/$i >> outfile.cabbage
  done
  #Insert END-OF-FILE string into <cabbage datafile>:
  ./Step0-END-FILE >> outfile.cabbage
fi

if [[ "$2" == "-type3" ]]
then
  anchor=`pwd`
  cd $1
  address=`pwd`
  cd $anchor
  #Convert all <pocket.pdb files> into individual <cabbage unit files>:
  #Concatenate <cabbage unit files> into <cabbage datafile>:
  for i in `ls $address`
  do
    ./Step0-cabbage-core $address/$i >> outfile.cabbage
  done
  #Insert END-OF-FILE string into <cabbage datafile>:
  ./Step0-END-FILE >> outfile.cabbage
fi

if [[ "$2" == "-decode" ]]
then
  anchor=`pwd`
  ./Step0-cabbage_decoder $1 > $1.txt
  mv $1.txt $anchor
fi

if [[ "$3" == "-runPM" ]]
then
  if [[ "$2" == "-type2" || "$2" == "-type3" ]]
  then
    ../Step3-PM_hybrid outfile.cabbage
  else
    echo "<Step0-cabbage.sh> error: '-runPM' compatible with '-type2' '-type3' only."
    exit
  fi
  exit
fi

  












