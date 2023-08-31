#!/bin/bash

# Script to launch on bootstrapping on different R instances.

Rout="Rout"
Rfiles="Jags"

echo

# Standard syntax.
for a in 1 2 3 4 
do
  rm -f $Rfiles/toRun.r
  echo "index <- $a" >> $Rfiles/toRun.r
  cat $Rfiles/bootstrapping.R >> $Rfiles/toRun.r
  R CMD BATCH $Rfiles/toRun.r $Rout/outBoot$a.txt
  echo -n "$a" 
done  

rm -f $Rfiles/toRun.r

echo "Bootstrapping done"
