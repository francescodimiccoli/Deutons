#!/bin/sh 

rm -rf file_check mc_summary livetime_summary

exe=" file_check mc_summary livetime_summary "

for i in $exe
do
  echo "Compiling $i ..."
  `root-config --cxx` -O3 -fPIC -std=c++11 `root-config --cflags` -I. -I../include GM_SubLibrary.C ${i}.cxx `root-config --libs` -lTreePlayer -lTMVA ../install/lib/libntpa.a -o ${i}
done

