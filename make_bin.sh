#!/bin/bash

#simple script to pack up a precompiled binary package, with the boost thread
# library statically linked in.

echo "packing up $1.tar.gz, using boost in $2, linking against $3"
mkdir $1
make clean
./configure --enable-intel64 --enable-vectorize --with-boost=$2 --with-boost-thread=$3
make
cp src/cufflinks $1
cp src/cuffcompare $1
cp src/cuffdiff $1

cp README $1
cp COPYING $1
cp AUTHORS $1

tar cvfz $1.tar.gz $1