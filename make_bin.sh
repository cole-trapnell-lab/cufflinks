#!/bin/bash
#simple script to pack up a precompiled binary package, with the boost thread
# library statically linked in; on x86_64 Linux, libstc++ and libgcc are linked statically also
if [[ -z "$2" ]]; then
 echo -e "Usage:\n./make_bin.sh <package_base_name> <boost_prefix> [<bam_prefix> [<Eigen_prefix>]]"
 exit 1
fi

echo "packing up $1.tar.gz, using boost in $2, linking against $3 and using BAM in $4, using Eigen in $5"
mkdir $1
#make clean
make distclean
if [[ $(uname -m) = "x86_64" ]]; then
echo "Linking statically on x86_64 (only for gcc 4.5+).."
export LDFLAGS="-static-libgcc -static-libstdc++"
fi 
l2="$2"
l3="$3"
if [[ -z "$l3" ]]; then
  l3="$l2"
fi
l4="$4"
if [[ -z "$l4" ]]; then
  l4="$l2"
fi
# l5="$5"
# if [[ -z "$l5" ]]; then
#   l5="$l2"
# fi


#./configure --enable-intel64 --with-boost=$l2 --with-boost-thread=$l3 --with-bam=$l4 --with-eigen=$l5
./configure --with-boost=$l2 --with-boost-thread=$l2/lib/libboost_thread.a --with-bam=$l3 --with-eigen=$l4
make
cp src/cufflinks $1
cp src/cuffcompare $1
cp src/cuffdiff $1
cp src/cuffmerge $1/cuffmerge
cp src/gffread $1
cp src/gtf_to_sam $1
cp README $1
cp LICENSE $1
cp AUTHORS $1

tar cvfz $1.tar.gz $1
