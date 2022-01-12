#!/bin/tcsh -f

if ( $#argv < 1 ) then 
	echo "Usage : $0 src_data_file src_dir";
	exit;
endif


set src_data_file = $1;
set src_dir = $2;

set name = `basename $0 .sh`; 
set idd = `basename $src_data_file`; 

echo  "start $idd $PWD $0 $* -> $name" ;

set q = '"';


cp -fv $src_data_file ./input.root


echo '########################################'
echo "Fun4All_ReadDST_eIDML.C"
echo '########################################'

cat Fun4All_ReadDST_eIDML.C


echo '########################################'
echo "Run"
echo '########################################'

echo ./input.root  > ./input.lst

set q = '"';

root -l -b -q "Fun4All_ReadDST_eIDML.C(10000000, ${q}./input.lst${q})" ;


echo '########################################'
echo "Copy out"
echo '########################################'

ls -lhv

cp -fv *_BECAL_*.root  ${src_data_file}_BECAL.root;
cp -fv *_EEMC_*.root  ${src_data_file}_EEMC.root;
cp -fv *_FEMC_*.root  ${src_data_file}_FEMC.root;

date ;
echo "end" ;
