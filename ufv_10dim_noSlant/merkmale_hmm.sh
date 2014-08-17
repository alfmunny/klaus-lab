#! /bin/bash

#
# NOTE: tools are found in $(ESMERALDA)/bin
#

SRC=$HOME/esmeralda/lab/klaus/lin 	#bearbeitete Bilder
DES=$HOME/esmeralda/lab/klaus/ufv
PGM=".pgm" #bearbeitete Bilder
#PGM=".pgm"
UFV=".ufv"

#
# frame extraction parameters
#
method=2 	# feature extraction method 
width=8	# frame width
x=4		# frame overlap :-(

cats=`(cd $SRC; ls -d [B]*)`
#cats=`ls $SRC/*.pgm` #bearbeitete Bilder
#cats=`ls $SRC/*.pgm`

for cat in $cats
	do
	if [ ! -d $cat ]; then
		mkdir $cat
	fi
	echo "processing category '" $cat "' ..."

	lines=`(cd $SRC/$cat; ls *$PGM)`

	for lin in $lines
	do
	base=`basename $lin $PGM`
	#echo $SRC/$cat/$lin

	(cat $SRC/$cat/$lin	\
		| im_filter binarize 1 	\
		| pen_fextract -f $method -w $width -x $x	\
		>$cat/$base.ufv) 2>/dev/null
	done
	done
