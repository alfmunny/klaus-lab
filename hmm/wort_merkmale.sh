#! /bin/bash

#
# NOTE: tools are found in $(ESMERALDA)/bin
#

SRC=$ESMERALDA/lab/12_Wort/lin 	#bearbeitete Bilder
DES=$ESMERALDA/lab/12_Wort/ufv
PGM="_slope_slant_scale.pgm" #bearbeitete Bilder
#PGM=".pgm"
UFV="_sss_ex.ufv"

#
# frame extraction parameters
#
method=2	# feature extraction method 
width=$1	# frame width
x=$2		# frame overlap :-(

cats=`ls $SRC/*_scale.pgm` #bearbeitete Bilder
#cats=`ls $SRC/*.pgm`

for lin in $cats
	do
	echo $lin
	base=`basename $lin $PGM`
	#echo $base
	pen_fextract -f $method -w $width -x $x $lin $DES/$base$UFV
	done 
