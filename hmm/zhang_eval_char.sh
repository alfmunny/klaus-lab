#! /bin/bash

eval=Charakter.test.lst
command="mm_align -b 75 -w 15 char/Charakter.$1.cl char/Charakter.$1.state char/Charakter.model ./Charakter.swu.def ./Charakter.clex.def"

cat $eval | sed "s/\$/	<CHARACTERS> + ;/" | eval "$command"
