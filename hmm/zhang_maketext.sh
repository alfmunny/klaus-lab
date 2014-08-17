#! /bin/bash

cat Charakter_train_pgm.lst \
| sed 's/BASchreiber.*/&\tBahn ;/' \
| sed 's/CSchreiber.*/&\tCorpus ;/' \
| sed 's/daSchreiber.*/&\tdas ;/' \
| sed 's/DSchreiber.*/&\tDelicti ;/' \
| sed 's/eSchreiber.*/&\tEnde ;/' \
| sed 's/GSchreiber.*/&\tGlueck ;/' \
| sed 's/ibSchreiber.*/&\tImbisswagen ;/' \
| sed 's/imSchreiber.*/&\tim ;/' \
| sed 's/MSchreiber.*/&\tMahlzeiten ;/' \
| sed 's/sSchreiber.*/&\tseinem ;/' \
| sed 's/USchreiber.*/&\tUnterhaltung ;/' \
| sed 's/vSchreiber.*/&\tVorfalls ;/' \
| sed '/^$/d' \
>Charakter_train_pgm.txt

