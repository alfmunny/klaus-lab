#! /bin/bash

#cat Character.Corpus.lst \
cat Character.Corpus_all.lst \
| sed 's/BASchreiber.*/&\t<space> Bahn <space> ;/' \
| sed 's/CSchreiber.*/&\t<space> Corpus <space> ;/' \
| sed 's/daSchreiber.*/&\t<space> das <space> ;/' \
| sed 's/DSchreiber.*/&\t<space> Delicti <space> ;/' \
| sed 's/eSchreiber.*/&\t<space> Ende <space> ;/' \
| sed 's/GSchreiber.*/&\t<space> Glueck <space> ;/' \
| sed 's/ibSchreiber.*/&\t<space> Imbisswagen <space> ;/' \
| sed 's/imSchreiber.*/&\t<space> im <space> ;/' \
| sed 's/MSchreiber.*/&\t<space> Mahlzeiten <space> ;/' \
| sed 's/sSchreiber.*/&\t<space> seinem <space> ;/' \
| sed 's/USchreiber.*/&\t<space> Unterhaltung <space> ;/' \
| sed 's/vSchreiber.*/&\t<space> Vorfalls <space> ;/' \
| sed '/^$/d' \
>Character.Corpus_all.train
#>Character.Corpus.train
