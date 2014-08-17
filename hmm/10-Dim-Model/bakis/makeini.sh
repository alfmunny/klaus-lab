#! /bin/bash
# die verfügbare .ini zu erzeugen

hmm=$1 # mit welcher HMM
name=klaus

paste -d '\t' ../../$name.Corpus.lst ../linear/$name.Corpus-$hmm.hyp \
| sed 's/<space1>/[space]/g' \
| sed 's/<space2>/[space]/g' \
| sed 's/<space3>/[space]/g' \
| grep -v '?' \
>$name.Corpus-$hmm.ini
