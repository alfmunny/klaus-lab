#! /bin/bash

echo '#' >testlist.tmp
for word in `cat word.a`
do
	grep "/"$word"_" klaus.Corpus.lst>>testlist.tmp
done

grep -v '#' testlist.tmp>klaus.test_85.lst
