#! /bin/sh

#this script is to make the path in the lists matched to a different $HOME

name=`echo $HOME | cut -d/ -f3`
echo "username is:"
echo $name

sed -i "s/matt/"$name"/" klaus.test_12_110.lst
sed -i "s/matt/"$name"/" klaus.Corpus_part_110.lst
sed -i "s/matt/"$name"/" klaus.Corpus_part_110.train
sed -i "s/matt/"$name"/" ./word/klaus.word_12_110.eval

echo "the paths have been matched with your username."
echo "make sure the folder esmeralda/ is direct under your $HOME"
