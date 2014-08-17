#! /bin/bash

paste -d ' ' Character.word.eval Character.Xeval-a.a1k-$1-b75-w15-l7-B20.words \
	| awk '{if($2==$4)right++}END{print right/NR}'
