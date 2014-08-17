#! /bin/bash

result=$1
num=$2

paste -d ' ' klaus.word_12_$num.eval $result \
	| awk '{if($2==$4)right++}END{print right/NR}'
