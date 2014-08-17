name=		klaus
#task=		t0
# only for ICFHR - not available in original IAM-DB
task=		a

SOURCE=		<wherever your original IAM-DB data was unpacked>
LINETRANS=	ascii/lines.txt
TRAINLST=	tasks/$(task)/trainset.txt

#all:	iam-db.Corpus-t0.txt iam-db.Corpus-t0.lex
#all:	iam-db.Corpus-a.txt iam-db.Corpus-a.lex
all:	klaus.Corpus-a.txt klaus.Corpus-a.lex

#$(name).Corpus-$(task).raw:	$(SOURCE)/$(TRAINLST) $(SOURCE)/$(LINETRANS)
#	@join $(SOURCE)/$(TRAINLST) $(SOURCE)/$(LINETRANS)	\
#	| cut -d" " -f1,9-	\
#	| sed "s/ /	/"	\
#	| tr " " "_"	\
#	| tr "|" " "	\
#	| sed "s#^\([^-]*\)-\([^-]*\)#$(name)/\1/\1-\2/\1-\2#"	\
#	>$@

$(name).Corpus-$(task).raw-lex:	$(name).Corpus-$(task).raw
	cut -f2 $?	\
	| tr " " "\012"	\
	| sort -u	\
	>$@

#
# WATCH OUT: conversion of punctuation symbols to words has to be applied
# 	twice to also cover adjacent substitutions (i.e. overlapping
# 	by a single blank :-(
# 
$(name).Corpus-$(task).txt:	$(name).Corpus-$(task).raw
	cat $?		\
	| sed -e "s/	/ /" -e 's/$$/ ;/'	\
	| sed -f punctuation2words.sed	\
	| sed -f punctuation2words.sed	\
	| sed "s/ /	/"	\
	>$@

$(name).Corpus-$(task).freq:	$(name).Corpus-$(task).txt
	cut -f2 $?	\
	| sed "s/ ;//"	\
	| tr " " "\012"	\
	| sort		\
	| uniq -c	\
	| sort -nr	\
	>$@

$(name).Corpus-$(task).lex:	$(name).Corpus-$(task).freq
	cat $?		\
	| sed "s/^  *//"	\
	| cut -d" " -f2	\
	| sort		\
	>$@
