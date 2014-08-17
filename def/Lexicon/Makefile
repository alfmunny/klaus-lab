MODULE=		Lexicon

name=104_Schreiber
#task=          t0
# only for ICFHR - not available in original IAM-DB
#task=           a

CORPUS_LEXICON=	../Corpus/$(name).Corpus.lex

#SUBMODULES=	iam-ondb20k.lex
#SUBMODULES=	iam-db-a1k.lex iam-db-b1k.lex
SUBMODULES=	104_Schreiber.Corpus.lex
ALLMODULES=	$(SUBMODULES)
GENDEFS=	$(ALLMODULES:.lex=.def)

CLEAN=		$(GENDEFS)
REMOVE=		Total.def Total.lex

.SUFFIXES:	.raw-lex .lex .def

#all:		Total.def Total.lex
all:		Total.lex $(SUBMODULES:.lex=.def)

.lex.def:
	echo $* "%= {" >$@
	cat $<	\
	| grep -v "^#"	\
	| sed	-e "s/$$/ |/"	\
	| sed	-f escape-chars.sed	\
	>>$@
	echo "<space> } + ;" >>$@

.raw-lex.lex:
	cat $<	\
	| sed -e "s/^/ /"	\
		-e "s/$$/ /"	\
	| sed -f punctuation2words.sed	\
	| tr -d " "	\
	| sort -u	\
	>$@

#Total.def:	$(PATCHMODULES) $(ALLMODULES:.lex=.def)
#	cat $^ >$@

Total.lex:	$(ALLMODULES) $(CORPUS_LEXICON)
	@cat $^	\
	| grep -v "^#"	\
	| sort -u	\
	| sed '/^<.*>$$/d'	\
	>$@
	@echo "the overall lexicon contains "`cat $@ | wc -l`" unique words."

clean:
	rm -f $(CLEAN)

remove:		clean
	rm -f $(REMOVE)
