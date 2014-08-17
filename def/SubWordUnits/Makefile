MODULE=		SubwordUnits
LANGUAGE=	German

BASICS=		../Basics/Total.def
#LEXICON=	../Lexicon/Total.lex
LEXICON=	../Lexicon/klaus-word_12.lex

#CORPUSFREQ=	../Corpus/Corpus.freq

#PMOD_OPTS=	-l1 -r1 -w
PMOD_OPTS=	-l0 -r0
PMOD_OPTS+=	-B ""	# no word-internal boundaries
PMOD_CHARS=	$(LANGUAGE).chars
PUNCTMAP=	punctuation.sed
ESCAPECHARS=	escape-chars.sed

#GENLEVELS=	SubwordUnits.uni.def SubwordUnits.bi.def
#GENLEVELS=	SubwordUnits.uni.def

CLEAN=		SubwordUnits.def SubwordUnits.lex \
		$(GENLEVELS) $(GENLEVELS:.def=.lex)
REMOVE=		Total.def $(LANGUAGE).chars

.SUFFIXES:	.lex .def

all:		$(LANGUAGE).chars Total.def

.def.lex:
	cut -f2 $<	\
	| sed	-e "s/ ;//" \
	| tr " " "\012"	\
	| sort -u	\
	>$@

#Total.def:	SubwordUnits.def $(GENLEVELS) $(BASICS) $(UNTRANSCRIBED)
Total.def:	SubwordUnits.def $(GENLEVELS) $(BASICS)
	cat $(BASICS) >$@
	cat $(GENLEVELS) SubwordUnits.def >>$@
#	cat $(UNTRANSCRIBED)	\
#	| sed "s/$$/ := <anyspeech> ;/"	\
#	>>$@

#
# ... wenn alle Woerter gleich behandelt werden
#
SubwordUnits.def:	$(LEXICON) $(PMOD_CHARS) $(PUNCTMAP) $(ESCAPECHARS)
	paste $(LEXICON) $(LEXICON)	\
	| grep -v "^<"	\
	| grep -v "^@"	\
	| mm_pmod $(PMOD_OPTS) $(PMOD_CHARS) \
	| sed -f $(PUNCTMAP)	\
	| sed -f $(ESCAPECHARS)	\
	> $@

SubwordUnits.uni.def:	SubwordUnits.lex
	grep "..*/.*/" $^ \
	| sed "s/^\(.*\)\(\/.*\/\)\(.*\)/\1\2\3 ==	\2 ;/" \
	> $@
	
#
# Characters valid for English words
#
German.chars:	letters.chars digits.chars punctuation.chars
	cat $^	\
	>$@


clean:
	rm -f $(CLEAN)

remove:		clean
	rm -f $(REMOVE)
