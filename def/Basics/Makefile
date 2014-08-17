MODULE=		Basics
#AUXMODULES=	Silence.def Human.def Noise.def
#AUXMODULES=	Silence.def
#SUBMODULES=	Phonemes.def
AUXMODULES=	Whitespace.def Punctuation.def

CLEAN=		$(SUBMODULES:.def=.lex)
REMOVE=		Total.def

all:		Total.def $(SUBMODULES:.def=.lex)

Total.def:	$(AUXMODULES) $(SUBMODULES)
	cat $^ >$@

#Phonemes.lex:	Phonemes.def
#	cpp -P $^	\
#	| grep "/[^*]*/ .="	\
#	| cut -d" " -f1	\
#	| sort -u	\
#	>$@
#	@echo "the are "`cat $@ | wc -l`" unique phonemes."

clean:
	rm -f $(CLEAN)

remove:		clean
	rm -f $(REMOVE)
