### this file has been modified for gnu makefile style 

SRC             = ./
CC              = gcc
YACC            = bison -y 
COPTS           = -ansi -O -w
EXTRA_DEF       = "-DRELEASE=\"v1.1\""
BINDIR          = ../bin/

SORC = fpp.c proc.c lex.c error.c rgram.c symtab.c service.c 
OBJ  = fpp.o proc.o lex.o error.o rgram.o symtab.o service.o
OBJ1 = y.tab.o


ALL_OPTS = $(COPTS) $(EXTRA_DEF) $(EXTRA_INC)

all:    $(BINDIR)/fpp
	

$(BINDIR)/fpp:	$(OBJ) $(OBJ1) $(EXTRA_LIB)
	$(CC) $(COPTS) -o $@ $(OBJ) $(OBJ1) $(EXTRA_LIB)

$(OBJ):	$(SORC)
	$(CC) -c $(ALL_OPTS) $(SORC)

#y.tab.o: y.tab.c
y.tab.o: 
	$(CC) -c $(ALL_OPTS) y.tab.c

#y.tab.c: $(SRC)/eval.y
#	$(YACC) $(SRC)/eval.y


clean:
	$(RM) fpp *.o $(BINDIR)/fpp
#	$(RM) fpp *.o y.tab.c $(BINDIR)/fpp

fpp.o:		$(SRC)/fpp.h 
fpp.o:		$(SRC)/service.h
service.o:	$(SRC)/fpp.h
service.o:	$(SRC)/service.h
proc.o:		$(SRC)/fpp.h
proc.o:		$(SRC)/service.h
proc.o:		$(SRC)/lex.h
proc.o:		$(SRC)/symtab.h
proc.o:		$(SRC)/rgram.h
proc.o:		$(SRC)/sb.h
lex.o:		$(SRC)/fpp.h
lex.o:		$(SRC)/service.h
lex.o:		$(SRC)/symtab.h
lex.o:		$(SRC)/rgram.h
lex.o:		$(SRC)/lex.h
lex.o:		$(SRC)/sb.h
symtab.o:	$(SRC)/fpp.h
symtab.o:	$(SRC)/service.h
symtab.o:	$(SRC)/symtab.h
error.o:	$(SRC)/fpp.h
y.tab.o:	$(SRC)/fpp.h
y.tab.o:	$(SRC)/symtab.h
y.tab.o:	$(SRC)/rgram.h
y.tab.o:	$(SRC)/service.h
y.tab.o:	$(SRC)/sb.h
rgram.o:	$(SRC)/fpp.h
rgram.o:	$(SRC)/rgram.h
