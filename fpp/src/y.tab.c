#ifndef lint
static char const 
yyrcsid[] = "$FreeBSD: src/usr.bin/yacc/skeleton.c,v 1.28 2000/01/17 02:04:06 bde Exp $";
#endif
#include <stdlib.h>
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYLEX yylex()
#define YYEMPTY -1
#define yyclearin (yychar=(YYEMPTY))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING() (yyerrflag!=0)
static int yygrowstack();
#define YYPREFIX "yy"
#define YYERRCODE 256
#define YYQWE 257
#define YYCOLON 258
#define YYOR 259
#define YYAND 260
#define YYNOT 261
#define YYEQV 262
#define YYNEQV 263
#define YYBITOR 264
#define YYBITAND 265
#define YYBITXOR 266
#define YYBITNOT 267
#define YYEQ 268
#define YYNE 269
#define YYLT 270
#define YYGT 271
#define YYLE 272
#define YYGE 273
#define YYLS 274
#define YYRS 275
#define YYADD 276
#define YYSUB 277
#define YYMUL 278
#define YYDIV 279
#define YYREM 280
#define YYDEG 281
#define YYLPAR 282
#define YYRPAR 283
#define YYNUM 284
#define YYCOMMA 285
#define YYSTOP 286
#define YYBADLEX 287
const short yylhs[] = {                                        -1,
    0,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    2,    2,    2,    2,    2,
};
const short yylen[] = {                                         2,
    2,    3,    3,    3,    3,    3,    3,    3,    3,    3,
    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
    3,    3,    5,    3,    1,    2,    2,    2,    3,    1,
};
const short yydefred[] = {                                      0,
    0,    0,    0,    0,   30,    0,    0,   25,   27,   28,
   26,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    1,   29,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,   24,
    0,    0,
};
const short yydgoto[] = {                                       6,
    7,    8,
};
const short yysindex[] = {                                    367,
  367,  367,  367,  367,    0,    0, -225,    0,    0,    0,
    0,  266,  367,  367,  367,  367,  367,  367,  367,  367,
  367,  367,  367,  367,  367,  367,  367,  367,  367,  367,
  367,  367,  367,  367,  367,    0,    0,  295,  346,  346,
  388,  388,  406,  406,  406,  422,  422,  361,  361,  361,
  361,   57,   57,  -16,  -16, -280, -280, -280, -280,    0,
  367,  322,
};
const short yyrindex[] = {                                      0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,  233,  241,
  203,  221,  131,  141,  171,  181,  211,   45,   58,   88,
  101,  -15,   15,  -75,  -45, -195, -165, -135, -105,    0,
    0, -250,
};
const short yygindex[] = {                                      0,
   -4,    1,
};
#define YYTABLESIZE 707
const short yytable[] = {                                      12,
   34,    9,   10,   11,   35,    0,    0,   23,   38,   39,
   40,   41,   42,   43,   44,   45,   46,   47,   48,   49,
   50,   51,   52,   53,   54,   55,   56,   57,   58,   59,
   60,   13,   23,   14,   15,   23,   16,   17,   18,   19,
   20,    0,   21,   22,   23,   24,   25,   26,   27,   28,
   29,   30,   31,   32,   33,   34,   62,    0,    0,   35,
   36,    3,    3,    3,    3,    0,    3,    3,    3,    3,
    3,    0,    3,    3,    3,    3,    3,    3,    3,    3,
    3,    3,    3,    3,    3,    0,    0,    3,    0,    0,
    3,    4,    4,    4,    4,    0,    4,    4,    4,    4,
    4,    0,    4,    4,    4,    4,    4,    4,    4,    4,
    4,    4,    4,    4,    4,    0,    0,    4,    0,    0,
    4,    5,    5,    5,    5,    0,    5,    5,    5,    5,
    5,    0,    5,    5,    5,    5,    5,    5,    5,    5,
    5,    5,    5,    5,    5,    0,    0,    5,    0,    0,
    5,    2,    2,    2,    2,    0,    2,    2,    2,    2,
    2,    0,    2,    2,    2,    2,    2,    2,    2,    2,
    2,    2,    2,    2,    2,    0,    0,    2,    0,    0,
    2,    6,    6,    6,    6,    0,    6,    6,    6,    6,
    6,    0,    6,    6,    6,    6,    6,    6,    6,    6,
    6,    6,    0,    0,    0,    0,    0,    6,    0,    0,
    6,    7,    7,    7,    7,    0,    7,    7,    7,    7,
    7,    0,    7,    7,    7,    7,    7,    7,    7,    7,
    7,    7,    0,    0,    0,    0,    0,    7,    0,    0,
    7,    8,    8,    8,    8,    0,    8,    8,    8,    8,
    8,    0,    8,    8,    8,    8,    8,    8,    8,    8,
    0,   31,   32,   33,   34,    0,    0,    8,   35,    0,
    8,    9,    9,    9,    9,    0,    9,    9,    9,    9,
    9,    0,    9,    9,    9,    9,    9,    9,    9,    9,
    0,    0,    0,    0,    0,    0,    0,    9,    0,    0,
    9,   10,   10,   10,   10,    0,   10,   10,   10,   10,
   10,    0,   10,   10,   12,   12,   12,   12,    0,   12,
   12,   12,   12,   12,    0,   12,   12,   10,    0,    0,
   10,    0,   29,   30,   31,   32,   33,   34,    0,    0,
   12,   35,    0,   12,   11,   11,   11,   11,    0,   11,
   11,   11,   11,   11,    0,   11,   11,   13,   13,   13,
   13,    0,   13,   13,   13,   13,   13,    0,   13,   13,
   11,    0,    0,   11,    0,    0,    0,    0,    0,    0,
    0,    0,    0,   13,    0,    0,   13,   18,   18,   18,
   18,    0,   18,   18,   18,   18,   18,   16,   16,   16,
   16,    0,   16,   16,   16,   16,   16,    0,    0,    0,
    0,    0,    0,   18,    0,    0,   18,    0,    0,    0,
    0,    0,    0,   16,    0,    0,   16,   17,   17,   17,
   17,    0,   17,   17,   17,   17,   17,   14,   14,   14,
   14,    0,   14,   14,   14,   14,   14,    0,    0,    0,
    0,    0,    0,   17,    0,    0,   17,    0,    0,   22,
   22,   22,   22,   14,   22,   22,   14,   15,   15,   15,
   15,    0,   15,   15,   15,   15,   15,   21,   21,   21,
   21,    0,   21,   21,    0,   22,    0,    0,   22,   20,
   20,   20,   20,   15,    0,    0,   15,   19,   19,   19,
   19,    0,    0,   21,    0,    0,   21,    0,    0,    0,
    0,    0,    0,    0,    0,   20,    0,    0,   20,    0,
    0,    0,   13,   19,   14,   15,   19,   16,   17,   18,
   19,   20,    0,   21,   22,   23,   24,   25,   26,   27,
   28,   29,   30,   31,   32,   33,   34,    0,   37,    0,
   35,   13,   61,   14,   15,    0,   16,   17,   18,   19,
   20,    0,   21,   22,   23,   24,   25,   26,   27,   28,
   29,   30,   31,   32,   33,   34,    0,    0,   13,   35,
   14,   15,    0,   16,   17,   18,   19,   20,    0,   21,
   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,
   32,   33,   34,    0,    0,    0,   35,   16,   17,   18,
   19,   20,    0,   21,   22,   23,   24,   25,   26,   27,
   28,   29,   30,   31,   32,   33,   34,    1,    0,    0,
   35,    0,    0,    2,   27,   28,   29,   30,   31,   32,
   33,   34,    0,    3,    0,   35,    0,    0,    4,    0,
    5,   18,   19,   20,    0,   21,   22,   23,   24,   25,
   26,   27,   28,   29,   30,   31,   32,   33,   34,    0,
    0,    0,   35,   21,   22,   23,   24,   25,   26,   27,
   28,   29,   30,   31,   32,   33,   34,    0,    0,    0,
   35,   23,   24,   25,   26,   27,   28,   29,   30,   31,
   32,   33,   34,    0,    0,    0,   35,
};
const short yycheck[] = {                                       4,
  281,    1,    2,    3,  285,   -1,   -1,  258,   13,   14,
   15,   16,   17,   18,   19,   20,   21,   22,   23,   24,
   25,   26,   27,   28,   29,   30,   31,   32,   33,   34,
   35,  257,  283,  259,  260,  286,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,  278,  279,  280,  281,   61,   -1,   -1,  285,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,  278,  279,  280,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,  278,  279,  280,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,  278,  279,  280,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,  278,  279,  280,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,   -1,   -1,   -1,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,   -1,   -1,   -1,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
   -1,  278,  279,  280,  281,   -1,   -1,  283,  285,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,  283,   -1,   -1,
  286,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  257,  258,  259,  260,   -1,  262,
  263,  264,  265,  266,   -1,  268,  269,  283,   -1,   -1,
  286,   -1,  276,  277,  278,  279,  280,  281,   -1,   -1,
  283,  285,   -1,  286,  257,  258,  259,  260,   -1,  262,
  263,  264,  265,  266,   -1,  268,  269,  257,  258,  259,
  260,   -1,  262,  263,  264,  265,  266,   -1,  268,  269,
  283,   -1,   -1,  286,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,  283,   -1,   -1,  286,  257,  258,  259,
  260,   -1,  262,  263,  264,  265,  266,  257,  258,  259,
  260,   -1,  262,  263,  264,  265,  266,   -1,   -1,   -1,
   -1,   -1,   -1,  283,   -1,   -1,  286,   -1,   -1,   -1,
   -1,   -1,   -1,  283,   -1,   -1,  286,  257,  258,  259,
  260,   -1,  262,  263,  264,  265,  266,  257,  258,  259,
  260,   -1,  262,  263,  264,  265,  266,   -1,   -1,   -1,
   -1,   -1,   -1,  283,   -1,   -1,  286,   -1,   -1,  257,
  258,  259,  260,  283,  262,  263,  286,  257,  258,  259,
  260,   -1,  262,  263,  264,  265,  266,  257,  258,  259,
  260,   -1,  262,  263,   -1,  283,   -1,   -1,  286,  257,
  258,  259,  260,  283,   -1,   -1,  286,  257,  258,  259,
  260,   -1,   -1,  283,   -1,   -1,  286,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,  283,   -1,   -1,  286,   -1,
   -1,   -1,  257,  283,  259,  260,  286,  262,  263,  264,
  265,  266,   -1,  268,  269,  270,  271,  272,  273,  274,
  275,  276,  277,  278,  279,  280,  281,   -1,  283,   -1,
  285,  257,  258,  259,  260,   -1,  262,  263,  264,  265,
  266,   -1,  268,  269,  270,  271,  272,  273,  274,  275,
  276,  277,  278,  279,  280,  281,   -1,   -1,  257,  285,
  259,  260,   -1,  262,  263,  264,  265,  266,   -1,  268,
  269,  270,  271,  272,  273,  274,  275,  276,  277,  278,
  279,  280,  281,   -1,   -1,   -1,  285,  262,  263,  264,
  265,  266,   -1,  268,  269,  270,  271,  272,  273,  274,
  275,  276,  277,  278,  279,  280,  281,  261,   -1,   -1,
  285,   -1,   -1,  267,  274,  275,  276,  277,  278,  279,
  280,  281,   -1,  277,   -1,  285,   -1,   -1,  282,   -1,
  284,  264,  265,  266,   -1,  268,  269,  270,  271,  272,
  273,  274,  275,  276,  277,  278,  279,  280,  281,   -1,
   -1,   -1,  285,  268,  269,  270,  271,  272,  273,  274,
  275,  276,  277,  278,  279,  280,  281,   -1,   -1,   -1,
  285,  270,  271,  272,  273,  274,  275,  276,  277,  278,
  279,  280,  281,   -1,   -1,   -1,  285,
};
#define YYFINAL 6
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 287
#if YYDEBUG
const char * const yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"YYQWE","YYCOLON","YYOR","YYAND",
"YYNOT","YYEQV","YYNEQV","YYBITOR","YYBITAND","YYBITXOR","YYBITNOT","YYEQ",
"YYNE","YYLT","YYGT","YYLE","YYGE","YYLS","YYRS","YYADD","YYSUB","YYMUL",
"YYDIV","YYREM","YYDEG","YYLPAR","YYRPAR","YYNUM","YYCOMMA","YYSTOP","YYBADLEX",
};
const char * const yyrule[] = {
"$accept : S",
"S : exp YYSTOP",
"exp : exp YYDEG exp",
"exp : exp YYMUL exp",
"exp : exp YYDIV exp",
"exp : exp YYREM exp",
"exp : exp YYADD exp",
"exp : exp YYSUB exp",
"exp : exp YYLS exp",
"exp : exp YYRS exp",
"exp : exp YYLT exp",
"exp : exp YYLE exp",
"exp : exp YYGT exp",
"exp : exp YYGE exp",
"exp : exp YYEQ exp",
"exp : exp YYNE exp",
"exp : exp YYBITAND exp",
"exp : exp YYBITXOR exp",
"exp : exp YYBITOR exp",
"exp : exp YYAND exp",
"exp : exp YYOR exp",
"exp : exp YYNEQV exp",
"exp : exp YYEQV exp",
"exp : exp YYQWE exp YYCOLON exp",
"exp : exp YYCOMMA exp",
"exp : term",
"term : YYSUB term",
"term : YYNOT term",
"term : YYBITNOT term",
"term : YYLPAR exp YYRPAR",
"term : YYNUM",
};
#endif
#ifndef YYSTYPE
typedef int YYSTYPE;
#endif
#if YYDEBUG
#include <stdio.h>
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH 10000
#endif
#endif
#define YYINITSTACKSIZE 200
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short *yyss;
short *yysslim;
YYSTYPE *yyvs;
int yystacksize;
#line 80 "eval.y"
#include "fpp.h"
#include "symtab.h"
#include "rgram.h"
#include "service.h"
#include "sb.h"

void
yyerror( const char *s ) {
	fppmess(ERR_IF,s);
}

int
yylex() {
	int	lmode;
	int	n;
	char	c, *s;
	SymPtr	sym;
	Token   *tokp;

	lmode = (mmode & ~(MOD_SPACE | MOD_CONT)) | MOD_IF;
	tokp = skipbl(lmode);
	switch (tokp->token[0]) {
	case '(':	return YYLPAR;
	case ')':	return YYRPAR;
	case ',':	return YYCOMMA;
	case '%':	return YYREM;
	case '/':	return YYDIV;
	case '*':	if (tokp->token[1] == '*')
				return YYDEG;
			else	return YYMUL;
	case '-':	return YYSUB;
	case '+':	return YYADD;
	case '~':	return YYBITNOT;
	case '^':	return YYBITXOR;
	case '?':	return YYQWE;
	case ':':	return YYCOLON;
	case '|':	
		if (tokp->token[1] == '|')
			return YYOR;
		else	return YYBITOR;
	case '&':
		if (tokp->token[1] == '&')
			return YYAND;
		else	return YYBITAND;
	case '!':
		if (tokp->token[1] == '=')
			return YYNE;
		else	return YYNOT;
	case '=':
		if (tokp->token[1] == '=')
			return YYEQ;
		else 	return YYSTOP;
	case '<':
		if (tokp->token[1] == '=')
			return YYLE;
		else if (tokp->token[1] == '<') 
			return YYLS;
		else 	return YYLT;
	case '>':
		if (tokp->token[1] == '=')
			return YYGE;
		else if (tokp->token[1] == '>') 
			return YYRS;
		else 	return YYGT;
	case '\n':	return YYSTOP;
	case '.':
		skipbl(lmode);
		sym = symget(tokp->token,CL_FOP);
		if (sym == NULL) return YYBADLEX;
		n = symvali(sym);
		switch (n) {
		case FTN_TRUE: n = YYNUM; yylval = 1;
			break;
		case FTN_FALSE: n = YYNUM; yylval = 0;
			break;
		case FTN_EQ : n = YYEQ;
			break;
		case FTN_NE : n = YYNE;
			break;
		case FTN_LT : n = YYLT;
			break;
		case FTN_LE : n = YYLE;
			break;
		case FTN_GT : n = YYGT;
			break;
		case FTN_GE : n = YYGE;
			break;
		case FTN_AND : n = YYAND;
			break;
		case FTN_OR : n = YYOR;
			break;
		case FTN_NEQV :
		case FTN_XOR : n = YYNEQV;	/* these two are the same */
			break;
		case FTN_EQV : n = YYEQV;
			break;
		case FTN_NOT : n = YYNOT;
			break;
		default:
			return YYBADLEX;
		}
		skipbl(lmode);
		if (tokp->token[0] != '.') 
			return YYBADLEX;
		else	return n;
	default: 
		if (tokp->val == TK_NAME) {
			if (!strcmp(tokp->token,"defined")) {
				int save;

				save = dosubstfl;
				dosubstfl = 0;
				tokp = skipbl(lmode);
				if (tokp->val == TK_NAME) {
					dosubstfl = save;
					if (symget(tokp->token,CL_NM))
						yylval = 1;
					else
						yylval = 0;
					return YYNUM;
				}	
				else if (tokp->token[0] != '(') {
					dosubstfl = save;
					return YYBADLEX;
				}
				tokp = skipbl(lmode);
				dosubstfl = save;
				if (tokp->val == TK_NAME && symget(tokp->token,CL_NM))
					yylval = 1;
				else
					yylval = 0;
				tokp = skipbl(lmode);
				if (tokp->token[0] != ')') 
					return YYBADLEX;
				return YYNUM;
			}
			else {
				if (sbfl) {
					sb_mref(tokp->token, 0, 0);
					sb_mref_end();
				}
				yylval = 0;
				return YYNUM;
			}
		}
		else if (tokp->val == TK_NUM) {
			s = tokp->token;
			while (c=*s++) {
				if (!is_num(c))
					return YYBADLEX;
			}
#if USE_C_HEX_CONST
			if (tokp->token[0] == '0') {
				/* C octal constant is allowed 
				 * in #if expression */
				strtoi(tokp->token,&yylval,8);
			}
			else
#endif /* USE_C_HEX_CONST */
				strtoi(tokp->token,&yylval,10);
			return YYNUM;
		}
		else if (tokp->val == TK_BOZ) {
			int	err;

			switch (lowcase(tokp->token[0])) {
			case 'b':
				tokp = skipbl(lmode);
				CHECK(tokp->token[0] == '\'' || tokp->token[0] == '"');

				/* Erase the trailing quote */
				tokp->token[--tokp->length] = 0;
				err = strtoi(tokp->token+1, &yylval, 2);
				if (err) 
					return YYBADLEX;
				else	return YYNUM;
				break;
			case 'o':
				tokp = skipbl(lmode);
				CHECK(tokp->token[0] == '\'' || tokp->token[0] == '"');
				tokp->token[--tokp->length] = 0;
				err = strtoi(tokp->token+1, &yylval, 8);
				if (err) 
					return YYBADLEX;
				else	return YYNUM;
				break;
			case 'x':
			case 'z':
				tokp = skipbl(lmode);
				CHECK(tokp->token[0] == '\'' || tokp->token[0] == '"');
				tokp->token[--tokp->length] = 0;
				err = strtoi(tokp->token+1, &yylval, 16);
				if (err) 
					return YYBADLEX;
				else	return YYNUM;
				break;
			case '\'':
			case '"':
				c = tokp->token[tokp->length - 1];
				if (lowcase(c) == 'o') {
					tokp->length -= 2;
					tokp->token[tokp->length] = 0;
					err = strtoi(tokp->token+1, &yylval, 8);
				}
				else if (lowcase(c) == 'x') {
					tokp->length -= 2;
					tokp->token[tokp->length] = 0;
					err = strtoi(tokp->token+1, &yylval, 16);
				}
				else return YYBADLEX;
				if (err) 
					return YYBADLEX;
				else	return YYNUM;
#if USE_C_HEX_CONST
			case '0':
				CHECK(lowcase(tokp->token[1]) == 'x');
				err = strtoi(tokp->token+2, &yylval, 16);
				if (err) 
					return YYBADLEX;
				else	return YYNUM;
				break;
#endif /* USE_C_HEX_CONST */
			default:
				CHECK(0);
				return YYBADLEX;
			}
		}
		else
			return YYBADLEX;
	}
/*	return YYBADLEX; */
}
#line 554 "y.tab.c"
/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack()
{
    int newsize, i;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = yystacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;
    i = yyssp - yyss;
    newss = yyss ? (short *)realloc(yyss, newsize * sizeof *newss) :
      (short *)malloc(newsize * sizeof *newss);
    if (newss == NULL)
        return -1;
    yyss = newss;
    yyssp = newss + i;
    newvs = yyvs ? (YYSTYPE *)realloc(yyvs, newsize * sizeof *newvs) :
      (YYSTYPE *)malloc(newsize * sizeof *newvs);
    if (newvs == NULL)
        return -1;
    yyvs = newvs;
    yyvsp = newvs + i;
    yystacksize = newsize;
    yysslim = yyss + newsize - 1;
    return 0;
}

#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab

#ifndef YYPARSE_PARAM
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG void
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif	/* ANSI-C/C++ */
#else	/* YYPARSE_PARAM */
#ifndef YYPARSE_PARAM_TYPE
#define YYPARSE_PARAM_TYPE void *
#endif
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG YYPARSE_PARAM_TYPE YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL YYPARSE_PARAM_TYPE YYPARSE_PARAM;
#endif	/* ANSI-C/C++ */
#endif	/* ! YYPARSE_PARAM */

int
yyparse (YYPARSE_PARAM_ARG)
    YYPARSE_PARAM_DECL
{
    register int yym, yyn, yystate;
#if YYDEBUG
    register const char *yys;

    if ((yys = getenv("YYDEBUG")))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    if (yyss == NULL && yygrowstack()) goto yyoverflow;
    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if ((yyn = yydefred[yystate])) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yysslim && yygrowstack())
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#if defined(lint) || defined(__GNUC__)
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#if defined(lint) || defined(__GNUC__)
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yysslim && yygrowstack())
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 1:
#line 34 "eval.y"
 {return yyvsp[-1];}
break;
case 2:
#line 36 "eval.y"
 { 
		int i, res = 1;
		if (yyvsp[0] < 0) res = 0;
		else {
			for (i=0; i < yyvsp[0]; i++) 
				res *= yyvsp[-2];
		}
		yyval = res; }
break;
case 3:
#line 44 "eval.y"
 { yyval = yyvsp[-2] * yyvsp[0]; }
break;
case 4:
#line 45 "eval.y"
 { 
		if (yyvsp[0] == 0) {
			fppmess(WARN_FPP_EXPR);
			yyval = 0;
		}
		else
			yyval = yyvsp[-2] / yyvsp[0]; }
break;
case 5:
#line 52 "eval.y"
 { yyval = yyvsp[-2] % yyvsp[0]; }
break;
case 6:
#line 53 "eval.y"
 { yyval = yyvsp[-2] + yyvsp[0]; }
break;
case 7:
#line 54 "eval.y"
 { yyval = yyvsp[-2] - yyvsp[0]; }
break;
case 8:
#line 55 "eval.y"
 { yyval = yyvsp[-2] << yyvsp[0]; }
break;
case 9:
#line 56 "eval.y"
 { yyval = yyvsp[-2] >> yyvsp[0]; }
break;
case 10:
#line 57 "eval.y"
 { yyval = yyvsp[-2] < yyvsp[0]; }
break;
case 11:
#line 58 "eval.y"
 { yyval = yyvsp[-2] <= yyvsp[0]; }
break;
case 12:
#line 59 "eval.y"
 { yyval = yyvsp[-2] > yyvsp[0]; }
break;
case 13:
#line 60 "eval.y"
 { yyval = yyvsp[-2] >= yyvsp[0]; }
break;
case 14:
#line 61 "eval.y"
 { yyval = yyvsp[-2] == yyvsp[0]; }
break;
case 15:
#line 62 "eval.y"
 { yyval = yyvsp[-2] != yyvsp[0]; }
break;
case 16:
#line 63 "eval.y"
 { yyval = yyvsp[-2] & yyvsp[0]; }
break;
case 17:
#line 64 "eval.y"
 { yyval = yyvsp[-2] ^ yyvsp[0]; }
break;
case 18:
#line 65 "eval.y"
 { yyval = yyvsp[-2] | yyvsp[0]; }
break;
case 19:
#line 66 "eval.y"
 { yyval = yyvsp[-2] && yyvsp[0]; }
break;
case 20:
#line 67 "eval.y"
 { yyval = yyvsp[-2] || yyvsp[0]; }
break;
case 21:
#line 68 "eval.y"
 { yyval = !yyvsp[-2] && yyvsp[0] || yyvsp[-2] && !yyvsp[0]; }
break;
case 22:
#line 69 "eval.y"
 { yyval = !yyvsp[-2] && !yyvsp[0] || yyvsp[-2] && yyvsp[0]; }
break;
case 23:
#line 70 "eval.y"
 { yyval = yyvsp[-4] ? yyvsp[-2] : yyvsp[0]; }
break;
case 24:
#line 71 "eval.y"
 { yyval = yyvsp[0]; }
break;
case 25:
#line 72 "eval.y"
 { yyval = yyvsp[0]; }
break;
case 26:
#line 74 "eval.y"
 { yyval = -yyvsp[0]; }
break;
case 27:
#line 75 "eval.y"
 { yyval = !yyvsp[0]; }
break;
case 28:
#line 76 "eval.y"
 { yyval = ~yyvsp[0]; }
break;
case 29:
#line 77 "eval.y"
 { yyval = yyvsp[-1]; }
break;
case 30:
#line 78 "eval.y"
 { yyval = yyvsp[0]; }
break;
#line 882 "y.tab.c"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yysslim && yygrowstack())
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
