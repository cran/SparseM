# Makefile for LaTeX files

LATEX	= latex
BIBTEX	= bibtex
DVIPS	= dvips
PSTOPDF = ps2pdf

SRC	= SparseM.tex

all 	: SparseM.pdf

%.pdf:	%.tex
	$(LATEX) $*
	$(DVIPS) $*
	$(PSTOPDF) $*.ps
	rm -f *.ps $*-* $*.dvi $*.log $*.bbl *.tex $*.aux $*.blg $*.log
