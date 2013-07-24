# Makefile for LaTeX files

LATEX	= pdflatex

SRC	= SparseM.tex

all 	: SparseM.pdf

%.pdf:	%.tex
	$(LATEX) $*
	rm -f *Rplots* *.ps $*-* $*.dvi $*.log $*.bbl *.tex $*.aux $*.blg $*.log
