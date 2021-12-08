.PHONY: all solution clean

all: solution clean

solution:
	pdflatex --shell-escape qmeesm.tex
	pdflatex --shell-escape qmeesm.tex
	pdflatex --shell-escape qmeesm.tex
	mv qmeesm.pdf ./obj/

clean:
	rm -f *.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.log
	rm -f *.out
	rm -f *.toc
	rm -f *.lot
	rm -f *.lof
