All: copynumber crlmmDownstream clean

copynumber: copynumber.tex
	cp -p ../scripts/copynumber.pdf .
	cp -p ../scripts/illumina_copynumber.pdf .

crlmmDownstream: crlmmDownstream.tex
	cp -p ../scripts/crlmmDownstream.pdf .

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
