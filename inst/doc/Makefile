All: copynumber illumina_copynumber crlmmDownstream clean

copynumber: copynumber.tex
#	cp  -p ../scripts/copynumber.pdf .
#	cp  -p ../scripts/illumina_copynumber.pdf .
	cp  ../scripts/copynumber.pdf .
	cp  ../scripts/illumina_copynumber.pdf .

illumina_copynumber: illumina_copynumber.tex
	cp ../scripts/illumina_copynumber.pdf .

crlmmDownstream: crlmmDownstream.tex
	cp -p ../scripts/crlmmDownstream.pdf .

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
