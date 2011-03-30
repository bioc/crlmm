All: copynumber crlmmDownstream clean

copynumber: copynumber.tex
	cp -p ../scripts/copynumber.pdf .
	cp -p ../scripts/AffymetrixPreprocessCN.pdf .
	cp -p ../scripts/IlluminaPreprocessCN.pdf .
	cp -p ../scripts/Infrastructure.pdf .

crlmmDownstream: crlmmDownstream.tex
	cp -p ../scripts/crlmmDownstream.pdf .

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
