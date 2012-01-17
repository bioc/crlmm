All: vignettes clean

## all pdfs must be generated from the first target of the Makefile
vignettes: AffyGW.tex
	cp ../scripts/AffyGW.pdf .
	cp ../scripts/IlluminaPreprocessCN.pdf .
	cp ../scripts/Infrastructure.pdf .
	cp ../scripts/crlmmDownstream.pdf .
	texi2dvi --pdf genotyping.tex
	texi2dvi --pdf CopyNumberOverview.tex

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
