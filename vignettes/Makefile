UNAME := $(shell uname)
All: vignettes clean

vignettes: AffyGW.tex
		cp ../inst/scripts/AffyGW.pdf .
		cp ../inst/scripts/IlluminaPreprocessCN.pdf .
		cp ../inst/scripts/Infrastructure.pdf .
		cp ../inst/scripts/gtypeDownstream.pdf .
ifeq ("$(UNAME)", "CYGWIN_NT-6.1-WOW64")
		cmd /C texi2dvi --pdf genotyping.tex
		cmd /C texi2dvi --pdf CopyNumberOverview.tex
else
		texi2dvi --pdf genotyping.tex
		texi2dvi --pdf CopyNumberOverview.tex
endif

clean:
		-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex