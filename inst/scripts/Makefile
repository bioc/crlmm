## Need to submit from a fat node!
all:	illumina affy 

cleanScripts:
	-rm *.pdf *.tex *.log *.aux 
	-rm *.bbl *.bbg
	-rm auto

cleanIllumina:
	-rm -rf /thumper/ctsa/snpmicroarray/rs/ProcessedData/crlmm/trunk/illumina_vignette
	mkdir /thumper/ctsa/snpmicroarray/rs/ProcessedData/crlmm/trunk/illumina_vignette

cleanAffy:
	-rm -rf /thumper/ctsa/snpmicroarray/rs/ProcessedData/crlmm/trunk/affy_vignette
	mkdir /thumper/ctsa/snpmicroarray/rs/ProcessedData/crlmm/trunk/affy_vignette

testing: testing.Rnw
	echo "Stangle(\"testing.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.templateRelease | perl -pe "s/Rprog/testing.R/" > testing.R.sh
	qsub -m e -r y -cwd -l mem_free=12G,h_vmem=16G testing.R.sh

testingFF: testingFF.Rnw
	echo "Stangle(\"testingFF.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.templateRelease | perl -pe "s/Rprog/testingFF.R/" > testingFF.R.sh
	qsub -m e -r y -cwd -l mem_free=12G,h_vmem=16G testingFF.R.sh

affy:	copynumber.Rnw
	echo "Stangle(\"copynumber.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.templateRelease | perl -pe "s/Rprog/copynumber.R/" > copynumber.R.sh
	qsub -m e -r y -cwd -l mem_free=12G,h_vmem=16G copynumber.R.sh

illumina:	illumina_copynumber.Rnw
	echo "Stangle(\"illumina_copynumber.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.templateRelease | perl -pe "s/Rprog/illumina_copynumber.R/" > illumina_copynumber.R.sh
	qsub -m e -r y -cwd -l mem_free=12G,h_vmem=16G illumina_copynumber.R.sh

texcopynumber:	copynumber.tex
	texi2dvi --pdf --clean --quiet copynumber.tex

texillumina_copynumber:	illumina_copynumber.tex
	texi2dvi --pdf --clean --quiet illumina_copynumber.tex		

