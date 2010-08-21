## Need to submit from a fat node!
all:	illumina_copynumber, copynumber

copynumber:	copynumber.Rnw
	echo "Stangle(\"copynumber.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/copynumber.R/" > copynumber.R.sh
	qsub -m e -r y -cwd -l mem_free=10G copynumber.R.sh

texcopynumber:	copynumber.tex
	texi2dvi --pdf --clean --quiet copynumber.tex

illumina_copynumber:	illumina_copynumber.Rnw
	echo "Stangle(\"illumina_copynumber.Rnw\")" | R --no-save --no-restore;	
	cat ~/bin/cluster.template | perl -pe "s/Rprog/illumina_copynumber.R/" > illumina_copynumber.R.sh
	qsub -m e -r y -cwd -l mem_free=10G illumina_copynumber.R.sh

texillumina_copynumber:	illumina_copynumber.tex
	texi2dvi --pdf --clean --quiet illumina_copynumber.tex		

clean:
	-rm *.tex *.aux *.sh *.Rout
