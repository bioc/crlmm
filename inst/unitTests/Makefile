TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
R="${R_HOME}/bin/R"

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	R CMD INSTALL --no-lock ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R --vanilla --slave < ${SUITE}

release: # Run unit tests
	cd ${TOP}/..;\
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/2.15-bioc-release/library R --vanilla --slave < ${SUITE}