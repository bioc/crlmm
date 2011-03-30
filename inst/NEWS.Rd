\name{NEWS}
\title{News for Package 'crlmm'}


\section{Changes in version 1.9}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Using NEWS.Rd

      \item batch slot in CNSet objects is class 'character'
      (previously, class was 'factor')

      \item the ff package is required for preprocessing and genotyping
      prior to copy number analyses with the crlmmCopynumber function.

      \item We have added several vignettes pertaining to copy number
      analyses in the crlmm package: CopyNumberOverview,
      AffymetrixPreprocessCN, IlluminaPreprocessCN, and Infrastructure.
      The 'copynumber' and 'Infrastructure' vignettes are applicable to
      both the Illumina and Affymetrix platforms. The CopyNumberOverview
      vignette gives a brief summary of the available vignettes for copy
      number analyses.

      \item Exporting function constructInf, preprocessInf, and
      genotypeInf for preprocessing and genotyping Illumina files prior
      to copy number analyses.  The genotype.Illumina function is now a
      wrapper for these functions.

      \item additional documentation for crlmm is provided in a
      compendium: http://www.biostat.jhsph.edu/~rscharpf/crlmmCompendium/index.html

    }
  }
}

\section{Changes in version 1.3}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item 3 new classes created:
      \itemize{
	\item 'ABset': container for quantile-normalized A and B
       	intensities for both SNP and copy number probes.  Required
       	assay data elements are 'A' and 'B'.  Extends eSet
       	directly.
	\enumerate{
          \item For nonpolymorphic probes, the quantile normalized
          intensity is stored in the 'A' assay data element.
          The corresponding row in the 'B' assay data element
          is NA.  This is a bit inefficient, but greatly
          simplifies downstream analyses.  In particular, '['
          works.
	}

        \item 'CrlmmSetList': container for results from preprocessing
	and genotyping.  This object is a list.  The first
	element of the list is an ABset.  The second element is a
	SnpSet containing genotype calls.  The two elements are
	required to have identical featureNames and sampleNames.
	\enumerate{
       	  \item added several methods for subsetting and accessing
 	  elements of this object, including featureNames,
 	  sampleNames, and '['.
	}

        \item 'CopyNumberSet': contains locus-level estimates of copy
  	number for SNPs and polymorphic probes.
	\enumerate{
	  \item Required assay data elements are 'CA' and 'CB',
          corresponding to the absolute copy number for allele
          A and B, respectively.

	  \item For nonpolymorphic probes, the total copy number is
          stored in the 'CA' slot and a NA is recorded for the
          corresponding row in the CB matrix.

	  \item Useful methods: 'copyNumber', 'ellipse', 'points'
	}
      }
      \item 'crlmmWrapper' function does preprocessing
      (quantile-normalization) and genotyping, saving an object of
      class CrlmmSetList for each chromosome

     \item 'computeCopynumber' now requires an object of class
     'CrlmmSetList' and returns an object of class 'CopyNumberSet'.
    }
  }

}
