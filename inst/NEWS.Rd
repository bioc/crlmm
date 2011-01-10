\name{NEWS}
\title{News for Package 'crlmm'}

\section{Changes in version 1.10}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Using NEWS.Rd
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
