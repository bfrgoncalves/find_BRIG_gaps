# find_BRIG_gaps

# Usage 

`find_brig_gaps.py [-h] [-x QUERYTAB] [-ib SEARCHBEGIN] [-ie SEARCHEND] [-n GAPSNAME] [-s MINIMUMLENGTH] [-f REFERENCEFASTA] [-a GFFFILE] [-o RESULTSFILENAME] [-m ISMULTIPLECONTIGCOMPARISON] [-c CDSONGAPS]`

# Description 

This program searches for alignment gaps in the scratch .tab files given as output of the BRIG software. BLAST hits on those files are ordered and then regions from the reference without coverage are classified as Gaps. It returns two files, one .fasta with the gap sequences and one .gff with the coordinates in the reference.

Arguments:
 
  -h show this help message and exit

  -x QUERYTAB (Required = True)
  			.tab file from the BRIG software

  -ib SEARCHBEGIN (Required = True)
  			Start region to start looking for gaps

  -ie GAPSNAME (Required = True)
  			Name to give to the gaps found

  -n SEARCHEND (Required = True)
        End region to start looking for gaps

  
  -s MINIMUMLENGTH (Required = True) 
            Minimum gap length
  
  -f REFERENCEFASTA (Required = True)
  			.fasta file with the reference sequence to retrieve the gap sequences

  -a GFFFILE (Required = False)
        .gff file to modify with the gaps found

  -o RESULTSFILENAME (Required = True)
  			Name to give to the results files (.fasta and .gff)

  -m ISMULTIPLECONTIGCOMPARISON (Required = True)
        True if is a multiple contig comparison, false if not

  -c CDSONGAPS (Required = False)
        Writes a file with the CDS inside the gaps

# Example of usage


`find_brig_gaps.py -x queryScratchFile.tab -ib 200000 -ie 250000 -s 1000 -f referenceGenome.fasta -o results`



#Dependencies

* Biopython http://www.biopython.org

* NumPy http://www.numpy.org/

* HTSeq https://pypi.python.org/pypi/HTSeq
