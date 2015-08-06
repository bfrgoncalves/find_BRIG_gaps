# find_BRIG_gaps

# Usage 

`find_brig_gaps.py [-h] [-x QUERYTAB] [-ib SEARCHBEGIN] [-ie SEARCHEND] [-s MINIMUMLENGTH] [-f REFERENCEFASTA] [-o RESULTSFILENAME]`

# Description 

This program searches for alignment gaps in the scratch .tab files given as output of the BRIG software. BLAST hits on those files are ordered and then regions from the reference without coverage are classified as Gaps.

Arguments:
 
  -h show this help message and exit

  -x QUERYTAB (Required = True)
  			.tab file from the BRIG software

  -ib SEARCHBEGIN (Required = True)
  			Start region to start looking for gaps

  -ie SEARCHEND (Required = True)
  			End region to start looking for gaps
  
  -s MINIMUMLENGTH (Required = True) 
            Minimum gap length
  
  -f REFERENCEFASTA (Required = True)
  			.fasta file with the reference sequence to retrieve the gap sequences

  -o RESULTSFILENAME (Required = True)
  			Name to give to the results file

# Example of usage


`find_brig_gaps.py -x queryScratchFile.tab -ib 200000 -ie 250000 -s 1000 -f referenceGenome.fasta -o results.txt`



#Dependencies

* Biopython http://www.biopython.org

* NumPy http://www.numpy.org/

* HTSeq https://pypi.python.org/pypi/HTSeq
