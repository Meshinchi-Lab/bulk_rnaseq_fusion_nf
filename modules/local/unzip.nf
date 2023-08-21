process UNZIP {

  input:
	  path zipped_file  

	output:
	  path "*", emit: UNZIPped_file

	script:
    """
    gUNZIP -f $zipped_file 
    """
}
