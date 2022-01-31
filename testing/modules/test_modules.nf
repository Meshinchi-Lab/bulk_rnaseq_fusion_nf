process echo_file  {
  input:
  path infile

  output:
  stdout 

  script:
  """
  echo $infile
  """
}

//Helper function to unzip files when needed 
process unzip {

  input:
	  path zipped_file  

	output:
	  path "*", emit: unzipped_file

  when: 
    check_gz != 'EMPTY'

	script:
    """
    gunzip -f $zipped_file 
    """
}