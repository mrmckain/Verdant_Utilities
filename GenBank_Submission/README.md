# GenBank Submission

<h1>Protocol for GenBank Submission Ready Files</h1>

The following steps can be used to quickly verify annotations and prepare GenBank upload files. These steps can be done on either a Linux distribution or a Mac.

<h2>Download tbl2asn program</h2>

<a href="https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/">tbl2asn</a> can be used to automate creation of sequence records for GenBank submission. It can be downloaded <a href="ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/">here</a>.

<h2>Create a submission template file</h2>

<a href="https://submit.ncbi.nlm.nih.gov/genbank/template/submission/">Create</a> GenBank submission template (a .sbt file) to be used for the project.  An .sbt file must be used with tbl2asn. If you are only checking annotations, you can create a dummy file with incomplete information. 

<h2>Download TBL files for project</h2>

On your Verdant project page, click the "DOWNLOAD PROJECT TBL DATA" button.  

<br>

<img src="https://github.com/mrmckain/Verdant_Utilities/blob/master/GenBank_Submission/images/Verdant_TBL.png" width="384" alt="verdant_tbl" align="middle">

<br>

Unzip the file, and navigate into the directory using a terminal.

<h2>Creating the upload files</h2>

Run the following commands from inside the tbl file directory downloaded from Verdant:

	mkdir all_files
	
	mkdir upload_files

	perl fix_fasta_header_for_tbl2asn.pl [This script in found in this directory.]

	mv *fsa all_files/.

	cp tbl_files/*tbl all_files/.

	tbl2asn -t YOURFILE.sbt -p all_files -V vb -r upload_files

These commands will create a directory called "upload_files" that contains *.val, *.sqn, and *.gbf files for each sample.  To verify protein coding genes from the annotation, look at the *.val files.  If there are problems (e.g. termination symbols in protein sequence, gap symbol at start of protein sequence), these files will tell you.  Termination symbols in protein sequence can either be an annotation issue or indication that this gene is pseudgenized or misassembled. It should be verified. Gap symbols at start of protein sequence are usually the result of RNA-edited genes.  For angiosperms, these are commonly rps19 and rpl2, though others exist for different lineages. These should be verifed as well.

If you are satisfied with the annotations, the "upload"files" can be submitted to GenBank directly.



