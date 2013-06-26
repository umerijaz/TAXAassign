TAXAassign
==========

 Script to annotate sequences at different taxonomic levels using  NCBI's taxonomy

* The taxonomic assignment is resolved using NCBI’s Taxonomy and running NCBI’s Blast against locally-installed NCBI’s nt database to minimize execution time.
* Version 0.3 has many orders of magnitude improvement in speed over 0.1. 
* To minimize the execution time, we use GNU Parallel, a shell tool for executing jobs in parallel on multicore computers. We split the sequence file into fixed size chunks and then run blastn in parallel on these chunks on separate cores. For a 16SrRNA dataset comprising 1000 most abundant OTU sequences, matching at most 100 reference sequences against a local NCBI’s NT database took 18.9 minutes on 45 cores. A speedup of 30 times or more is achieved this way.
* To find Genbank ID to Taxa ID mapping, one can use gi_taxid_nucl.dmp.gz from ftp://ftp.ncbi.nih.gov/pub/taxonomy/  which is updated every Monday around 2am EST on NCBI’s website. However, searching Taxa IDs through this beastly 4.4GB+ unzipped file is time consuming and slows down the pipeline as noticed in the previous version. To improve the performance, instead of using this file, we use the latest version of NCBI blast 2.28+ which also gives Taxa ID for each hit and so we skip this step altogether in this version.
* To improve the time for finding parent Taxa IDs for a given Taxa ID, instead of using the flat file taxdump.tar.gz from the above ftp site, we use BioSQL and host NCBI’s taxonomy data in a local MySQL server. Furthermore, we use GNU Parallel to run multiple SQL queries in parallel on multiple records of blast output file thus reducing the execution time significantly
* Both NCBI’s taxonomy database on MySQL server/sqlite3 and local nt database can be updated frequently by submitting a cron job on the server scheduled to run when the server is less busy i.e. at night time and thus the information does not get outdated.
* All the time consuming steps in the TAXAassign pipeline are not repeated on reruns. Thus, if pipeline has finished processing the data and you want to run it again with a different taxonomic level, it will skip blasting the fasta file. To see what you have already done, output gets logged in TAXAassign.log file in the current folder. This is useful for debugging, should a problem arise.
* If you have the OTUs abundance table (in csv or tsv format) along with OTUs sequences, the output file generated from TAXAassign can then be used withhttp://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/convIDs.pl to annotate the table.
* For better understanding of BioSQL, refer to my tutorialhttp://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/BIOSQL_tutorial.pdf 
