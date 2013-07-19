TAXAassign v0.4
===============

TAXAassign is useful for annotating nucleotide sequences (contigs from assemblies, reads from whole-shot gun sequencing, 16S rRNA sequences, etc.) at different taxonomic levels (Phylum, Class, Order, Family, Genus, and Species) using NCBI's Taxonomy. 
The first public release (version 0.3) was developed and completed by [Umer Zeeshan Ijaz](http://userweb.eng.gla.ac.uk/umer.ijaz/index.htm) under the supervision of [Christopher Quince](http://userweb.eng.gla.ac.uk/christopher.quince) in the following hackathon:

**Event title:** ProBin: Probabilistic binning for metagenome contigs  
**Location:** Instituto Gulbenkian De CiêNcia, Lisbon, Portugal  
**Dates:** from 24-06-2013 to 28-06-2013  

organized by **European Union's Earth System Science and Environmental Management ES1103 COST Action** ("[Microbial ecology & the earth system: collaborating for insight and success with the new generation of sequencing tools](http://www.cost.eu/domains_actions/essem/Actions/ES1103)").  
This work would not have been possible, were it not for the useful discussions with other participants of the hackathon, namely,  

**[Nick Loman](http://pathogenomics.bham.ac.uk/blog/author/nick/)**  
**[Joshua Quick](http://pathogenomics.bham.ac.uk/clinicogenomics/)**  
**[Brynjar Smari Bjarnason](http://is.linkedin.com/pub/brynjar-sm%C3%A1ri-bjarnason/b/a7b/964)**  
**[Johannes Alneberg](http://se.linkedin.com/pub/johannes-alneberg/49/bab/84)**  
**[Anders Anderson](http://www.scilifelab.se/index.php?content=research_groups&id=2)**  
**[Ino de Bruijn](http://se.linkedin.com/in/deknappeinodebruijn)**  

There are several softwares out there that perform similar classification of high throughput sequencing and do a very good job, e.g., [MEGAN](http://ab.inf.uni-tuebingen.de/software/megan/), [CREST](http://apps.cbu.uib.no/crest/index), [RDP Classifier](http://sourceforge.net/projects/rdp-classifier/files/) etc., but we believe that TAXAassign stands apart as it is simple, has ease-of-use, and offers more control over filtering out unwanted assignments.
The taxonomic assignment is resolved using NCBI’s Taxonomy and by running NCBI’s Blastn against locally-installed NCBI’s nt database. To minimize execution time, we use GNU Parallel, a shell tool for executing jobs in parallel on multicore computers. The sequence file is split into fixed size chunks and are run through Blastn in parallel on separate cores.
For example, when a given 16S rRNA dataset comprising 1000 most abundant OTU sequences and matching at most 100 reference sequences was run using GNU parallel, it took 18.9 minutes on 45 cores. 
A speedup of 30 times or more is achieved this way.  
To find Genbank ID to Taxa ID mapping, one can use **gi_taxid_nucl.dmp.gz** from ftp://ftp.ncbi.nih.gov/pub/taxonomy/ (It gets updated every Monday around 2am EST on NCBI’s website). 
However, searching Taxa IDs through this beastly 4.4GB+ unzipped file is time consuming and slows down the pipeline as noticed in the previous version. 
To improve the performance, instead of using this file, we use the latest version of NCBI's Blast 2.28+ which also gives Taxa ID for each hit in it's output
and so this step can then be skipped altogether.  
To improve the time for finding parent Taxa IDs for a given Taxa ID, instead of using the flat file **taxdump.tar.gz** from the above ftp site, we use BioSQL and host NCBI’s taxonomy data in a local MySQL server or an sqlite3 database. 
Furthermore, we use GNU Parallel to run multiple SQL queries in parallel on multiple records of Blastn output file thus reducing the execution time significantly.
Both NCBI’s taxonomy database on MySQL server or sqlite3, and local nt database can be updated frequently by submitting a cron job on the server scheduled to run when the server is less busy i.e. at night time and thus the information does not get outdated.
sqlite3 database (it can be downloaded from Umer's website mentioned below in the "Dependencies" section, or can alternatively be generated from the MySQL server by exporting the schema and populating the records) is offered as an alternative in those scenarios when the user does not have the resources to run a MySQL server on his computer.   
All the time consuming steps in the TAXAassign pipeline are not repeated on re-runs. Thus, if pipeline has finished processing the data and you want to run it again with different set of parameters, just remove the files that need to be generated again while keeping record of the order in which they are created. 
To see what you have already done, output gets logged in TAXAassign.log file in the current folder. This is useful for keeping track of the parameters used when running the software.
If you have the OTUs abundance table (in csv or tsv format) along with OTUs sequences, the output file generated from TAXAassign can then be used withhttp://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/convIDs.pl to annotate the table.
To gain better understanding of BioSQL, refer to the [BioSQL tutorial](http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/BIOSQL_tutorial.pdf).   
In addition to finding taxonomic assignments, TAXAassign also generates taxa abundances for given set of sequences which may then be collated with taxa abundances from other sets of sequences to form a frequency table (akin to OTUs table) that can be analyzed further in an environmental context using **[TAXAenv](http://quince-srv2.eng.gla.ac.uk:8080)** for multi-variate statistical analysis.


# Dependencies
### Installing GNU Parallel

You can download GNU Parallel from http://www.gnu.org/software/parallel/

###Installing Blastn 2.28+ and NCBI's nt database

The current version of TAXAassign only works with Blastn 2.28+ as previous versions dont give Taxa IDs in output format.
You can download it from http://www.ncbi.nlm.nih.gov/books/NBK1763/. 

For nt database, go to your local Blastn installation folder and use update_blastdb.pl as
```
update_blastdb.pl --showall
update_blastdb.pl nt 
```

Once downloaded, edit **TAXAassign.sh**, and set **BLASTN_DIR** and **BLASTDB_DIR** to appropriate paths by finding the following section:

```
# = Parameters to set ============== #
LOGFILE="`pwd`/TAXAassign.log" # Where to save the log
BLASTN_DIR="/home/opt/ncbi-blast-2.2.28+/bin"; # Path where blastn is installed
BLASTDB_DIR="/home/opt/ncbi-blast-2.2.28+/db"; # Path where nt is installed
FASTA_FILE="" # This field should be empty
PARALLELIZE_FLAG=0
NUMBER_OF_CORES=10
NUMBER_OF_REFERENCE_MATCHES=10
MINIMUM_PERCENT_IDENT=97
MINIMUM_QUERY_COVERAGE=97
CONSENSUS_THRESHOLD=90
TAXONOMIC_LEVELS_THRESHOLD=""
# =/Parameters to set ============== #
```

### Installing BioSQL

Download BioSQL from http://biosql.org/DIST/biosql-1.0.1.tar.gz. Once the software is installed, setup a database and import the BioSQL schema. The following command line should create a database called bioseqdb on your computer belonging to the root user account:
```
mysqladmin -u root create bioseqdb
```

We can then tell MySQL to load the BioSQL scheme we downloaded above. Change to the scripts subdirectory from the unzipped BioSQL download, then use the command:
```
mysql -u root bioseqdb < biosqldb-mysql.sql
```
To update the NCBI taxonomy, change to the scripts subdirectory from the unzipped BioSQL download, then use the command (output is also shown):
```
./load_ncbi_taxonomy.pl --dbname bioseqdb --driver mysql --dbuser root --download true
Loading NCBI taxon database in taxdata:
 ... retrieving all taxon nodes in the database
 ... reading in taxon nodes from nodes.dmp
 ... insert / update / delete taxon nodes
 ... (committing nodes)
 ... rebuilding nested set left/right values
 ... reading in taxon names from names.dmp
 ... deleting old taxon names
 ... inserting new taxon names
 ... cleaning up
Done.
```

You can also use sqlite3 to store the database locally in the installation folder if you don't want to go for the centralized MySQL server option.
Last time I checked BioSQL didn't have any option to load database schema in sqlite3 directly or loading data with **load_ncbi_taxonomy.pl** script.
A work around is to dump your MySQL database to sqlite3 and place the database as **db.sqlite** in the database folder.
You can then edit the parameters section of **blast_concat_taxon.py** and set **use_MySQL=False**. The section is as follows

```
 # Parameters #########################################
 use_MySQL=True
 #MySQL server settings for BioSQL
 MySQL_server='localhost'
 MySQL_user='root'
 MySQL_password=''
 MySQL_database='bioseqdb'
 #sqlite3 database
 sqlite3_database=os.getcwd()+"/../database/db.sqlite"
 ######################################################
```

There are several conversion scripts out there to export data from MySQL to sqlite3 but not all of them work. The only thing that worked for me
was ruby-gem. The commands are as follows:
```
sudo gem install sequel
sudo gem install sqlite3
sudo gem install mysql
sequel mysql://root@localhost/bioseqdb -C sqlite://db.sqlite
```

Make sure that you have development version of both sqlite3 and MySQL installed.


### Downloading sqlite3 database

```
cd <TAXAassign_directory>
mkdir database
cd database
wget http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/db.sqlite.gz
gunzip sqlite.db.gz
```

### Additional packages required

Python: MySQLdb, sqlite3 (use easy_install or pip install)  
Perl: Getopt::Long (use perl -MCPAN -e shell)  

## Walk-through with the included test.fasta file

**Step 1:** To test TAXAassign, we will use a small dataset test.fasta comprising 10 unknown sequences. To do so, create a folder and copy test.fasta provided in the installation directory
```
[uzi@quince-srv2 ~/]$ mkdir check_TAXAassign; cd check_TAXAassign 
[uzi@quince-srv2 ~/check_TAXAassign]$ cp ~/<TAXAassign_directory>/data/test.fasta .
[uzi@quince-srv2 ~/check_TAXAassign]$ ls
test.fasta
```
**Step 2:** You may execute TAXAassign.sh without any arguments to look at the syntax. If you are running TAXAassign on a multicore server, you may turn the parallel processing option on by -p switch and then set the number of cores that you want to use using -c switch.
You can specify the total number of reference matches per query sequence using -r switch, percentage identity or similarity using -m switch, and query sequence coverage using -c switch.
TAXAassign takes a consensus of multiple hits, so if you specify -t 70 then it means the sequence will be classified only if more than 70% of the hits agree on the same taxonomic level assignment. To specify thresholds for each taxonomic levek, use -a switch.
```
[uzi@quince-srv2 ~/check_TAXAassign]$ bash ~/TAXAassign_v0.4/TAXAassign.sh 
Usage:
    bash TAXAassign.sh -f <fasta_file.fasta> [options]
Options:
    -p Turn parallel processing on
    -c Number of cores to use (Default: 10)
    -r Number of reference matches (Default: 10)
    -m Minimum percentage identity in blastn (Default: 97)
    -q Minimum query coverage in blastn (Default: 97)
    -a Threshold at different taxonomic levels (Default:"-m,-m,-m,-m,-m,-m" where -m is the minimum percentage identity argument)
       The order is as follows: Phylum,Class,Order,Family,Genus,Species
       For example, -a "60,70,80,95,95,97"
    -t Consensus threshold (Default: 90)
```

**Step 3:** You can then run the file test.fasta as follows:

```    
[uzi@quince-srv2 ~/check_TAXAassign]$ ~/TAXAassign_v0.4/TAXAassign.sh -p -c 10 -t 70 -m 60 -a "60,70,80,95,95,97" -f test.fasta 
[2013-07-19 10:17:23] TAXAassign v0.4. Copyright (c) 2013 Computational Microbial Genomics Group, University of Glasgow, UK
[2013-07-19 10:17:23] Using  /home/opt/ncbi-blast-2.2.28+/bin/blastn
[2013-07-19 10:17:23] Using  /home/uzi/TAXAassign_v0.4/scripts/blast_concat_taxon.py
[2013-07-19 10:17:23] Using  /home/uzi/TAXAassign_v0.4/scripts/blast_gen_assignments.pl
[2013-07-19 10:17:23] Using  parallel
[2013-07-19 10:17:23] Blast against NCBI's nt database with minimum percent ident of 60%, maximum of 10 reference sequences, and evalue of 0.0001 in blastn.
[2013-07-19 10:18:06] blastn using GNU parallel took 43 seconds for test.fasta.
[2013-07-19 10:18:06] test_B.out generated successfully!
[2013-07-19 10:18:06] Filter blastn hits with minimum query coverage of 97%.
[2013-07-19 10:18:06] test_BF.out generated successfully!
[2013-07-19 10:18:06] Annotate blastn hits with NCBI's taxonomy data.
[2013-07-19 10:18:08] test_BFT.out generated successfully!
[2013-07-19 10:18:08] Generate taxonomic assignment tables from blastn hits with consensus threshold of 70%.
[2013-07-19 10:18:08] test_ASSIGNMENTS.csv, test_PHYLUM.csv, test_CLASS.csv, test_ORDER.csv, test_FAMILY.csv, test_GENUS.csv, and test_SPECIES.csv generated successfully!
[2013-07-19 10:18:08] Sequences assigned at phylum level: 10/10 (100.00000%)
[2013-07-19 10:18:08] Sequences assigned at class level: 10/10 (100.00000%)
[2013-07-19 10:18:08] Sequences assigned at order level: 10/10 (100.00000%)
[2013-07-19 10:18:08] Sequences assigned at family level: 10/10 (100.00000%)
[2013-07-19 10:18:08] Sequences assigned at genus level: 10/10 (100.00000%)
[2013-07-19 10:18:08] Sequences assigned at species level: 2/10 (20.00000%)
```
**Step 4:** Check the contents of the assignment file. The columns are arranged in the following order: Sequence ID, Phylum, Class, Order, Family, Genus, Species.
```
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_ASSIGNMENTS.csv
seq6,Proteobacteria,Gammaproteobacteria,Pseudomonadales,Pseudomonadaceae,Pseudomonas,__Unclassified__
seq3,Proteobacteria,Gammaproteobacteria,Pseudomonadales,Pseudomonadaceae,Pseudomonas,__Unclassified__
seq7,Proteobacteria,Gammaproteobacteria,Pseudomonadales,Pseudomonadaceae,Pseudomonas,__Unclassified__
seq9,Proteobacteria,Betaproteobacteria,Burkholderiales,Burkholderiaceae,Ralstonia,Ralstonia solanacearum
seq2,Actinobacteria,Actinobacteria,Actinomycetales,Nocardiaceae,Rhodococcus,__Unclassified__
seq10,Actinobacteria,Actinobacteria,Actinomycetales,Microbacteriaceae,Microbacterium,Microbacterium laevaniformans
seq8,Proteobacteria,Alphaproteobacteria,Caulobacterales,Caulobacteraceae,Brevundimonas,__Unclassified__
seq1,Proteobacteria,Alphaproteobacteria,Caulobacterales,Caulobacteraceae,Brevundimonas,__Unclassified__
seq4,Proteobacteria,Deltaproteobacteria,Desulfovibrionales,Desulfovibrionaceae,Desulfovibrio,__Unclassified__
seq5,Bacteroidetes/Chlorobi group,Sphingobacteriia,Sphingobacteriales,Sphingobacteriaceae,Pedobacter,__Unclassified__
```
**Step 5:** Also check the abundances of taxa at different taxonomic levels: 
```
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_PHYLUM.csv
Actinobacteria,2
Bacteroidetes/Chlorobi group,1
Proteobacteria,7
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_CLASS.csv
Actinobacteria,2
Alphaproteobacteria,2
Betaproteobacteria,1
Deltaproteobacteria,1
Gammaproteobacteria,3
Sphingobacteriia,1
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_ORDER.csv
Actinomycetales,2
Burkholderiales,1
Caulobacterales,2
Desulfovibrionales,1
Pseudomonadales,3
Sphingobacteriales,1
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_FAMILY.csv
Burkholderiaceae,1
Caulobacteraceae,2
Desulfovibrionaceae,1
Microbacteriaceae,1
Nocardiaceae,1
Pseudomonadaceae,3
Sphingobacteriaceae,1
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_GENUS.csv
Brevundimonas,2
Desulfovibrio,1
Microbacterium,1
Pedobacter,1
Pseudomonas,3
Ralstonia,1
Rhodococcus,1
[uzi@quince-srv2 ~/check_TAXAassign]$ cat test_SPECIES.csv
Microbacterium laevaniformans,1
Ralstonia solanacearum,1
__Unclassified__,8
[uzi@quince-srv2 ~/check_TAXAassign]$ 

```
