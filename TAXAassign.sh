#!/bin/bash
# ***************************************************************
# Name:    	TAXAassign.sh
# Purpose: 	Script to annotate sequences at different taxonomic levels using  NCBI's taxonomy
#	   	
# Dependencies:  
# 		GNU Parallel 
#    			Install: http://www.gnu.org/software/parallel/
# 		Blastn 2.28+ (The previous versions dont give Taxa IDs)
#    			Software: http://www.ncbi.nlm.nih.gov/books/NBK1763/
#    				For nt database, go to your local blastn db folder and use update_blastdb.pl: 
#        			update_blastdb.pl --showall
#        			update_blastdb.pl nt 
#               	Blastn output format used:
#				blastn -db <path_to_nt> -query <fasta_file>
#               		-out <blast_file> -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids" <other options>
#		BioSQL
#               	Download BioSQL from http://biosql.org/DIST/biosql-1.0.1.tar.gz. Once the software is installed, 
#               	setup a database and import the BioSQL schema. The following command line should create a new database
#               	on your own computer called bioseqdb, belonging to the root user account:
#                       	mysqladmin -u root create bioseqdb
#               	We can then tell MySQL to load the BioSQL scheme we downloaded above. Change to the scripts subdirectory 
#               	from the unzipped BioSQL download, then use the command:
#                       	mysql -u root bioseqdb < biosqldb-mysql.sql
#               	To update the NCBI taxonomy, change to the scripts subdirectory from the unzipped BioSQL download, then use
#               	the command (output is also shown):
#                        	./load_ncbi_taxonomy.pl --dbname bioseqdb --driver mysql --dbuser root --download true
#                               	Loading NCBI taxon database in taxdata:
#                                       	        ... retrieving all taxon nodes in the database
#                                              	 	... reading in taxon nodes from nodes.dmp
#                                               	... insert / update / delete taxon nodes
#                                               	... (committing nodes)
#                                               	... rebuilding nested set left/right values
#                                               	... reading in taxon names from names.dmp
#                                               	... deleting old taxon names
#                                               	... inserting new taxon names
#                                               	... cleaning up
#                               	Done.
#
#			You can also use sqlite3 to store the database in case you don't want to go for MySQL server option.
#			Last time I checked BioSQL didn't have any option to load database schema  in sqlite directly or loading data with load_ncbi_taxonomy.pl script
#			A work around is to dump your MySQL database to sqlite3 and place the database as db.sqlite in the database folder.
#			You can then edit the parameters section in blast_concat_taxon.py and set use_MySQL=False. The section is as follows
#
#				    	# Parameters #########################################
#    					use_MySQL=True 
#    					#MySQL server settings for BioSQL
#    					MySQL_server='localhost'
#   		 			MySQL_user='root'
#    					MySQL_password=''
#    					MySQL_database='bioseqdb'
#    					#sqlite3 database
#    					sqlite3_database=os.getcwd()+"/../database/db.sqlite"
#    					#####################################################
#
#			There are several conversion scripts out there to export data from MySQL to sqlite3  but not all of them work. The only thing that worked for me
#			was ruby-gem. The commands are as follows:
#
#			sudo gem install sequel
#			sudo gem install sqlite3
#			sudo gem install mysql
#			sequel mysql://root@localhost/bioseqdb -C sqlite://db.sqlite
#
#			Make sure that you have development version of both sqlite3 and MySQL installed.
#
# Version:   0.4
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Last modified:   2013-07-18
# License:   Copyright (c) 2013 Computational Microbial Genomics Group, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/     

HELPDOC=$( cat <<EOF
Script to annotate sequences at different taxonomic levels using  NCBI's taxonomy

Usage:
    bash `basename $0` -f <fasta_file.fasta> [options]
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

EOF
)

set -o errexit

# = Parameters to set ============== #
LOGFILE="`pwd`/TAXAassign.log" # Where to save the log
BLASTN_DIR="/home/opt/ncbi-blast-2.2.28+/bin"; # Path where blastn is installed
BLASTDB_DIR="/home/opt/ncbi-blast-2.2.28+/db"; # Path where nt is installed
FASTA_FILE=""   # This field should be empty
PARALLELIZE_FLAG=0
NUMBER_OF_CORES=10
NUMBER_OF_REFERENCE_MATCHES=10
MINIMUM_PERCENT_IDENT=97
MINIMUM_QUERY_COVERAGE=97
CONSENSUS_THRESHOLD=90
TAXONOMIC_LEVELS_THRESHOLD=""
# =/Parameters to set ============== #

CURRENT_DIR=`pwd`

# = Enable FP support ============== #
# By default, there is limited capability in bash to handle floating point
# operations. In this script bc is used to calculate the floating point operations.
# $float_scale parameter specifies the precision of the floating point.
# Reference: http://www.linuxjournal.com/content/floating-point-math-bash

float_scale=5
# Evaluate a floating point number expression.

function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}


# Evaluate a floating point number conditional expression.

function float_cond()
{
    local cond=0
    if [[ $# -gt 0 ]]; then
        cond=$(echo "$*" | bc -q 2>/dev/null)
        if [[ -z "$cond" ]]; then cond=0; fi
        if [[ "$cond" != 0  &&  "$cond" != 1 ]]; then cond=0; fi
    fi
    local stat=$((cond == 0))
    return $stat
}

function ceil () {
  echo "define ceil (x) {if (x<0) {return x/1} \
        else {if (scale(x)==0) {return x} \
        else {return x/1 + 1 }}} ; ceil($1)" | bc;
 }

# =/Enable FP support ============== #


# Create directories if they don't exist yet
function create_dirs() {
    local dir
    for dir in "$@"
    do
        if [ ! -d "$dir" ]; then
            mkdir "$dir"
        fi
    done
}

# Check if files exist
function check_prog() {
    local prog
    for prog in "$@"
    do
    if which $prog >/dev/null; then
        TAXAassign_print 'Using ' $prog
    else
        echo "$prog not in your path" >&2; exit 1;
    fi

    done
}

function skip_gen_file(){
    if [ -f "$1" ]; then
    echo "true"
    else
    echo "false"
    fi
}

function skip_gen_dir(){
    if [ -d "$1" ]; then
        echo "true"
    else
        echo "false"
    fi
}

function TAXAassign_print() {
    echo [`date "+%Y-%m-%d %H:%M:%S"`] "$@" | tee -a $LOGFILE
}


# Parse options
while getopts ":phc:r:m:f:t:q:a:" opt; do
    case $opt in
        p)
            PARALLELIZE_FLAG=1
            ;;
        f)
            FASTA_FILE=$OPTARG
            ;;
	m)
	    MINIMUM_PERCENT_IDENT=$OPTARG
	    ;;	
        c)
            NUMBER_OF_CORES=$OPTARG
            ;;
        r)
            NUMBER_OF_REFERENCE_MATCHES=$OPTARG
            ;;
	t)
	    CONSENSUS_THRESHOLD=$OPTARG
	    ;;
	q)
	    MINIMUM_QUERY_COVERAGE=$OPTARG
	    ;;
	a)
	    TAXONOMIC_LEVELS_THRESHOLD=$OPTARG
	    ;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        \?)
            echo "$HELPDOC"
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done
if [ -z $FASTA_FILE ] 
then
        echo "$HELPDOC"
        exit 1
fi

if [ -z "$TAXONOMIC_LEVELS_THRESHOLD" ]; then
        TAXONOMIC_LEVEL_THRESHOLD="$MINIMUM_PERCENT_IDENT,$MINIMUM_PERCENT_IDENT,$MINIMUM_PERCENT_IDENT,$MINIMUM_PERCENT_IDENT,$MINIMUM_PERCENT_IDENT,$MINIMUM_PERCENT_IDENT"
fi

OIFS=$IFS;
IFS=",";
TLTArray=($TAXONOMIC_LEVELS_THRESHOLD);
IFS=$OIFS;

if [ "${#TLTArray[@]}" != "6" ]; then
        echo "$HELPDOC"
        exit 1
fi

for i in "${TLTArray[@]}"
do
   :
   if ! [[ $i =~ ^-?[0-9]+$ ]]; then
        echo "$HELPDOC"
        exit 1
   fi
done


if ! [[ $MINIMUM_PERCENT_IDENT =~ ^-?[0-9]+$ ]] || ! [[ $NUMBER_OF_CORES =~ ^-?[0-9]+$ ]] || ! [[ $NUMBER_OF_REFERENCE_MATCHES =~ ^-?[0-9]+$ ]] || ! [[ $CONSENSUS_THRESHOLD =~ ^-?[0-9]+$ ]] || ! [[ $MINIMUM_QUERY_COVERAGE =~ ^-?[0-9]+$ ]]; then
	echo "$HELPDOC"
	exit 1
fi

# Using /usr/bin/dirname to get the full path to this script location
# without being affected from where this script was invoked
export TAXAASSIGN_DIR=$(cd "$(dirname "$0")"; pwd)
TAXAassign_print "TAXAassign v0.4. Copyright (c) 2013 Computational Microbial Genomics Group, University of Glasgow, UK"
 
check_prog $BLASTN_DIR/blastn
check_prog $TAXAASSIGN_DIR/scripts/blast_concat_taxon.py
check_prog $TAXAASSIGN_DIR/scripts/blast_gen_assignments.pl
if [ $PARALLELIZE_FLAG -eq 1 ]
then
    check_prog parallel
fi


fileName=`echo "$(basename $FASTA_FILE)" | cut -d'.' -f1`

# Run blastn
TAXAassign_print "Blast against NCBI's nt database with minimum percent ident of $MINIMUM_PERCENT_IDENT%, maximum of $NUMBER_OF_REFERENCE_MATCHES reference sequences, and evalue of 0.0001 in blastn."
# Format for blastn
blastOutFmt="\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids\""
blastFileName=$fileName'_B'
if [ "$(skip_gen_file $blastFileName'.out')" == "true" ];then
   TAXAassign_print $blastFileName'.out' already exists. Skipping this step.
elif [ $PARALLELIZE_FLAG -eq 1 ]; then

    # Get the file size in KB
    sizeFileBytes=$(du -b ${FASTA_FILE} | sed 's/\([0-9]*\)\(.*\)/\1/')
    sizeChunks=$(ceil $(float_eval "$sizeFileBytes / ($NUMBER_OF_CORES * 1024)"))
    sizeChunksString="${sizeChunks}k"
    startTime=`date +%s`

    cat $FASTA_FILE | parallel --block $sizeChunksString --recstart '>' --pipe $BLASTN_DIR/blastn -perc_identity $MINIMUM_PERCENT_IDENT -evalue 0.00001 -dust no -num_threads 1 -outfmt $blastOutFmt -max_target_seqs $NUMBER_OF_REFERENCE_MATCHES -db $BLASTDB_DIR'/nt' -query - > $blastFileName'.out'

    TAXAassign_print "blastn using GNU parallel took $(expr `date +%s` - $startTime) seconds for $FASTA_FILE".
    TAXAassign_print $blastFileName'.out' generated successfully!
else
    startTime=`date +%s`
    $BLASTN_DIR/blastn -db  $BLASTDB_DIR'/nt' -query $FASTA_FILE -perc_identity $MINIMUM_PERCENT_IDENT -out $blastFileName'.out' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids" -max_target_seqs $NUMBER_OF_REFERENCE_MATCHES -evalue 0.00001 -dust no -num_threads 1
    TAXAassign_print "blastn took $(expr `date +%s` - $startTime) seconds for $FASTA_FILE".
    TAXAassign_print $blastFileName'.out' generated successfully!
fi


TAXAassign_print "Filter blastn hits with minimum query coverage of $MINIMUM_QUERY_COVERAGE%."
blastFilteredFileName=$fileName'_BF'
if [ "$(skip_gen_file $blastFilteredFileName'.out')" == "true" ];then
   TAXAassign_print $blastFilteredFileName'.out' already exists. Skipping this step.
else 
   cat $blastFileName'.out' | awk -F"\t" -v pattern=$MINIMUM_QUERY_COVERAGE '$13>pattern{print $0}' > $blastFilteredFileName'.out'	
   TAXAassign_print $blastFilteredFileName'.out' generated successfully!	
fi


blastFileNameWithTaxonomy=$fileName'_BFT'

TAXAassign_print "Annotate blastn hits with NCBI's taxonomy data."

if [ "$(skip_gen_file $blastFileNameWithTaxonomy'.out')" == "true" ];then
   TAXAassign_print $blastFileNameWithTaxonomy'.out' already exists. Skipping this step.
elif [ $PARALLELIZE_FLAG -eq 1 ]; then
   cat $blastFilteredFileName'.out' | parallel -j $NUMBER_OF_CORES python $TAXAASSIGN_DIR/scripts/blast_concat_taxon.py -b {} > $blastFileNameWithTaxonomy'.out'
   TAXAassign_print $blastFileNameWithTaxonomy'.out' generated successfully!
else
   python $TAXAASSIGN_DIR/scripts/blast_concat_taxon.py -b $blastFilteredFileName'.out' > $blastFileNameWithTaxonomy'.out'
   TAXAassign_print $blastFileNameWithTaxonomy'.out' generated successfully!	
fi


TAXAassign_print "Generate taxonomic assignment tables from blastn hits with consensus threshold of $CONSENSUS_THRESHOLD%."
if [[ "$(skip_gen_file $fileName'_ASSIGNMENTS.csv')" == "true" || "$(skip_gen_file $fileName'_PHYLUM.csv')" == "true" || "$(skip_gen_file $fileName'_CLASS.csv')" == "true" || "$(skip_gen_file $fileName'_ORDER.csv')" == "true" || "$(skip_gen_file $fileName'_FAMILY.csv')" == "true" || "$(skip_gen_file $fileName'_GENUS.csv')" == "true" || "$(skip_gen_file $fileName'_SPECIES.csv')" == "true" ]];then
	TAXAassign_print Assignment files  already exists. Skipping this step.
else
        perl $TAXAASSIGN_DIR/scripts/blast_gen_assignments.pl -b $blastFileNameWithTaxonomy'.out' -c $CONSENSUS_THRESHOLD -a "$TAXONOMIC_LEVELS_THRESHOLD" | grep -v "^HASH(" >  $fileName'_ASSIGNMENTS.csv'
	totalReads=$(grep -c ">" $FASTA_FILE)
	phylumLevelAssignments=$(cut -d, -f2 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | wc -l)
	classLevelAssignments=$(cut -d, -f3 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | wc -l)
	orderLevelAssignments=$(cut -d, -f4 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | wc -l)
	familyLevelAssignments=$(cut -d, -f5 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | wc -l)
	genusLevelAssignments=$(cut -d, -f6 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | wc -l)
	speciesLevelAssignments=$(cut -d, -f7 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | wc -l)

	cut -d, -f2 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | sed 's/ /_____/g'| sort | uniq -c | awk '{gsub("_____"," ",$2);print $2","$1}' > $fileName'_PHYLUM.csv'
	if [ "$totalReads" -ne "$phylumLevelAssignments" ]; then 
		echo "__Unclassified__,$(float_eval "$totalReads - $phylumLevelAssignments ")" >> $fileName'_PHYLUM.csv'
	fi
	cut -d, -f3 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | sed 's/ /_____/g'| sort | uniq -c | awk '{gsub("_____"," ",$2);print $2","$1}' > $fileName'_CLASS.csv'
	if [ "$totalReads" -ne "$classLevelAssignments" ]; then
		echo "__Unclassified__,$(float_eval "$totalReads - $classLevelAssignments ")" >> $fileName'_CLASS.csv'
	fi
	cut -d, -f4 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | sed 's/ /_____/g'| sort | uniq -c | awk '{gsub("_____"," ",$2);print $2","$1}' > $fileName'_ORDER.csv'
	if [ "$totalReads" -ne "$orderLevelAssignments" ]; then
		echo "__Unclassified__,$(float_eval "$totalReads - $orderLevelAssignments ")" >> $fileName'_ORDER.csv'
	fi
	cut -d, -f5 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | sed 's/ /_____/g'| sort | uniq -c | awk '{gsub("_____"," ",$2);print $2","$1}' > $fileName'_FAMILY.csv'
	if [ "$totalReads" -ne "$familyLevelAssignments" ]; then
		echo "__Unclassified__,$(float_eval "$totalReads - $familyLevelAssignments ")" >> $fileName'_FAMILY.csv'
	fi
	cut -d, -f6 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | sed 's/ /_____/g'| sort | uniq -c | awk '{gsub("_____"," ",$2);print $2","$1}' > $fileName'_GENUS.csv'
	if [ "$totalReads" -ne "$genusLevelAssignments" ]; then
		echo "__Unclassified__,$(float_eval "$totalReads - $genusLevelAssignments ")" >> $fileName'_GENUS.csv'
	fi
	cut -d, -f7 $fileName'_ASSIGNMENTS.csv' | grep -v "__Unclassified__" | sed 's/ /_____/g'| sort | uniq -c | awk '{gsub("_____"," ",$2);print $2","$1}' > $fileName'_SPECIES.csv'
	if [ "$totalReads" -ne "$speciesLevelAssignments" ]; then
		echo "__Unclassified__,$(float_eval "$totalReads - $speciesLevelAssignments ")" >> $fileName'_SPECIES.csv'
	fi

	TAXAassign_print $fileName'_ASSIGNMENTS.csv', $fileName'_PHYLUM.csv', $fileName'_CLASS.csv', $fileName'_ORDER.csv', $fileName'_FAMILY.csv', $fileName'_GENUS.csv', and $fileName'_SPECIES.csv' generated successfully!
	TAXAassign_print "Sequences assigned at phylum level: $phylumLevelAssignments/$totalReads ($(float_eval "($phylumLevelAssignments / $totalReads) * 100")%)"
        TAXAassign_print "Sequences assigned at class level: $classLevelAssignments/$totalReads ($(float_eval "($classLevelAssignments / $totalReads) * 100")%)"
        TAXAassign_print "Sequences assigned at order level: $orderLevelAssignments/$totalReads ($(float_eval "($orderLevelAssignments / $totalReads) * 100")%)"
        TAXAassign_print "Sequences assigned at family level: $familyLevelAssignments/$totalReads ($(float_eval "($familyLevelAssignments / $totalReads) * 100")%)"
        TAXAassign_print "Sequences assigned at genus level: $genusLevelAssignments/$totalReads ($(float_eval "($genusLevelAssignments / $totalReads) * 100")%)"
        TAXAassign_print "Sequences assigned at species level: $speciesLevelAssignments/$totalReads ($(float_eval "($speciesLevelAssignments / $totalReads) * 100")%)"
fi

