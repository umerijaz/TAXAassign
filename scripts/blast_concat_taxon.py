#!/usr/bin/python
# ***************************************************************
# Name:      	blast_concat_taxon.py
# Purpose:   	This scripts takes a blast output file, extracts the
#         	gid and appends taxonomic path
#          
#
#            	This script can filter blast files generated through the
#            	the following command:
#
#            	blastn -db <path_to_nt> -query <fasta_file>
#            	-out <blast_file> -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids" <other options>
#
# Dependencies: BioSQL
# 		Download BioSQL from http://biosql.org/DIST/biosql-1.0.1.tar.gz. Once the software is installed, 
#		setup a database and import the BioSQL schema. The following command line should create a new database
#		on your own computer called bioseqdb, belonging to the root user account:
#			mysqladmin -u root create bioseqdb
#		We can then tell MySQL to load the BioSQL scheme we downloaded above. Change to the scripts subdirectory 
#		from the unzipped BioSQL download, then use the command:
#			mysql -u root bioseqdb < biosqldb-mysql.sql
#		To update the NCBI taxonomy, change to the scripts subdirectory from the unzipped BioSQL download, then use
#		the command (output is also shown):
#			 ./load_ncbi_taxonomy.pl --dbname bioseqdb --driver mysql --dbuser root --download true
#				Loading NCBI taxon database in taxdata:
#        					... retrieving all taxon nodes in the database
#        					... reading in taxon nodes from nodes.dmp
#        					... insert / update / delete taxon nodes
#        					... (committing nodes)
#        					... rebuilding nested set left/right values
#        					... reading in taxon names from names.dmp
#        					... deleting old taxon names
#        					... inserting new taxon names
#        					... cleaning up
#				Done.
#
#                       You can also use sqlite3 to store the database in case you don't want to go for MySQL server option.
#                       Last time I checked BioSQL didn't have any option to load database schema  in sqlite directly or loading data with load_ncbi_taxonomy.pl script
#                       A work around is to dump your MySQL database to sqlite3 and place the database as db.sqlite in the database folder.
#                       You can then edit the parameters section below and set use_MySQL=False. The section is as follows
#
#                                       # Parameters #########################################
#                                       use_MySQL=True 
#                                       #MySQL server settings for BioSQL
#                                       MySQL_server='localhost'
#                                       MySQL_user='root'
#                                       MySQL_password=''
#                                       MySQL_database='bioseqdb'
#                                       #sqlite3 database
#                                       sqlite3_database=os.getcwd()+"/../database/db.sqlite"
#                                       #####################################################
#
#                       There are several conversion scripts out there to export data from MySQL to sqlite3  but not all of them work. The only thing that worked for me
#                       was ruby-gem. The commands are as follows:
#
#                       sudo gem install sequel
#                       sudo gem install sqlite3
#                       sudo gem install mysql
#                       sequel mysql://root@localhost/bioseqdb -C sqlite://db.sqlite
#
#                       Make sure that you have development version of both sqlite3 and MySQL installed.
#
# Version:   0.3
# History:   0.1 to 0.2 Full paths instead of taxonomic levels
#	     0.2 to 0.3 sqlite3 supported as well	
#		
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Created:   2013-06-25
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

import re,urllib,sys,time,getopt
from xml.dom import minidom
import MySQLdb as mdb
import os
import sqlite3 as lite

def print_record(record,cur):
    ncbi_taxon_id=[]
    ncbi_taxon_ids=record.split("\t")[13]
    
    for ncbi_taxon_id in ncbi_taxon_ids.split(";"):
    	taxonomy=[]
    
    	#We only have the ncbi_taxon_id so use it to get the first bioSQL taxon_id
    	#then iterate up through the taxonomy
    	status=cur.execute("""SELECT taxon_name.name, taxon.node_rank, taxon.parent_taxon_id FROM taxon, taxon_name WHERE taxon.taxon_id=taxon_name.taxon_id AND taxon_name.name_class='scientific name' AND taxon.ncbi_taxon_id = %s""" % (ncbi_taxon_id,))
    	if status:
    		name, rank, parent_taxon_id = cur.fetchone()
		taxonomy.append(name+":"+rank)
    		taxon_id = parent_taxon_id
    		while taxon_id:
			cur.execute("""SELECT taxon_name.name, taxon.node_rank, taxon.parent_taxon_id FROM taxon, taxon_name WHERE taxon.taxon_id=taxon_name.taxon_id AND taxon_name.name_class='scientific name' AND taxon.taxon_id = %s""" % (taxon_id,))
			name, rank, parent_taxon_id = cur.fetchone()
			if taxon_id == parent_taxon_id:
				break
        		taxonomy.insert(0,name+":"+rank)
			taxon_id=parent_taxon_id
    	print record.rstrip('\n')+"\t"+";".join(taxonomy)
def usage():
    print 'Usage:'
    print '\tpython blast_concat_taxon.py -b <blast_file> > <output_file>'
def main(argv):
   
    # Parameters #########################################
    use_MySQL=True	
    #MySQL server settings for BioSQL
    MySQL_server='localhost'
    MySQL_user='root'
    MySQL_password=''
    MySQL_database='bioseqdb'
    #sqlite3 database
    sqlite3_database=os.getcwd()+"/../database/db.sqlite"
    #####################################################


    #either a blast file or a single blast record
    blast_file=''

    try:
        opts, args = getopt.getopt(argv,"hb:",["blast_file="])
    except getopt.GetoptError:
        usage()
        exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            exit()
        elif opt in ("-b", "--blast_file"):
            blast_file = arg
    if (blast_file==""):
        usage()
        exit()         
   
    #required for MySQL    
    con = None
    cur = None	

    if use_MySQL:
    	#connect to MySQL server
    	try:
		con = mdb.connect(MySQL_server,MySQL_user,MySQL_password,MySQL_database); 
   		cur = con.cursor()
    	except mdb.Error, e:
    		print "Error %d: %s" % (e.args[0],e.args[1])
    		sys.exit(1)
    else:
	con=lite.connect(sqlite3_database)
	cur=con.cursor()

	
    #check if it is a single record or a file
    if re.findall(r'gi\|(\w.*?)\|',blast_file):
        print_record(blast_file,cur)
    else:
        ins=open(blast_file,"r")    
	for line in ins:
            print_record(line,cur)
        ins.close()
    if con:    
    	con.close();


if __name__ == '__main__':
    main(sys.argv[1:])
