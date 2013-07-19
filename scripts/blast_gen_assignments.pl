#!/usr/bin/perl
## ***************************************************************
# Name:      blast_gen_assignments.pl
# Purpose:   This program takes the blast file with taxonomic path appended as last column
# 	     and gives consensus assignment based on a cutoff value while making sure different taxonomic assignments are above a minimum threshold
# 	     
# 	     This script depends on the blast output file generated from blast_concat_taxon.py script with order as follows:
# 	      fasta_file
# 	        |
# 	        |
# 	        V
# 	      blastn -db <path_to_nt> -query <fasta_file>
# 	      -out <blast_file> -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids" <other options>
#		|
#		|
#		V
#	      blast_concat_taxon.py
#	      	|
#	      	|
#	      	V
#	      blast_gen_assignments.pl	
# Version:   0.1
# Known bug: Need to fix this: The key,value pair in %ASSIGNMENTS_hash generates a "HASH(.." key for which temporary hack is to pipe the output through grep -v "HASH("
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

use strict;
use Getopt::Long;

my %opts; #store the input arguments
GetOptions(\%opts,
        'blast_file|b=s',
	'consensus_cutoff|c=i',
	'taxonomic_levels_threshold|a=s',
);

if(not defined $opts{"blast_file"})
        {
print <<EOF;
Usage:
        perl blast_gen_assignments.pl -b <blast_file> -i <consensus_cutoff> -a <taxonomic_levels_threshold>
		Default for -a : "97,97,97,97,97,97" (order is as follows: Phylum,Class,Order,Family,Genus,Species)
		Default for -c : 90
EOF
        exit;
        }


my $blast_file = $opts{"blast_file"};
unless(-e $blast_file)
        {
        print "Error:",$blast_file, " not found!\n";
        exit;
        }
my $consensus_cutoff=90;
unless(not defined $opts{"consensus_cutoff"})
        {$consensus_cutoff=$opts{"consensus_cutoff"};}

my $taxonomic_levels_threshold="97,97,97,97,97,97";
unless(not defined $opts{"taxonomic_levels_threshold"})
        {$taxonomic_levels_threshold=$opts{"taxonomic_levels_threshold"};}


# Store assignments in a hash
my %ASSIGNMENTS_hash={};
my $line;
my @tokens;
my $id;
my $path;

# Specify which columns (numbering starting from zero) contain the required data
my $id_column=0;
my $ident_column=2;
my $path_column=14;

@tokens=split(/,/,$taxonomic_levels_threshold);
my $phylum_threshold=$tokens[0];
my $class_threshold=$tokens[1];
my $order_threshold=$tokens[2];
my $family_threshold=$tokens[3];
my $genus_threshold=$tokens[4];
my $species_threshold=$tokens[5];


open(FILE, $blast_file) or die "Can't open $blast_file\n";
while(my $line=<FILE>){
	chomp($line);
	@tokens=split(/\t/,$line);
	$id=$tokens[$id_column];
	$path=$tokens[$path_column].";".$tokens[$ident_column].":perc_ident";
        if (not defined @{$ASSIGNMENTS_hash{$id}})
        	{
                push @{$ASSIGNMENTS_hash{$id}}, $path;
                }
        else    {
                unless($path ~~ @{$ASSIGNMENTS_hash{$id}}) #enter unique path
                      {
                      push @{$ASSIGNMENTS_hash{$id}}, $path;
                      }
                }
}
my $del=",";
my $debug_mode=0;

# print recovered tree ####################################################
my($key, $value);
while ( ($key, $value) = each(%ASSIGNMENTS_hash) ) {
    print "$key".$del;
       my @phylum_assignments=();my $phylum_unknowns=0;
       my @class_assignments=();my $class_unknowns=0;
       my @order_assignments=();my $order_unknowns=0;
       my @family_assignments=();my $family_unknowns=0;
       my @genus_assignments=();my $genus_unknowns=0;
       my @species_assignments=();my $species_unknowns=0;
       print "\n" if $debug_mode;	
       foreach (@{$value}){
               print $del if $debug_mode;
	       my @assignment;
	       my @percentage_identity=split(/:/,join(/,/,grep {/perc_ident/} split(/;/,$_)));
	       @assignment=split(/:/,join(/,/,grep {/phylum/} split(/;/,$_)));
	       if (scalar(@assignment)>0) {
			print $assignment[0].$del if $debug_mode;

			if (($assignment[0]!~/^(uncultured|unclassified|unidentified)/) and ($percentage_identity[0]>=$phylum_threshold))
				{
				push @phylum_assignments, $assignment[0];
				}
			else	{
				$phylum_unknowns++;
				}
			} 
               else {
                        print "__Unknown__\t" if $debug_mode;
			$phylum_unknowns++;
		    }
	       @assignment=split(/:/,join(/,/,grep {/class/} split(/;/,$_)));
	       if (scalar(@assignment)>0) {
			print $assignment[0].$del if $debug_mode;
                        if (($assignment[0]!~/^(uncultured|unclassified|unidentified)/) and ($percentage_identity[0]>=$class_threshold))
                                {
				push @class_assignments, $assignment[0];
				}
                        else    {
                                $class_unknowns++;
                                }		
			} 
	       else {
			print "__Unknown__\t" if $debug_mode;
			$class_unknowns++;
		    }
	       @assignment=split(/:/,join(/,/,grep {/order/} split(/;/,$_)));
	       if (scalar(@assignment)>0) {
			print $assignment[0].$del if $debug_mode;
                        if (($assignment[0]!~/^(uncultured|unclassified|unidentified)/) and ($percentage_identity[0]>=$order_threshold))
                                {
                                push @order_assignments, $assignment[0];
                                }       
                        else    {
                                $order_unknowns++;
                                }
			} 
	       else {
			print "__Unknown__\t" if $debug_mode;
			$order_unknowns++;
                    }
               @assignment=split(/:/,join(/,/,grep {/family/} split(/;/,$_)));
	       if (scalar(@assignment)>0) {
			print $assignment[0].$del if $debug_mode;
                        if (($assignment[0]!~/^(uncultured|unclassified|unidentified)/) and ($percentage_identity[0]>=$family_threshold))
                                {
                                push @family_assignments, $assignment[0];
                                }       
                        else    {
                                $family_unknowns++;
                                }
			} 
	       else {
			print "__Unknown__\t" if $debug_mode;
			$family_unknowns++;
               }
               @assignment=split(/:/,join(/,/,grep {/genus/} split(/;/,$_)));
	       if (scalar(@assignment)>0) {
			print $assignment[0].$del if $debug_mode;
                        if (($assignment[0]!~/^(uncultured|unclassified|unidentified)/) and ($percentage_identity[0]>=$genus_threshold))
                                {
                                push @genus_assignments, $assignment[0];
                                }       
                        else    {
                                $genus_unknowns++;
                                }			
			} 
	       else {
			print "__Unknown__\t" if $debug_mode;
			$genus_unknowns++;
		    }
               @assignment=split(/:/,join(/,/,grep {/species/} split(/;/,$_)));
	       if (scalar(@assignment)>0) {
			print $assignment[0]."\n" if $debug_mode;
                        if (($assignment[0]!~/^(uncultured|unclassified|unidentified)/) and ($percentage_identity[0]>=$species_threshold))
                                {
                                push @species_assignments, $assignment[0];
                                }       
                        else    {
                                $species_unknowns++;
                                }			
			} 
	       else {
			print "__Unknown__\n" if $debug_mode;
			$species_unknowns++;
	       }
  
       }
	
    my %arr_counts;
    my @qualifying_terms;
    print "Consensus:\t" if $debug_mode;

    for (@phylum_assignments) { $arr_counts{$_}++ };
    my @qualifying_terms = grep { $arr_counts{$_} > (($consensus_cutoff/100)*(scalar(@{$value})-$phylum_unknowns)) } keys %arr_counts;
    if (scalar(@qualifying_terms)>0) {print $qualifying_terms[0].$del;} else {print "__Unclassified__".$del;};

    %arr_counts={};
    @qualifying_terms=();
    for (@class_assignments) { $arr_counts{$_}++ };
    my @qualifying_terms = grep { $arr_counts{$_} > (($consensus_cutoff/100)*(scalar(@{$value})-$class_unknowns)) } keys %arr_counts;
    if (scalar(@qualifying_terms)>0) {print $qualifying_terms[0].$del;} else {print "__Unclassified__".$del;};


    %arr_counts={};
    @qualifying_terms=();
    for (@order_assignments) { $arr_counts{$_}++ };
    my @qualifying_terms = grep { $arr_counts{$_} > (($consensus_cutoff/100)*(scalar(@{$value})-$order_unknowns)) } keys %arr_counts;
    if (scalar(@qualifying_terms)>0) {print $qualifying_terms[0].$del;} else {print "__Unclassified__".$del;};


    %arr_counts={};
    @qualifying_terms=();
    for (@family_assignments) { $arr_counts{$_}++ };
    my @qualifying_terms = grep { $arr_counts{$_} > (($consensus_cutoff/100)*(scalar(@{$value})-$family_unknowns)) } keys %arr_counts;
    if (scalar(@qualifying_terms)>0) {print $qualifying_terms[0].$del;} else {print "__Unclassified__".$del;};

    %arr_counts={};
    @qualifying_terms=();
    for (@genus_assignments) { $arr_counts{$_}++ };
    my @qualifying_terms = grep { $arr_counts{$_} > (($consensus_cutoff/100)*(scalar(@{$value})-$genus_unknowns)) } keys %arr_counts;
    if (scalar(@qualifying_terms)>0) {print $qualifying_terms[0].$del;} else {print "__Unclassified__".$del;};

    %arr_counts={};
    @qualifying_terms=();
    for (@species_assignments) { $arr_counts{$_}++ };
    my @qualifying_terms = grep { $arr_counts{$_} > (($consensus_cutoff/100)*(scalar(@{$value})-$species_unknowns)) } keys %arr_counts;
    if (scalar(@qualifying_terms)>0) {print $qualifying_terms[0];} else {print "__Unclassified__";};

    print "\n";
}

close(FILE);


