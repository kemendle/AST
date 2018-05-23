#! /usr/bin/perl -w

###############################################################################
#Copyright 2013 Chan Zhou
###############################################################################
#This program is distributed under the terms of the GNU General Public License.
###############################################################################
#   AST.pl is a free free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

#Program: Automated Sampling representative sequences over Taxa (AST)
#website: http://csbl.bmb.uga.edu/~zhouchan/AST.php
#Author: Chan Zhou
#Contact: zhouchan99@gmail.com
#
#Reference:
#Chan Zhou, Fenglou Mao, Yanbin Yin, Jinling Huang, JohannPeter Gogarten, Ying Xu,
#AST: an automated sequence-sampling method for improving the taxonomic diversity of gene phylogenetic trees, 2013
#(submitted)
#Last updated: 2/2013
#Pls read README before using this script.

our $VERSION="version 0.4";

sub distributeM;
use strict;
use POSIX;
use warnings; #02/04/13
use Getopt::Std;#02/04/13
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %opt;
my $USAGE =<<EOF;
Usage:
	[perl] $0 [-hvul] -f -m 
Syntax:
	-h	This help
	-v	Print version and quit
	-u	Flags: 0, or 1, or 2. 
		0: exclude all taxa whose annotations have "unclassified" as prefix;
		1 [default]: include taxa annotated as "unclassified xxx"; but exclude the taxonomy "unclassified"(taxonomic id: 12908);
		2: include all taxa whose annotations have "unclassified" as prefix
	-l	boolean flags: 0 or 1.
		0: exclude taxa annotated as "environmental samples (xxx)"
		1 [default]: include taxa annotated as "environmental samples (xxx)"
	-f	input file name. It is composed of the list of gene/protein IDs, taxonomic IDs, score of non redundant hits.Pls ref. to the README for details.
	-m 	integer number of sampled sequences

Examples :
	perl $0  -u 1 -l 1 -f test -m 30 
	or
	perl $0 -f test -m 28

NOTE: 
	1. pls put the nodes.dmp & names.dmp files in the current directory.
	2. the version of nodes.dmp should be consistent with the version of gi2tax which is used to map the gi to the corresponding taxonomic gi in the input file (gi2tax: gi_taxid_prot.dmp for proteins / gi_taxid_nucl.dmp for nucleotide sequences)\n";
EOF

sub main::HELP_MESSAGE {
	print $USAGE, $/;
	exit 1;
}

sub main::VERSION {
	print $VERSION,$/;
	exit 1;
}

getopt('hvulf:m:', \%opt);#02/04/13 $opt_u & $opt_e: boolean flags, store their values

if (exists ($opt{h})) {
	main::HELP_MESSAGE;
}

if (exists ($opt{v})) {
	main::VERSION;
}

if ((not exists $opt{f}) or (not exists $opt{m})) {
	main::HELP_MESSAGE;
}

my $unclassified=1;
my $environment_sample=1;
if(defined $opt{u}){
	$unclassified=$opt{u};
}

if(defined $opt{l}){
	$environment_sample=$opt{l};
}

################################################################################################
#Step -(1)Specify taxonomic groups (T1,...,TG), start from level 2 of NCBI taxonomy by default.
#     -(2)Calculate the distribution of n hits over G taxonomic groups (n1,...,nG)
#     -(3)Calculate (m1,...,mG)such that
#        -(a)mi=1, keep it(sample it)
#        -(b)mi>1, -if Ti has no subgroup, then sample all mi hits in group Ti,& stop
#                  -if Ti has subgroups, then take Ti's subgroups as groups and take mi as m.
#     -(4)Sample sequences
#################################################################################################

my $m = $opt{m};#02/05/2013

my %taxon2m = ();
my %taxon2n = ();   #:=taxon2counts // #record how many hits in each taxonomy.
my %stop_T;         #stop searching  taxonomy in this node---includes %sample_T
my %sample_T;       #sample hits in taxonomy T group & stop
my %taxon2sons = ();    #hash-array
my %taxon2name = ();    #hash

my $level = 1;          #test

open NETHITS, $opt{f}; #02/05/2013
chomp( my @nethits = <NETHITS> );
my $numnethits = @nethits;
close NETHITS;

#2/6/13 ##########################################################
my $getue_state=`grep -P 'unclassified|environmental' names.dmp > names_ue.dmp`;
CHECK_GREP: if(not $getue_state){#indicate ---- finish "grep" command
	#warn "finish extracting the unclassified & environmental taxa from names.dmp\n";
}else{
	sleep(60);
	goto CHECK_GREP;
}

open UENAMEDMP, "<names_ue.dmp";
my $uenamedmp_line;
my %taxon2uenames;
my %exclude_taxa;
while($uenamedmp_line=<UENAMEDMP>){
	my @fields=split /\|/, $uenamedmp_line;
	my $taxon=$fields[0];
	#warn "ue-taxon=$taxon\n";#test
	my $uename=$fields[1];
	#warn "ue-names=$uename\n"; #test
	if (not exists $taxon2uenames{$taxon}){
		$taxon2uenames{$taxon}=$uename;
	}

	if($unclassified == 0){#exclude all
		if(($uename=~/unclassified/) and 
			(not exists $exclude_taxa{$taxon}) )
		{ 
			$exclude_taxa{$taxon}=$uename;
		}
	}elsif(($unclassified == 1) and ($taxon=~/12908/)){#exclude partial
			$exclude_taxa{$taxon}=$uename;
	}

	if($environment_sample == 0){#exclude all
		if(($uename=~/environmental\ssamples/) 
				and (not exists $exclude_taxa{$taxon}))
		{
			$exclude_taxa{$taxon}=$uename;
		}
	}
}#-end-for(my $i=0;$i<@uenames;$i++){
close UENAMEDMP;


##################################################################################
#Step -(1) Specify taxonomy groups (T1,...,TG) in each taxonomic level, 
#start from level 2 of NCBI taxonomy by default.
##################################################################################

if ( $m > $numnethits ) {    ###output all hits in $non_redundant_file_name
    open SAMPOF, ">out.$m.sample";
    for ( my $i = 0 ; $i < $numnethits ; $i++ ) {
        print SAMPOF "$nethits[$i]\n";
    }
    close SAMPOF;

}
else {                       #$m<=$numnethits

    open NODES, "<nodes.dmp";#2/5/13
    my $node_line;
#    chomp( my @nodes = <NODES> );
#    my $numnodes = @nodes;
    my %taxon2parent = ();
#    for ( my $i = 0 ; $i < $numnodes ; $i++ ) {
    while($node_line=<NODES>){
        if ( $node_line =~ /^(\d+)\s+\|\s+(\d*)\s+\|\s+(\w+)/ ) {
            $taxon2parent{$1} = $2;
            push( @{ $taxon2sons{$2} }, $1 );
            $taxon2name{$1} = $3;
        }
    }
    close NODES;

#exclude leaf-nethit which is the offsprings of any taxonomy in %exclude_taxa
#2/6/13
    for ( my $k = 0 ; $k < $numnethits ; $k++ ) {
		my @fields   = split /\t/, $nethits[$k];
		my $taxid    = $fields[1];
	        my $parent = "NA";
	        my $depth  = 1;                     #test
		my $taxon=$taxid; #tip taxon of $hitid
	        until ( ($taxon eq 1) or (exists $exclude_taxa{$taxon} ) ){ 
	            if ( exists $taxon2parent{$taxon} ) {
                	$parent = $taxon2parent{$taxon};
	                $taxon = $parent;
		    } else {
	                die "Cannot find taxonomy $taxon in this nodes.dmp. NOTE: the versions of gi2tax version & nodes.dmp should be consistent.(gi2tax: gi_taxid_prot.dmp / gi_taxid_nucl.dmp\n";
	            }
		}#end-until
		if(exists $exclude_taxa{$taxon}){
			delete $nethits[$k];
			splice (@nethits,$k,1);
		}
    } #end-for -$k

    my %nethit2taxon    = ();
    my %nethit2bitscore = ();
    my %taxon2nethit    = ();    #hash of array

    for ( my $i = 0 ; $i < $numnethits ; $i++ ) {
        my @fields   = split /\t/, $nethits[$i];
        my $hitid    = $fields[0];
	my $taxid    = $fields[1];#12/20/2012
        my $bitscore = $fields[2]; #12/20/2012
	
	if (not exists $exclude_taxa{$taxid}){ #-u -e option
	        $nethit2taxon{$hitid} = $taxid; #$fields[0]:=gi; $fields[2]:=taxonomy id
       		$nethit2bitscore{$hitid} = $bitscore;

       		 if ( exists $taxon2nethit{$taxid} ) {
	       	    push @{ $taxon2nethit{$taxid} }, $hitid;
       	 	 } else { #-u -e options
        	    $taxon2nethit{$taxid}[0] = $hitid;
       		 } 
	}#end-if-not-exists $exclude_taxa

    }#end-for-$i=0

    my $nethit;
    my $taxon;


#################################################################################
#Step - (2) calculate the distribution of n non-redundant hits in each taxonomy,
#then do the iteration for each taxonomic level, 
#e.g. level 3 - child of level 2 ( with taxonomic id: 131567-cellular organism)
#################################################################################
    open PRO, ">cal.process";#2/5/13; intermediate values and parameters

    foreach $nethit ( keys %nethit2taxon ) {
        $taxon = $nethit2taxon{$nethit};    #deepest leaf-node
        my $parent = "NA";
        my $depth  = 1;                     #test
        until ( $taxon eq 1  ) { 
            if   ( exists $taxon2n{$taxon} ) { $taxon2n{$taxon}++; }
            else                             { $taxon2n{$taxon} = 1; }

            if ( exists $taxon2nethit{$taxon} ) {
                push @{ $taxon2nethit{$taxon} }, $nethit;
            }
            else { $taxon2nethit{$taxon}[0] = $nethit; }

            if ( exists $taxon2parent{$taxon} ) {
                $parent = $taxon2parent{$taxon};
                $depth++;                   #test
                    print PRO "depth=$depth\ttaxon=$taxon\tparent=$parent\n"; #test
                $taxon = $parent;
            }
            else {
                die "Cannot find taxonomy $taxon in this nodes.dmp. NOTE: the versions of gi2tax version & nodes.dmp should be consistent.(gi2tax: gi_taxid_prot.dmp / gi_taxid_nucl.dmp\n";
            }
        }    #end-until
	#....
    }    #end-foreach $nethit

###################################################################################
#Step - (3) Calculate the distribution of (m1,m2, ..., mG) in each taxonomic level 
#           start from level 2 -- "cellular organism" with taxonomic id: 131567; 
#           then use %taxon2nethit to find the hits
###################################################################################

    $taxon = 131567 ; #start from level-2-node--taxonomy-id:131567, "cellular organism";
    my @taxons; #array of taxons install the data required to go into searching in the current level.
    $taxons[0] = $taxon;    #initial
    my $rarray_taxons;      #ref to array of @array_of_subtaxons;
    my $r_dis_m;            #ref to @dis_m == @m;

    my @array_of_taxons;    #array of @subtaxons//array of arrays
    my @m;

    ( $rarray_taxons, $r_dis_m ) = &distributeM( @taxons, $m );    #starting point
    @array_of_taxons = @$rarray_taxons;
    @m               = @$r_dis_m;
    my $num_array_of_taxons = @array_of_taxons;

    until ( $num_array_of_taxons == 0 ) {
        my @tmp_array_of_taxons = ();
        my @tmp_m               = ();

        for ( my $i = 0 ; $i < $num_array_of_taxons ; $i++ ) {
            @taxons = @{ $array_of_taxons[$i] };
            if ( @taxons > 0 ) {    #in case m==1, clear its returned sub@taxons
                ( $rarray_taxons, $r_dis_m ) = &distributeM( @taxons, $m[$i] );
                if ( defined @$rarray_taxons ) {
                    @tmp_array_of_taxons =
                      ( @tmp_array_of_taxons, @$rarray_taxons );
                    @tmp_m = ( @tmp_m, @$r_dis_m );
                }
            }   # end-if(@taxons > 0)
        }
        @array_of_taxons     = @tmp_array_of_taxons;
        @m                   = @tmp_m;
        $num_array_of_taxons = @array_of_taxons;
    }  # end-until


########################################################################################
#Step - (4) sample sequences after computing the distribution (m_i) over (T_i).
#           take m_i non-redundant hits from T_i according to their homolog score with the query seq.
#########################################################################################

    open TAX2M, ">out.$m.taxon2m";
    print TAX2M "Taxonomic_ID\tm\tn\ttaxonomic_name\n";
    open SAMPOF, ">out.$m.sample";
    print SAMPOF "Seq_ID\tTaxonomic_ID\tScore\n";

    my $sT;
    foreach $sT ( sort keys %sample_T ) {
        print TAX2M "$sT\t$taxon2m{$sT}\t$taxon2n{$sT}\t$taxon2name{$sT}\n";

        my $num_sT2nethit = @{ $taxon2nethit{$sT} };
        my @nethit_sort;    #descending according to their bitscores

        my @nethits_in_T      = @{ $taxon2nethit{$sT} };
        my %nethit_in_T2score = ();
        for ( my $w = 0 ; $w < @nethits_in_T ; $w++ ) {
            $nethit_in_T2score{ $nethits_in_T[$w] } =
              $nethit2bitscore{ $nethits_in_T[$w] };
        }
        my $j = 0;
        my $T;
        foreach $T ( sort { $nethit_in_T2score{$b} <=> $nethit_in_T2score{$a} }
            keys %nethit_in_T2score )
        {
            $nethit_sort[$j] = $T;
            $j++;
        }

        if ( $taxon2m{$sT} == 1 ) {
            my $sample_hit;
            if ( $num_sT2nethit == 1 ) {
                $sample_hit = $taxon2nethit{$sT}[0];
            }
            elsif ( $num_sT2nethit > 1 ) {

		    #sort-hash-by-value
                $sample_hit = $nethit_sort[0];
            }
            elsif ( $num_sT2nethit == 0 ) {
                print PRO "Error: no hit in taxonomy $sT ?\n";
            }
            print SAMPOF "$sample_hit\t$sT\t$nethit2bitscore{$sample_hit}\n";
        }
        elsif ( $taxon2m{$sT} > 1 ) {
            my @sample_hit;
            for ( my $i = 0 ; $i < $taxon2m{$sT} ; $i++ ) {
                $sample_hit[$i] = $nethit_sort[$i];
                print SAMPOF
                  "$sample_hit[$i]\t$sT\t$nethit2bitscore{$sample_hit[$i]}\n";
            }

        }   # end-if-elsif($taxon2m-$sT
    }   # end-foreach(%sample_T){
    close SAMPOF;
    close TAX2M;

}  # end-if($m>$numnethits)-else-{}



#####################################################################################
#----subroutine---------------------------------------------------------------------#
#----pass two values: @taxons-(contains the information of G),m---------------------#
#----return the references of two arrays: @array_of_subtaxons, @m ------------------#
#-----------------(value required to be distributed in the next level---------------#
#####################################################################################

sub distributeM {

    #It deals with three cases: (1) G = 1; (2) G = 2; (3) G > 2

    my $m      = pop(@_);
    my @taxons = @_;
    print PRO "subroutine \@taxons=@taxons\n";
    print PRO "subroutine \$m=$m\n";    #test
    my $G = @taxons;
    my $taxon;

    my @T;
    my @sort_T;
    my @array_of_e_sub_taxons_T;    #array of @e_sub_taxons_T;
    my @dis_m;    #distribute $m over @Ts, ($m0,$m1,...$m...)--->@dis_m;
    my $refarray_taxons = \@array_of_e_sub_taxons_T;
    my $ref_dis_m       = \@dis_m;

    print PRO "G=$G\n";    #test
    my %c_taxon2n = ();    #current taxonomy 2 n
    for ( my $i = 0 ; $i < @taxons ; $i++ ) {
        $c_taxon2n{ $taxons[$i] } = $taxon2n{ $taxons[$i] };
    }

    #--------begin-test-----------
    my $key;                  
    print PRO "\%c_taxon2n:\n";    
    foreach $key ( sort keys %c_taxon2n ) {    
        print PRO  "$key\t$c_taxon2n{$key}\n";       
    }
    #--------end-of-test----------
    
    #Case 1#######################################################################
    if ( $G == 1 ) {
        $taxon = $taxons[0];
        my @e_subtaxons = ()
          ; # effective sub taxons are these which have hits//== undef @e_subtaxons

        if ( $taxon2n{$taxon} < $m or $taxon2n{$taxon} == $m ) {    #$n<=$m;
            $taxon2m{$taxon}  = $taxon2n{$taxon};
            $sample_T{$taxon} = 0;
            $stop_T{$taxon}   = 0;
            undef @e_subtaxons;           #revised on 02/24/2011
            print PRO  "sample_T: $taxon\n";    #02/24/2011
        }
        else {                            #$n > $m
            $taxon2m{$taxon} = $m;
            if ( not( exists $taxon2sons{$taxon} ) ) {
                $sample_T{$taxon} = 0;
                $stop_T{$taxon}   = 0;
                undef @e_subtaxons;
                print PRO  "sample_T: $taxon\n";    #02/24/2011
            }
            else {                            #existing ...
                my @subtaxons = @{ $taxon2sons{$taxon} };
                for ( my $i = 0 ; $i < @subtaxons ; $i++ ) {
                    if ( exists $taxon2n{ $subtaxons[$i] } ) {
                        push @e_subtaxons, $subtaxons[$i];
                    }
                }

                if ( $m == 0 ) {
                    $stop_T{$taxon} = 0;
                    undef @e_subtaxons;
                }
                elsif ( $m == 1 ) {
                    $sample_T{$taxon} = 0;
                    $stop_T{$taxon}   = 0;
                    undef @e_subtaxons;
                    print PRO  "sample_T: $taxon\n";    #02/24/2011
                }
                elsif ( $m > 1 and @e_subtaxons == 0 )
                {    #the current taxon has no effective sub taxons
                    $sample_T{$taxon} = 0;
                    $stop_T{$taxon}   = 0;
                    undef @e_subtaxons;
                    print PRO  "sample_T: $taxon\n";    #02/24/2011
                }
                elsif ( $m < 0 ) { die "err \$m=$m < 0\n"; }    #$m<0
            }    #end-if-not-exists $taxon2sons{$taxon}

        }    #end-if($taxon2n{$taxon}<$m  or == )
        push @array_of_e_sub_taxons_T, [@e_subtaxons];
        $dis_m[0] = $m;
        my $nextG = @e_subtaxons;
        print PRO  "For G=1, nextG=$nextG\n";

    }

    #Case 2##############################################################################
    elsif ( $G == 2 ) {
        if ( $taxon2n{ $taxons[0] } < $taxon2n{ $taxons[1] } ) {
            $T[0] = $taxons[0];
            $T[1] = $taxons[1];
        }
        else {
            $T[0] = $taxons[1];
            $T[1] = $taxons[0];
        }

        my @n;
        $n[0] = $taxon2n{ $T[0] };
        $n[1] = $taxon2n{ $T[1] };    #n0<n1
        print PRO  "n0=$n[0]\tn1=$n[1]\n";  #test

        my @array_of_sub_taxons_T;    #array of @sub_taxons_Ts;
        if ( exists $taxon2sons{ $T[0] } ) {
            @{ $array_of_sub_taxons_T[0] } = @{ $taxon2sons{ $T[0] } };
        }
        else {                        #no sub-taxons
            @{ $array_of_sub_taxons_T[0] } = ();    #==undef
        }

        if ( exists $taxon2sons{ $T[1] } ) {
            @{ $array_of_sub_taxons_T[1] } = @{ $taxon2sons{ $T[1] } };
        }
        else {
            @{ $array_of_sub_taxons_T[1] } = ();    #==undef
        }

        for ( my $j = 0 ; $j < @array_of_sub_taxons_T ; $j++ ) {
            my @sub_taxons_T = @{ $array_of_sub_taxons_T[$j] };
            if ( @sub_taxons_T > 0 ) {
                for ( my $i = 0 ; $i < @sub_taxons_T ; $i++ ) {
                    if ( exists $taxon2n{ $sub_taxons_T[$i] } ) {
                        push @{ $array_of_e_sub_taxons_T[$j] },
                          $sub_taxons_T[$i];
                    }
                }
            }
            else {    #no-sub-taxons
                @{ $array_of_e_sub_taxons_T[$j] } = ();    #==undef
            }
        }

        #compute $dis_m[0], $dis_m[1];
        my @nextGs = ( '0', '0' );
        if ( defined @{ $array_of_e_sub_taxons_T[0] } ) {
            $nextGs[0] = @{ $array_of_e_sub_taxons_T[0] };
        }    #02/15/2011
        if ( defined @{ $array_of_e_sub_taxons_T[1] } ) {
            $nextGs[1] = @{ $array_of_e_sub_taxons_T[1] };
        }                             #02/15/2011
        print PRO  "\@nextGs=@nextGs\n";    #02/24/2011
        if ( ( $n[0] + $n[1] ) < $m or ( $n[0] + $n[1] ) == $m ) {    #(3-0)
            $dis_m[0]          = $taxon2n{ $T[0] };
            $dis_m[1]          = $taxon2n{ $T[1] };
            $sample_T{ $T[0] } = 0;
            $sample_T{ $T[1] } = 0;
            $stop_T{ $T[0] }   = 0;
            $stop_T{ $T[1] }   = 0;
            undef @{ $array_of_e_sub_taxons_T[0] };
            undef @{ $array_of_e_sub_taxons_T[1] };
            print PRO  "sample_T: $T[0]\t$T[1]\n";    #02/24/2011
        }
        else {
            if (    ( $n[0] < $m / 2 or $n[0] == $m / 2 )
                and ( $m / 2 < $n[1] or $n[1] == $m / 2 ) )
            {                                   #(3-1)
                $dis_m[0] = $n[0];
                $dis_m[1] = $m - $dis_m[0];
            }
            elsif ( $m / 2 < $n[0] or $m / 2 == $n[0] ) {    #(3-2)
                print PRO  "subroutine if(floor($m/2) == $m/2\n";    #test
                if ( floor( $m / 2 ) == $m / 2 ) {             #integer
                    $dis_m[0] = $m / 2;
                    $dis_m[1] = $m - $dis_m[0];                #02/15/2011
                }
                else {                                         #not integer
                    if (   ( $n[0] == $n[1] and $nextGs[0] < $nextGs[1] )
                        or ( $n[0] < $n[1] ) )
                    {
                        $dis_m[0] = floor( $m / 2 );
                        $dis_m[1] = $m - $dis_m[0];
                    }
                    elsif ( $n[0] == $n[1] and $nextGs[0] > $nextGs[1] ) {
                        $dis_m[1] = floor( $m / 2 );
                        $dis_m[0] = $m - $dis_m[1];
                    }
                    elsif ( $n[0] == $n[1] and $nextGs[0] == $nextGs[1] )
                    {                                          #02/15/2011
                        my $random_num = rand();
                        if ( $random_num < 0.5 ) {
                            $dis_m[0] = floor( $m / 2 );
                            $dis_m[1] = $m - $dis_m[0];
                        }
                        else {
                            $dis_m[1] = floor( $m / 2 );
                            $dis_m[0] = $m - $dis_m[1];
                        }
                    }
                }    #end-if(floor....
            }    #end-if-elsif(($n0<$m/2....

            print PRO  "m:\t@dis_m\n";    #test
            print PRO  "\@nextlevel_taxons_array[0]=@{$array_of_e_sub_taxons_T[0]}\n"
              ;                     #test
            print PRO  "\@nextlevel_taxons_array[1]=@{$array_of_e_sub_taxons_T[1]}\n"
              ;                     #test
        }    #end-if-else-(($n0+$n1)<$m or ($n0+$n1)== $m){#(3-0)
        $taxon2m{ $T[0] } = $dis_m[0];
        $taxon2m{ $T[1] } = $dis_m[1];

        for ( my $j = 0 ; $j < $G ; $j++ ) {
            print PRO  "\$j==$j\n";    #test
            if ( $dis_m[$j] == 0 ) {
                $stop_T{ $T[$j] } = 0;
                undef @{ $array_of_e_sub_taxons_T[$j] };
            }
            elsif (( $dis_m[$j] == 1 )
                or ( $dis_m[$j] > 1 and $nextGs[$j] == 0 )
                or ( $taxon2m{ $T[$j] } == $taxon2n{ $T[$j] } ) )
            {
                $sample_T{ $T[$j] } = 0;
                $stop_T{ $T[$j] }   = 0;
                undef @{ $array_of_e_sub_taxons_T[$j] };
                print PRO  "sample_T: $T[$j]\n";    #02/24/2011
            }
        }    #end-for-$j

    }

    #Case 3####################################################################################
    elsif ( $G > 2 ) {    #i.e.$G>=3

        #sorting..(n)i ascendingly
        my @m;
        my $k = -1;

        #ascendlingly sort by hash value "n"
        my $z;
        foreach $z ( sort { $c_taxon2n{$a} <=> $c_taxon2n{$b} }
            keys %c_taxon2n )
        {
            $k++;
            $sort_T[$k] = $z;
        }

        print PRO  "\@sort_T=@sort_T\n";    #test

	#------------------rank taxa in the next level if this taxa has sub-taxa-----------
        my @array_of_sub_taxons_T;

	#array of arrays: T-sorted taxonomy; value:-array-son taxonomys of the key-taxonomy
        for ( $k = 0 ; $k < $G ; $k++ ) {
            if ( exists $taxon2sons{ $sort_T[$k] } ) {
                @{ $array_of_sub_taxons_T[$k] } =
                  @{ $taxon2sons{ $sort_T[$k] } };
            }
            else {                    #no sub-taxons
                @{ $array_of_sub_taxons_T[$k] } = ();
            }
        }
        for ( $k = 0 ; $k < $G ; $k++ ) {
            my @sub_taxons_T = @{ $array_of_sub_taxons_T[$k] };
            if ( @sub_taxons_T > 1 ) {
                for ( my $i = 0 ; $i < @sub_taxons_T ; $i++ ) {
                    if ( exists $taxon2n{ $sub_taxons_T[$i] } ) {
                        push @{ $array_of_e_sub_taxons_T[$k] },
                          $sub_taxons_T[$i];
                    }
                }
            }
            else {    #no sub-taxons
                @{ $array_of_e_sub_taxons_T[$k] } = ();    #== undef
            }
        }    #end-for-$k

        #compare the distribution of @{$array_of_e_sub_taxons_T[$k]};
        my @nextGs;
        my $N = 1000000;
        my @nextG_v;    #variation: G-->N*G+index; keep rank
                        #to avoid that some $nextGs[$i] have the same value.
        for ( $k = 0 ; $k < $G ; $k++ ) {
            if ( defined @{ $array_of_e_sub_taxons_T[$k] } ) {
                $nextGs[$k] = @{ $array_of_e_sub_taxons_T[$k] };
            }
            else { $nextGs[$k] = 0; }    # it's possible that $nextGs[$k] == 0
            $nextG_v[$k] = $nextGs[$k] * $N + $k;
        }
        my @sort_nextGs  = sort { $b <=> $a } @nextGs;    #descending
        my @sort_nextG_v = sort { $b <=> $a } @nextG_v;

        my %index;
        @index{@nextG_v} = ( 0 .. $#nextG_v );

        #choose (m-G'*[m/G']) num of Ti which have hits distributing no less
        #narrow over the subgroups than others
        my @top_k;
        for ( my $i = 0 ; $i < $G ; $i++ ) {
            push @top_k, $index{ $sort_nextG_v[$i] };
        }
        $level++;                                         #test
        print PRO  "level=$level\n";                            #test
        print PRO  "nextlevel_taxonomy:\n";                     #test
        #------------------end---rank the  of next level_taxa-------------------------------
        #------------------------------------------------------------------------------------
        #####-----test-----#####
        print PRO  "\@nextG=@nextGs\n";
        print PRO  "\@sort_nextGs=@sort_nextGs\n";
        print PRO  "\@sort_nextG_v=@sort_nextG_v\n";
        print PRO  "\%index=%index\n";
        print PRO  "\@top_k=@top_k\n";

        #
        my $sum_n = 0;
        for ( $z = 0 ; $z < @taxons ; $z++ ) {
            $sum_n += $taxon2n{ $taxons[$z] };
        }

        print PRO  "sum_n=$sum_n\tcurrent_m=$m\n";
        if ( ( $sum_n < $m ) or ( $sum_n == $m ) ) {

            #(4-0)
            for ( my $l = 0 ; $l < @sort_T ; $l++ ) {
                $taxon2m{ $sort_T[$l] }  = $taxon2n{ $sort_T[$l] };
                $dis_m[$l]               = $taxon2n{ $sort_T[$l] };
                $stop_T{ $sort_T[$l] }   = 0;
                $sample_T{ $sort_T[$l] } = 0;
                undef @{ $array_of_e_sub_taxons_T[$l] };
                print PRO  "sample_T: $sort_T[$l]\n";    #02/24/2011
            }
        }
        elsif (( $c_taxon2n{ $sort_T[0] } > $m / $G )
            or ( $c_taxon2n{ $sort_T[0] } == $m / $G ) )
        {

            #(4-1)
            if ( floor( $m / $G ) == $m / $G ) {    #integer
                for ( my $l = 0 ; $l < $G ; $l++ ) {
                    $taxon2m{ $sort_T[$l] } = $m / $G;
                    $dis_m[$l] = $m / $G;
                }
            }
            else {                                  #non-integer
                my $big_mi       = floor( $m / $G ) + 1;
                my $num_big_mi   = $m - $G * floor( $m / $G );
                my $small_mi     = floor( $m / $G );
                my $num_small_mi = $G - $num_big_mi;

                for ( my $x = 0 ; $x < $num_big_mi ; $x++ ) {
                    $taxon2m{ $sort_T[ $top_k[$x] ] } = $big_mi;
                    $dis_m[ $top_k[$x] ] = $big_mi;
                }
                for ( my $y = $num_big_mi ; $y < $G ; $y++ ) {
                    $taxon2m{ $sort_T[ $top_k[$y] ] } = $small_mi;
                    $dis_m[ $top_k[$y] ] = $small_mi;
                }
            }
            ###end--(4-1)
        }
        else {    #find index t,such that n_{t} <= m/G<n_{t+1}, 0<=t<=(G-2)
                  #(4-2)
            my $t;
            for ( $k = 0 ; $k < ( $G - 1 ) ; $k++ ) {
                if (
                    (
                           ( $c_taxon2n{ $sort_T[$k] } < $m / $G )
                        or ( $c_taxon2n{ $sort_T[$k] } == $m / $G )
                    )
                    and ( $m / $G < $c_taxon2n{ $sort_T[ $k + 1 ] } )
                  )
                {    #revised on 02/22/2011
                    $t = $k;
                    last;
                }
            }

            #-----take m_i=n_i for 0<=i<=t;
            for ( my $i = 0 ; $i <= $t ; $i++ ) {
                $dis_m[$i] = $c_taxon2n{ $sort_T[$i] };
                $taxon2m{ $sort_T[$i] } = $dis_m[$i];
            }

            my $G_plus = $G - ( $t + 1 );    #previous index: 0,1,2,...,t
            my @T_plus;
            my $m_plus = $m;                 #initial value
            my $u;
            for ( $u = 0 ; $u <= $t ; $u++ ) {
                $m_plus = $m_plus - $dis_m[$u];
            }
            for ( $u = ( $t + 1 ) ; $u < $G ; $u++ ) {
                push @T_plus, $sort_T[$u];
            }

            #---------------------
            if ( $G_plus > 0 and $m_plus > 0 ) {
                print PRO  "plus: \@T_plus=@T_plus\t\$m_plus=$m_plus\n";
                my ( $rarray_T_plus, $r_dis_m_plus ) =
                  &distributeM( @T_plus, $m_plus );

                for ( my $i = $t + 1 ; $i < $G ; $i++ ) {
                    my $w = $i - ( $t + 1 );
                    print PRO  "\@T_plus to m ($i, $w): $taxon2m{$T_plus[$w]}\n";
                    $dis_m[$i] = $taxon2m{ $T_plus[$w] };
                }

            }
            elsif ( $m_plus < 0 or $m_plus == 0 ) {
                die "debug: \$m_plus = $m_plus, it should be >0\n";
            }
            for ( my $j = 0 ; $j < $G ; $j++ ) {
                $taxon2m{ $sort_T[$j] } = $dis_m[$j];
            }    # revised on 02/23/2011
        }    #end-if($sum_n)-elsif...

        for ( my $j = 0 ; $j < $G ; $j++ ) {

            #	print PRO  "\$j==$j\n";#test
            if ( $dis_m[$j] == 0 ) {
                $stop_T{ $sort_T[$j] } = 0;
                undef @{ $array_of_e_sub_taxons_T[$j] };
            }
            elsif (( $dis_m[$j] == 1 )
                or ( $dis_m[$j] > 1 and $nextGs[$j] == 0 )
                or ( $taxon2m{ $sort_T[$j] } == $taxon2n{ $sort_T[$j] } ) )
            {
                $sample_T{ $sort_T[$j] } = 0;
                $stop_T{ $sort_T[$j] }   = 0;
                undef @{ $array_of_e_sub_taxons_T[$j] };
                print PRO  "sample_T: $sort_T[$j]\n";    #02/24/2011
            }
        }    #end-for-$j;

    }    #end-if($G == 1)...$G==2, $G>2

    #---begin test----
    if    ( $G == 1 ) { print PRO  "\@sort_T=@taxons\n"; }
    elsif ( $G == 2 ) { print PRO  "\@sort_T=$T[0]\t$T[1]\n"; }
    else              { print PRO  "\@sort_T=@sort_T\n"; }
    print PRO  "\@dis_m=@dis_m\n\n\n";    
    #---end test----

    return ( $refarray_taxons, $ref_dis_m );  
}    #end-subroutine
