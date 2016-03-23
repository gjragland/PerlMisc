#!usr/bin/perl -w
# GJR 10/3/2014
#full annotation of Trinity-based transcriptome assembly using .xml blastx output from searches of:
  #1. Uniprot swiss prot
  #2. Uniprot trmbl
  #3. Flybase translation database
#uses web query to retrieve additional details (e.g., GO) from uniprot online database
#can't figure out programmatic access to flybase, so need to batch-download information on hits, then use that file for flybase annotation

use strict;
use LWP::Simple;

my $sprotfile="allSprot.xml";
#my $sprotfile="sprottest.xml";
my $trmblfile="allTrmbl.xml";
#my $trmblfile="trmbltest.xml";
my $flybasefile="flybaseResults.all.xml";
#my $flybasefile="flybasetest.xml";
my $uniprotSpecList="speclist.txt";
my $flybaseAnnoFile="FlybaseAnnoFromWeb.txt";
my $fastafile="Trinity.fasta";
my $outfile="Trinity.fasta.anno.txt";

my %nonEukaryotes;
open IN, "<$uniprotSpecList";
while (<IN>) {
  if (/^\S+\s+(\S)\s+\d+\:\s+N\=(.*$)/) {
    my $taxa=$1;
    my $species=$2;
    $species =~ s/\(.*\)//;
    chomp $species;
    $nonEukaryotes{$species}=$taxa if $taxa !~ m/E/;
  }
}
close IN;

my %flybaseAnno;
my $flybaseHeader;
open IN, "<$flybaseAnnoFile";
while (<IN>) {
  if ($.==1) {
    $flybaseHeader=$_;
    $flybaseHeader =~ s/^\#//;
    chomp $flybaseHeader;
    next;
  }
  $_ =~ s/<newline>/\;/g;
  my @vals=split "\t";
  $flybaseAnno{$vals[0]}=$_;
}
close IN;


my $sprotHitsRef=hits_from_blast($sprotfile,"uniprot");
my $trmblHitsRef=hits_from_blast($trmblfile,"uniprot");
my $flybaseHitsRef=hits_from_blast($flybasefile,"flybase");

my @ids;
open IN, "<$fastafile";
while (<IN>) {
  if (/^>(\S+)/) {
    push @ids, $1;
  }
}
close OUT;


open OUT, ">$outfile";
print OUT "Flag\tQuery\tUniprotEval\tUniprotAcc\tUniprotName\tStatus\tProteinNames\tGeneNames\tOrganism\tLength\tGeneOntology\tGO-IDs\tPathway\tInterPro\tflybaseEval\t$flybaseHeader";
for my $id (@ids) {
  my $out;
  my $nonEukScore=0;
  my $idRef;
  my $print=0;
  if ($idRef=$sprotHitsRef->{$id}) {
    $print=1;
    my($eval,$hitId, $species)=split("\t",shift(@$idRef));
    $out=retrieve_uniprot($hitId);
    $out="$id\t$eval\t$out";
    $nonEukScore=2 if exists $nonEukaryotes{$species};
    while ( my $info = shift @$idRef ) {
      my($toss,$toss2, $spec)=split("\t",$info);
      $nonEukScore++ if exists $nonEukaryotes{$species};
    }
  } elsif ($idRef=$trmblHitsRef->{$id}) {
    $print=1;
    my($eval,$hitId, $species)=split("\t",shift(@$idRef));
    $out=retrieve_uniprot($hitId);
    $out="$id\t$eval\t$out";
    $nonEukScore=2 if exists $nonEukaryotes{$species};
    while ( my $info = shift @$idRef ) {
      my($toss,$toss2, $spec)=split("\t",$info);
      $nonEukScore++ if exists $nonEukaryotes{$species};
    }
  } else {$out="$id".("\tNA") x 12}
  #replace missing data with NA
  #need match loop to account for overlapping matches (can't use /g)
  $out =~ s/\t\t/\tNA\t/ while $out =~ m/\t\t/;
  $out =~ s/\t\n/\tNA/g;
  chomp $out;
  if ($idRef=$flybaseHitsRef->{$id}) {
    $print=1;
    my($eval,$hitId)=split("\t",shift(@$idRef));
    if (my $flyinfo=$flybaseAnno{$hitId}) {
      $flyinfo =~ s/\t\-/\tNA/ while $flyinfo =~ m/\t\-/;
      $flyinfo =~ s/\t\n/\tNA/g;
      chomp $flyinfo;
      $out=$out."\t$eval\t$flyinfo";
    } else {$out=$out.("\tNA") x 10}  
  } else {$out=$out.("\tNA") x 10}
  my $flag="none";
  $flag="nonEuk" if $nonEukScore >= 2;
  print OUT "\n$flag\t$out" if $print==1;
}
close OUT;

#my $out=retrieve_uniprot('Q05893');


#-----subs-------------

#query uniprot database via http
sub retrieve_uniprot {
  my $query=shift;
  my $base='http://www.uniprot.org/uniprot/';
  #my $query='Q05893';
  my $format='tab';
  #available annotation values available at http://www.uniprot.org/faq/28#batch_retrieval_of_entries
  my $cols='id,entry%20name,reviewed,protein%20names,genes,organism,length,go,go-id,pathway,interpro';
  #form http query
  $query="$base\?query=accession%3a$query\&format=$format\&columns=$cols";
  #perform query through LWP::Simple function 'get'
  my $content = get $query;
  if (defined $content) {
    #remove first line, which is header
    $content =~ s/^[^\n]+?\n//s;
  } else {$content="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"}
  return($content);
}


#reads the blastx outfile chunk by chunk (iteration by iteration), each chunk ending in the terminator
# takes blastx xml file name
# returns reference to hash mapping query id to subject hit id and e-value
# added extra paramater to designate database origin, which will change capturing lines
sub hits_from_blast {
  my $file=shift;
  my $database=shift;
  my $term = $/;
  my %idHash;
  open(XMLFILE, "<$file");
  $/ = "<\/Iteration>";	
  while(<XMLFILE>) {
    #capture query id
    if (/<Hit_num>/ and /<Iteration_query\-def>(\S+).+?</) {
      my $QueryId=$1;
      my $line;
      my $acc;
      my $eval;
      my $species;
      if ($database =~ m/uniprot/) {
	while (/.*?<Hit_def>.*?OS\=(.*?\=).*?<Hit_accession>(\S+?)<\/Hit_accession>.*?<Hsp_evalue>(\S+)<\/Hsp_evalue>/sg) {
	  $species=$1;
	  $acc=$2;
	  $eval=$3;
	  $species =~ s/\S\S\=//;
	  $species =~ s/\(.*\)//;
	  $species =~ s/\s+$//;
	  $idHash{$QueryId}=[] unless exists $idHash{$QueryId};
	  push @{$idHash{$QueryId}},"$eval\t$acc\t$species";
	}
      }
      if ($database =~ m/flybase/) {
	while (/.*?<Hit_def>.*?(FBgn\d+).*?<Hsp_evalue>(\S+)<\/Hsp_evalue>/sg) {
	  $acc=$1;
	  $eval=$2;
	  $idHash{$QueryId}=[] unless exists $idHash{$QueryId};
	  push @{$idHash{$QueryId}},"$eval\t$acc";
	}
      }
    }
  }
  close XMLFILE;
  $/ = $term;
  return(\%idHash);
}


#http://www.uniprot.org/uniprot/?query=accession%3aQ05893&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,go,go-id,pathway,interpro
#http://www.uniprot.org/uniprot/?query=accession%3aQ05893&format=tab&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,go,go-id,pathway,interpro
