#!/usr/bin/perl


#given a bam file and a haplotype database file
#drop from the DB file entries (either in the seq dict
#or in the data itself records that don't appear in the bam's dict)

#DEFINE COMMAND-LINE ARGS FOR INPUT 
my $N_ARGS=$#ARGV ;
$N_ARGS++;
if($N_ARGS!=3) {
	die "Usage perl $0 <IN_BAM> <IN_HAPLO> <OUT_FILTERED_HAPLO>" ;
	}
#get input data
my $in_bam=$ARGV[0];
my $in_hap=$ARGV[1];
my $out_hap=$ARGV[2];

#load the seq dict items into a hash
my $seq_dict_cmd="samtools view -H $in_bam|grep -P '^\@SQ'|cut -f2|grep -Po ':.*'|tr -d \":\"";
#print "Using command $seq_dict_cmd\n";
my @seq_dict=`$seq_dict_cmd` ;
my %seq_hash;
for my $s (@seq_dict)
	{
	$s=trim($s);
	#print "got key $s\n" ;
	$seq_hash{$s}=1;
	}

#read the haplo file and
# for Seq dict entries (^@SQ) drop them if they're in the dict
# for data lines (not beginning with @) drop them if the seq is in the hash

open(HAPLOREAD,"<",$in_hap) or die "Error in opening $in_hap for reading !" ;
open(HAPLOWRITE,">",$out_hap) or die "Error in opening $out_hap for writing !" ;
while(<HAPLOREAD>)
	{
	my $line=$_;
	chomp($line);
	#print "got $line\n" ;
	my $emit_flag=1;
	if($line=~m/^\@SQ/)
		{
		my @pieces=split('\t',$line);
		my $seq_piece=$pieces[1];
		if($seq_piece=~m/^[^:]+:(.*)$/)
			{
			#print "from line $line got seq_piece $seq_piece\n";
			my $sn=$1;
			$sn=trim($sn);
			#print "got sn $sn \n" ;
			if(exists $seq_hash{$sn})
				{
				#print "hash exists\n";
				}
			else
				{
				#print "hash doesn't exists!\n" ;
				#NOT in the bam dict!
				$emit_flag=0;
				}
			}
		else
			{

			}
		}
	elsif(!($line=~m/^@/) && !($line=~m/^#/))
		{
		#if the line doesn't begin with at sign, then it's data!
		#print "this line is data : '$line'\n";
		my @pieces=split('\t',$line);
		my $seq=$pieces[0];
		$seq=trim($seq);
		#print "seq is $seq\n";
		if(exists $seq_hash{$seq})
			{
			#is in bam so go ahead and emit
			}
		else
			{
			#not in BAM !
			$emit_flag=0;
			}
		}

	print HAPLOWRITE $line."\n" unless $emit_flag==0;
	}

#close the IO channels
close(HAPLOWRITE);
close(HAPLOREAD);



sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
