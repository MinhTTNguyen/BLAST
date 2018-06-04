=pod
30th April 2012
This script is to report the BLASTP result in tabular format as followed:
Query	Query_length	Hit_ID	Target_length	E_value	%Id	%Si	%query_coverage	%target_coverage	Query_start	Target_start	Query_end	Target_end
*new: Query sequences include those from Spoth2p4 and the collected AARSs
*new: add query coverage and target coverage

#Modified on January 10th 2014

#Modified on January 27th 2015
only print out first hits
use for any BLASTN v2.2.29 output file
=cut


#! C:\Perl\bin -w
use strict;

print "\nInput common path:";
my $path=<STDIN>;
chomp($path);

print "\nInput BLASTN folder:";
my $folder_in=<STDIN>;
chomp($folder_in);

my $folder_out=$folder_in."_best_hit";
mkdir "$path\\$folder_out";
chdir "$path\\$folder_out";

opendir(DIR, "$path\\$folder_in") || die "Cannot open folder: $folder_in";
my @files=readdir(DIR);
closedir(DIR);

shift(@files);shift(@files);

foreach my $filein (@files)
{
	my $fileout=substr($filein,0,-4);
	$fileout=$fileout."_best_hit.txt";

	open(In,"<$path\\$folder_in\\$filein") || die "Cannot open file $filein";
	open(Out,">$fileout") || die "Cannot open file $fileout";
	print Out "query\tquery_len\ttarget_id\ttarget_len\tquery_strand\ttarget_strand\te_value\tid\tgap\tquery_cov\tquery_start\tquery_end\ttarget_cov\ttarget_start\ttarget_end\n";
	my $query="";
	my $query_len="";
	my $target_id="";
	my $target_len="";
	my $e_value="";
	my $id="";
	my $gap="";
	
	my $query_start="";
	my $query_end="";
	my $target_start="";
	my $target_end="";
	my $query_cov="";
	my $target_cov="";
	
	my $query_strand="";
	my $target_strand="";
	my %hash_query_print_status=();
	while (<In>)
	{
		chomp($_);
		#Query= sp|P50475|SYAC_RAT Alanine--tRNA ligase, cytoplasmic OS=Rattus norvegicus GN=Aars PE=1 SV=3
		# or Query= Spoth2p4_006882_Arg_cyt
		if ($_=~/^Query\=/)
		{
			if ($query)
			{
				if ($target_id)
				{
					if ($query_len>0) 
					{
						if ($query_strand eq "Plus"){$query_cov=($query_end-$query_start+1)/$query_len;}
						else{$query_cov=($query_start-$query_end+1)/$query_len;}
					}else  {print "$query\t$query_len\t$target_id\t$target_len";exit;}
					
					if ($target_len>0)
					{
						if ($target_strand eq "Plus"){$target_cov=($target_end-$target_start+1)/$target_len;}
						else{$target_cov=($target_start-$target_end+1)/$target_len;}
					}else{print "$query\t$query_len\t$target_id\t$target_len";exit;}
					
					unless($hash_query_print_status{$query}){print Out "$query\t$query_len\t$target_id\t$target_len\t$query_strand\t$target_strand\t$e_value\t$id\t$gap\t$query_cov\t$query_start\t$query_end\t$target_cov\t$target_start\t$target_end\n";$hash_query_print_status{$query}=1;}
				}else {unless($hash_query_print_status{$query}){print Out "$query\t$query_len\tNo hits\n";$hash_query_print_status{$query}=1;}}
			}
			$query="";
			$query_len="";
			$target_id="";
			$target_len="";
			$e_value="";
			$id="";
			$gap="";
			$query_start="";
			$query_end="";
			$target_start="";
			$target_end="";
			$query_cov="";
			$target_cov="";
			
			$query_strand="";
			$target_strand="";
			
			if ($_=~/^Query\=\s*(.*)/){$query=$1;$query=~s/\s*//g;}
			
			else {print "Query is not as described!!!\n$_\n";}
		}
	
		#>Spoth2p4_011163
		#Length = 1118
		if ($_=~/^\>\s*(.*)/)
		{
			if ($target_id)
			{
				if ($query_len>0) 
				{
					if ($query_strand eq "Plus"){$query_cov=($query_end-$query_start+1)/$query_len;}
					else{$query_cov=($query_start-$query_end+1)/$query_len;}
				}else  {print "$query\t$query_len\t$target_id\t$target_len";exit;}
				
				if ($target_len>0)
				{
					if ($target_strand eq "Plus"){$target_cov=($target_end-$target_start+1)/$target_len;}
					else{$target_cov=($target_start-$target_end+1)/$target_len;}
				}else{print "$query\t$query_len\t$target_id\t$target_len";exit;}
					
				unless($hash_query_print_status{$query}){print Out "$query\t$query_len\t$target_id\t$target_len\t$query_strand\t$target_strand\t$e_value\t$id\t$gap\t$query_cov\t$query_start\t$query_end\t$target_cov\t$target_start\t$target_end\n";$hash_query_print_status{$query}=1;}
			}
			$target_id="";
			$target_len="";
			$e_value="";
			$id="";
			$query_start="";
			$query_end="";
			$target_start="";
			$target_end="";
			$query_cov="";
			$target_cov="";
			$query_strand="";
			$target_strand="";
			
			$target_id=$1;
			#print "Target ID: $target_id";exit;
		}
		#old BLASTP version
		#Query= Spoth2p4_006882_Arg_cyt
		#         (1176 letters) 
		# Below is the new version
		#Query= Spoth2p4_001239_Ala_cyt
		#Length=959
		#>Spoth2p4_011163
		#Length = 1118
		if ($_=~/^Length=(\d+)/)
		{
			my $temp=$1;
			if ($target_id) {$target_len=$temp;}
			else {$query_len=$temp;}
		}
		
		# Score =  512 bits (1318), Expect = e-145,   Method: Compositional matrix adjust. (old version)
		# Score = 1393 bits (3606),  Expect = 0.0, Method: Compositional matrix adjust. (new version)
		if ($_=~/.*\,\s*Expect\s*\=\s*(.+)/)
		{
			if ($e_value)
			{
				if ($query_len>0) 
				{
					if ($query_strand eq "Plus"){$query_cov=($query_end-$query_start+1)/$query_len;}
					else{$query_cov=($query_start-$query_end+1)/$query_len;}
				}else  {print "$query\t$query_len\t$target_id\t$target_len";exit;}
					
				if ($target_len>0)
				{						
					if ($target_strand eq "Plus"){$target_cov=($target_end-$target_start+1)/$target_len;}
					else{$target_cov=($target_start-$target_end+1)/$target_len;}
				}else{print "$query\t$query_len\t$target_id\t$target_len";exit;}
					
				unless($hash_query_print_status{$query}){print Out "$query\t$query_len\t$target_id\t$target_len\t$query_strand\t$target_strand\t$e_value\t$id\t$gap\t$query_cov\t$query_start\t$query_end\t$target_cov\t$target_start\t$target_end\n";$hash_query_print_status{$query}=1;}
			}
			$e_value="";
			$id="";
			$query_start="";
			$query_end="";
			$target_start="";
			$target_end="";
			$query_cov="";
			$target_cov="";
			$query_strand="";
			$target_strand="";
			$query_strand="";
			$target_strand="";
			
			$e_value=$1;
		}
	
		# Identities = 16/64 (25%), Positives = 26/64 (40%), Gaps = 6/64 (9%)
		# Identities = 10/19 (52%), Positives = 13/19 (68%)
		# new version:
		#  Identities = 511744/511744 (100%), Gaps = 0/511744 (0%)
		if ($_=~/^\s*Identities\s*\=\s*.*\((.+)\%\)\,\s*Gaps\s*\=\s*.*\((.+)\%\)/){$id=$1;$gap=$2;}
	
		# Strand=Plus/Plus
		if ($_=~/^\s*Strand\=(\w+)\/(\w+)/)
		{
			$query_strand=$1;
			$target_strand=$2;
		}
		
		# Query: 90  LVAVGDHASKQMVKFAANINKESIVDVEGVVRKVNQKIGS-CTQQDVELHVQKIYVISLA 148
		# Sbjct: 68  FVSLGDGSSLAPLQALVQADDAKDLAVGAAVRLTGSWVSSPGVAQSHELHVSRVEVLGPS 127
		# new:
		#Query  3    EIEWTGARVRKTFLDFFAERGHSIVPSSSVVPHNDPTLLFTNAGMNQFKPIFLGTIGKTE  62
		#Sbjct  6    EHRWSAPRVRQAFLDFFSQKEHTIVPSSSVVPHNDPTLLFTNAGMNQFKPVFLGTVAQSD  65
		if ($_=~/^Query\s*(\d+)\s*.*\s+(\d+)/) 
		{	
			unless ($query_start){$query_start=$1;}
			$query_end=$2;
		}	
	
		if ($_=~/^Sbjct\s*(\d+)\s*.*\s+(\d+)/) 
		{	
			unless ($target_start){$target_start=$1;}
			$target_end=$2;
		}
	}
	if ($target_id)
	{
		if ($query_len>0) 
		{
			if ($query_strand eq "Plus"){$query_cov=($query_end-$query_start+1)/$query_len;}
			else{$query_cov=($query_start-$query_end+1)/$query_len;}
		}else  {print "$query\t$query_len\t$target_id\t$target_len";exit;}
					
		if ($target_len>0)
		{
			if ($target_strand eq "Plus"){$target_cov=($target_end-$target_start+1)/$target_len;}
			else{$target_cov=($target_start-$target_end+1)/$target_len;}
		}else{print "$query\t$query_len\t$target_id\t$target_len";exit;}
		
		unless($hash_query_print_status{$query}){print Out "$query\t$query_len\t$target_id\t$target_len\t$query_strand\t$target_strand\t$e_value\t$id\t$gap\t$query_cov\t$query_start\t$query_end\t$target_cov\t$target_start\t$target_end\n";$hash_query_print_status{$query}=1;}
	}else {unless($hash_query_print_status{$query}){print Out "$query\t$query_len\tNo hits\n";$hash_query_print_status{$query}=1;}}
	close(In);
	close(Out);
}
