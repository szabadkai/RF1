use 5.014;

open INTRON ,'last_intron_3utr.txt' or die('INTRON');
open CARYYA ,'Nb_caryya.txt'or die('CARYYA');

my %caryya;
my %intron;

while(<CARYYA>){
	if (/(comp\d+?_c\d+?_seq\d+?)\s(.+)/){
		$caryya{$1}=$2;
	}
}

while(<INTRON>){
	if (/(comp\d+?_c\d+?_seq\d+?)(.+)/){
		$intron{$1}=$2;
	}
}

for my $key (keys %caryya){
	
	print "$key\t$caryya{$key} intron: $intron{$key}\n";
	
}
