use 5.014;
my %utr;
while(<>){
	if(/.+?\sexon\s(\d+?)\s(\d+?)\s.+?Target=(comp\d+?_c\d+?_seq\d+?)\s(\d+?)\s(\d+?)\s/){
		${$utr{$3}}{$4}=$5;
		}
	else{
		next;
	}
}

for my $id (keys %utr){
	if (keys %{$utr{$id}}>1){
		my @a; #exon vegek
		for my $start (sort { $a <=> $b } keys %{$utr{$id}}){
			push(@a,$utr{$id}{$start});
		}
		pop(@a);
		print("$id $a[-1]\n");
	}
}
say"";
