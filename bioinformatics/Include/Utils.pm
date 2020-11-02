package Utils;


use strict 'vars';

##################################################
# split a bed file in n bed files based on chrom
# bed file format: chrom  start  end  rest
sub split_bed_by_chrom {
    my ($fileBed, $fileOut);
    my ($line, $chrom, $start, $end, $rest);
    my ($n_in, $n_out, $n_out_total);
    my ($last_chrom);

    $last_chrom = "";
    
    $fileBed = @_;

    print "open bed file $fileBed\n";
    open(FHI, $fileBed) or die("cannot open $fileBed:$!\n");
    while ($line = <FHI>) {
	next unless ($line =~ /^chr/);
	$n_in++;
	$line =~ s/\012//;
	$line =~ s/\015//;
	($chrom, $start, $end, $rest) = split(/\t/, $line);
	if ($chrom ne $last_chrom) {
	    if ($last_chrom ne "") {
		close(FHO);
		print "close $fileOut ($n_out lines)\n";
	    }
	    $n_out = 0;
	    $fileOut = $fileBed . "_$chrom";
	    print "open file $fileOut\n";
	    open(FHO, ">$fileOut") || die("cannot creat $fileOut");
	    $last_chrom = $chrom;
	}
	print FHO $chrom."\t".$start."\t".$end."\n";
	$n_out++;
	$n_out_total++;
    }
    close(FHI);
    print "in $n_in out $n_out_total\n";
    close(FHO);


}


1;
