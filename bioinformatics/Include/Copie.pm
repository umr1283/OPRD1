package Copie;

sub run {
    my($dir_from, $from, $dir_to, $to) = @_;

    check($dir_from, $from, $dir_to, $to);

    #my @commands = ( "sleep 5 && echo step 1 done",
    #                 "sleep 10 && echo step 2 done",
    #                 "sleep 7 && echo step 3 done" );

    my @commands = ( "rsync -a $dir_from/$from/ $dir_to/$to/raw --log-file=$dir_to/$to/raw/rsync.log");

    my @pids;
    foreach my $cmd( @commands ) {
    my $pid = fork;
    if ( $pid ) {
        # parent process
        push @pids, $pid;
        next;
    }

    # now we're in the child
    system( $cmd );
    exit;            # terminate the child
    }

    wait for @pids;   # wait for each child to terminate

    print localtime()." : all done.\n";

    #print "$dir_from, $from, $dir_to, $to";
}

sub run_old_system {
    my($dir_from, $from, $dir_to, $to) = @_;

    check($dir_from, $from, $dir_to, $to);
    # count_all($dir_from, $from);
#    system("date");
    cp_all($dir_from, $from, $dir_to, $to);
    #print "copy from $dir_from/$from to $dir_to/$to\n";
#    system("date");
}

sub check {
    my ($dir_from, $from, $dir_to, $to) = @_;

#    print "check $dir_from, $from, $dir_to, $to\n";
    die("no dir_from $dir_from\n") unless (-d $dir_from);
    die("no dir_to $dir_to\n") unless (-d $dir_to);
    die("no from $dir_from/$from\n") unless (-d "$dir_from/$from");
}

sub cp_all {
    my ($dir_from, $from, $dir_to, $to) = @_;
    my ($tot);

#    print "cp_all $dir_from, $from, $dir_to, $to\n";
    mkdir "$dir_to" unless (-d $dir_to);
#   chmod(0775, $dir_to) or die "Couldn't chmod $dir_to: $!";
    mkdir "$dir_to/$to" unless (-d "$dir_to/$to");
#   chmod(0775, "$dir_to/$to") or die "Couldn't chmod $dir_to/$to: $!";
    mkdir "$dir_to/$to/raw" unless (-d "$dir_to/$to/raw");
#   chmod(0775, "$dir_to/$to/raw") or die "Couldn't chmod $dir_to/$to/raw: $!";
    $tot = cp_dir("$dir_from/$from", 1, "$dir_to/$to/raw");
}

sub cp_dir {
    my ($dir, $indent, $dest) = @_;
    my ($tot_file, $tot_dir);
    my ($tot_byte, $tot_K, $tot_M, $tot_G);
    my ($node);
    my (@tmp, $file_size);
    my ($prefix);
    my (@size);
    my ($dh);
     
    $prefix = "." x $indent;

    # print "count_dir $dir\n";
    $tot_byte = 0;
    $tot_file = 0;
    $tot_dir = 0;
    opendir $dh, "$dir" or die("cannot open $dir :$!\n");
    while ($node = readdir $dh) {
	next if ($node eq "..");
	next if ($node eq ".");
	next if ($node eq "Processed");
	next if ($node =~ /.cif/);
	next if ($node =~ /.cnf/);
	# next if ($node =~ /s_.*tif/);
	if ($node =~ /s_.*tif/) {
	    if ("$dir/$node" =~ /Images\/L\d\d\d\/[CP]\d+\.\d\/s_\d.*\.tif/) {
		next;
	    } elsif ("$dir/$node" =~ /FindFocus_.*\.tif.gz/) {
		next; # image de calibration
	    } else {
#		print "keep tif file  $dir/$node\n";
	    }
	}
	if (-d "$dir/$node") {
	    # print "$prefix work on dir $dir/$node\n";
	    # print "mkdir $dest/$node", "\n";
	    mkdir "$dest/$node" unless (-d "$dest/$node");
	    $tot_byte += cp_dir("$dir/$node", $indent + 1, "$dest/$node");
	    $tot_dir++;
	} else {
#	    print "$prefix work on file $dir/$node\n";
	    if (-f "$dest/$node"){
#		print ("skip already exist $dest/$node\n");}
	    } else {
#	    print ("cp $dir/$node $dest/$node\n");
	    system("cp $dir/$node $dest/$node");}
	    $tot_file ++;
	    @tmp = stat("$dir/$node");
	    $file_size = $tmp[7];
	    $tot_byte += $file_size;
	}
    }
    closedir $dh;
    $tot_K = sprintf("%d", $tot_byte / 1024);
    $tot_M = sprintf("%d", $tot_K / 1024);
    $tot_G = sprintf("%d", $tot_M / 1024);
    if ($indent < 4) {
#	print "$prefix $dir : $tot_dir dirs, $tot_file files, $tot_byte bytes, $tot_K K, $tot_M M, $tot_G G\n";
    }
    return ($tot_byte);
}

sub count_all {
    my ($dir_from, $from) = @_;
    my ($tot);

#    print "count_all $dir_from, $from\n";
    $tot = count_dir("$dir_from/$from", 1);
#    print "count_all $from : $tot\n";
}

sub count_dir {
    my ($dir, $indent) = @_;
    my ($tot_file, $tot_dir);
    my ($tot_byte, $tot_K, $tot_M, $tot_G);
    my ($node);
    my (@tmp, $file_size);
    my ($prefix);
    my (@size);
    my ($dh);
     
    $prefix = "." x $indent;

    # print "count_dir $dir\n";
    $tot_byte = 0;
    $tot_file = 0;
    $tot_dir = 0;
    opendir $dh, "$dir" or die("cannot open $dir :$!\n");
    while ($node = readdir $dh) {
	next if ($node eq "..");
	next if ($node eq ".");
	next if ($node =~ /.cif/);
	next if ($node =~ /.cnf/);
	# next if ($node =~ /s_.*tif/);
	if ($node =~ /s_.*tif/) {
	    if ("$dir/$node" =~ /Images\/L\d\d\d\/[CP]\d+\.\d\/s_\d.*\.tif/) {
		next;
	    } else {
#		print "keep tif file  $dir/$node\n";
	    }
	}
	if (-d "$dir/$node") {
	    # print "$prefix work on dir $dir/$node\n";
	    $tot_byte += count_dir("$dir/$node", $indent + 1);
	    $tot_dir++;
	} else {
	    # print "$prefix work on file $dir/$node\n";
	    $tot_file ++;
	    @tmp = stat("$dir/$node");
	    $file_size = $tmp[7];
	    $tot_byte += $file_size;
	}
    }
    closedir $dh;
    $tot_K = sprintf("%d", $tot_byte / 1024);
    $tot_M = sprintf("%d", $tot_K / 1024);
    $tot_G = sprintf("%d", $tot_M / 1024);
    if ($indent < 4) {
#	print "$prefix $dir : $tot_dir dirs, $tot_file files, $tot_byte bytes, $tot_K K, $tot_M M, $tot_G G\n";
    }
    return ($tot_byte);
}

1;
