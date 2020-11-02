package Folder;

sub check_run {
	my($run) = @_;
	my($warning);

	if(-d $run){$warning = $warning."Le run existe déjà sur le disque, il s'agit d'une ré-analyse\n";}

	return $warning;
}

sub create_run {
	my($run) = @_;

	mkpath("$run");
	#system("chmod -R 775 $run");

	return 1;
}

### Création d'un chemin
sub mkpath {
    my($path) = @_;
    my($dir);
    my($tmp) = "";

    $tmp = ($path =~ /^\.\//) ? "." : "";
    for $dir (split('/', $path)) {
	next if ($dir eq "");
	$tmp .= "/$dir";
	checkfolder($tmp);
    }
}
sub checkfolder {
    my($dir) = @_;

    unless ( -d $dir ) { system ("mkdir -m 777 $dir"); }
}

1;
