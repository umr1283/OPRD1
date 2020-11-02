package Projet;

use File::Basename qw(dirname);
use Cwd qw(abs_path);

### Chargement des fichiers de configuration du projet
sub load_config {
	my(@file);
	my(%config);
	my(%hash);
	
	my $racine = dirname(abs_path($0));
	@file = listFile("$racine/Config");
	for $file (@file){
		%hash = file2hash("$racine/Config/$file");
		@config{ keys %hash } = values %hash;
	}
	return %config;
}

### Chargement des fichiers de configuration du run
sub load_run {
	my($file) = @_;
	my(%hash);
	
	%hash = file2hash("$file");
	
	return %hash;
}

### Envoi du mail de fin d'analyse
sub notification {

    use MIME::Lite;
	my($email,$run) = @_;
	
      # Creation d'un message avec plusieurs parties (multipart)
     my $msg = new MIME::Lite
                 From    =>"<XXX\@XXX>",
                 To      =>"$email",
                 Subject =>"Analyse de $run",
                 Type    =>'multipart/mixed';
     # Ajout des pieces attachees
     # Chaque "attach" a les memes arguments que "new":
     attach $msg
                 Type =>'TEXT',
                 Data =>"Bonjour, \n\nL'analyse de $run est terminee.\n\nHave a good day \nRobot BIBS \n";

     $msg->send('smtp','localhost',Debug=>0);


} 

### Vérification du fichier de configuration du run
sub check_config {
	my(%run) = @_;
	my($erreur);
	
	if($run{name} eq ""){$erreur = $erreur."Merci de renseigner un nom de run/analyse dans le fichier de configuration\n";}
	if($run{email} eq ""){$erreur = $erreur."Merci de renseigner une adresse email dans le fichier de configuration\n";}
	
	return $erreur;
}

############################################################################
############################################################################
############################################################################

### Liste les dossiers d'un dossier
sub listDir {
    my($dir) = @_;
    my($fhdir) = 'FHDIR';
    my(@liFi);
    my($fic);

    opendir($fhdir, $dir);
    while ($fic = readdir($fhdir)) {
       next if ( $fic eq "." );
       next if ( $fic eq ".." );
       next unless -d "$dir/$fic";
       push @liFi, $fic;
    }
    closedir($fhdir);
    return @liFi;
}

### Liste les fichiers d'un dossier
sub listFile {
    my($dir) = @_;
    my($fhdir) = 'FHDIR';
    my(@liFi);
    my($fic);

    opendir($fhdir, $dir);
    while ($fic = readdir($fhdir)) {
       next if ( $fic eq "." );
       next if ( $fic eq ".." );
       next unless -f "$dir/$fic";
       push @liFi, $fic;
    }
    closedir($fhdir);
    return @liFi;
}

### Charge un fichier en table de hash
sub file2hash {
    my($file) = @_;
    my($in, $line, $key, $val);
    my(%hash);

    $in = file2var("$file");
    $in = "" unless (defined($in));
    return 0 unless ($in ne "");
    for $line (split('\n', $in)) {
	next if ($line =~ /^\s*$/);
	next if ($line =~ /^\#/);
    #\015 retour chariot
    #\012 nouvelle ligne
	$line =~ s/\015//g;
	$line =~ s/\012//g;
	my (@val);
	($key, @val) = split(/\t/, $line);
	$val = join(" ", @val);
	$hash{$key} = $val;
    }
    return %hash;
}

### Charge un fichier en variable
sub file2var {
    my($file) = @_;
    my($save, $in);
    my($fh1) = 'FHf2v';

    $in = "";
    unless (-f $file) {
	return $in;
    }
    $save = $/;
    undef $/;
    open($fh1, $file) || Die("$0 $NAME file2var : cannot open $file $!");
    $in = <$fh1>;
    close($fh1);
    $/ = $save;
    return $in;
}


1;
