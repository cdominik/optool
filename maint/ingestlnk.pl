#!/usr/bin/perl

# Take all the link files on the command line and turn them into
# subroutines of a fortaran 90 code spit out to STDOUT. The names of
# the lnk files must follow narrow conventions - see the User guide
# for more information.

# Usage: ingestlnk path/to/lnk/file.lnk [...]

unless (@ARGV) {
  &usage();
  exit(1);
}

# ------------------------------------------------------------------
# Format the helptext from lnk_data/lnk-help.txt
# also, extract the association from quick keys to keys
# from that text
# ------------------------------------------------------------------

$helpblock = &makehelpblock();
while ($helpblock =~ /(\S+)\s*->\s*(\S+)/g) {
  $quick = $1;
  $key = $2;
  $key_to_quick{$key} = $quick;
  print STDERR "Establishing quick key by looking at lnk-help.txt: $quick -> $key\n";
}

# ------------------------------------------------------------------
# Get the template for ref_ind.f90 into the code variable
# ------------------------------------------------------------------
open my $fh, '<', 'maint/ref_ind.template' or die "Can't open file $!";
$code = do { local $/; <$fh> };
close $fh;

# ------------------------------------------------------------------
# Loop over the files
# ------------------------------------------------------------------
push @allkeys,"ice-w"; # hardcoded for backward compatibility
while ($file = shift) {
  # Extract the material key and the Reference
  if ($file =~ /(.*\/)?(.+)-(\w+[0-9]+).lnk$/) {
    $key = $2;
    $ref = $3;
    die "Material key $key must start with a letter" unless $key=~/^[a-zA-Z]/;
    die "Duplicate material key: $key" if $seen_raw_key{$key}++;
    # Set the name of the subroutine
    $sbrname = $key;
    push @allkeys,$sbrname;
    $sbrname =~ s/[^a-zA-Z0-9]/_/g;
    push @allkeys,$sbrname;
    die "Duplicate FORTRAN name $sbrname results from key $key\n   This means that there are two material keys in your file list\n   that result into the same FORTRAN name after turning all forbidden\n   characters into an underscore.\n"
      if $seen_processed_key{$sbrname}++;
  } else {
    die "File $file does not have a name like key-reference2017.lnk\n";
  }

  # open the file
  open($fh, "<" ,$file) or die "No such file $file\n";
  $line_number = 0;
  # Extract the comment lines
  $comments = "";
  while (<$fh>) {
    $line_number++;
    last unless /^\s*[#!*]/;
    next if  /^\s*[#] First line is N_lam/; # Dont need this line in the comments.
    s/^\s*[#!*]/  !/;
    $comments .= $_;
  }

  # Check the first line, which should be number of wavelengths and the density only
  ($nlam,$rho,$dummy) = (split(' ',$_));
  if ($dummy) {
    # There is a third value on the first non-comment line
    die "File seems to be missing the npoints rho line: $file\n";
  }
  if ($rho <= 0.) {
    die "rho=$rho in $file, must be larger than 0\n";
  }

  # Read the data, line by line
  @ll = (); @nn = (); @kk = ();
  foreach $i (1..$nlam) {
    $line_number++;
    $line = <$fh>;
    ($l,$n,$k) = (split(' ',$line));
    if ($n<=0) {
      print STDERR "WARNING: non-positive value of n at lambda=$l in file $file\n";
    }
    if ($k<=0) {
      print STDERR "WARNING: non-positive value of k at lambda=$l in file $file\n";
    }
    if (not ($l and $n and $k)) {
      $line =~ s/[ \t]+\n?$//;
      die "Error in line $line_number in file $file\nThe line is: $line\nMore lnk data was expected.";
    }
    push @ll,$l; push @nn,$n; push @kk,$k;
  }
  die "Extra data after expected last line $line_number in file $file\n" if ($_=<$fh> and /[^\s]/);

  # Make sure lambda is increasing in the data
  if ($ll[0] > $ll[-1]) {
    @ll = reverse(@ll);
    @nn = reverse(@nn);
    @kk = reverse(@kk);
  }
  
  # Is there q quick key for this dataset?
  $quick = $key_to_quick{$key};
  if ($quick and $quick ne $key) {
    push @allkeys,$quick;
    $quick = ",'" . $quick . "'";
  } else {
    $quick = "";
  }

  # Add to the externals and cases lists
  $externals .= "external $sbrname\n  ";
  $cases     .= "case('$key'$quick)\n     call ${sbrname}(x,y1,y2,n0,rho)\n  ";

  # Template for the subroutine
  $sbr = <<"SUBROUTINE_TEMPLATE";

subroutine ${sbrname}(l_lnk,n_lnk,k_lnk,nlam,rho)
  ! Created by ingesting $file
COMMENTS
  IMPLICIT NONE
  integer nlam,j
  real*8 l_lnk(*),n_lnk(*),k_lnk(*)
  real*8 l0($nlam),n0($nlam),k0($nlam)
  real*8 rho
  nlam = $nlam
  rho  = $rho
  data (l0(j),j=1,$nlam) /LAMBODY
  data (n0(j),j=1,$nlam) /NBODY
  data (k0(j),j=1,$nlam) /KBODY
  l_lnk(1:nlam) = l0(1:nlam)
  n_lnk(1:nlam) = n0(1:nlam)
  k_lnk(1:nlam) = k0(1:nlam)
  return
end subroutine $sbrname

SUBROUTINE_TEMPLATE

  # Insert comments and data into the template
  $sbr =~ s/COMMENTS/$comments/;
  $sbr =~ s/LAMBODY/&makebody(\@ll)/e;
  $sbr =~ s/NBODY/&makebody(\@nn)/e;
  $sbr =~ s/KBODY/&makebody(\@kk)/e;

  # Add the new subroutine to the code string
  $code .= $sbr;
}

# Add the externals, case statements and help text
$externals = substr($externals,0,-3);
$cases    .= "case('ice-w') ! backward compatibility\n     call h2o_w(x,y1,y2,n0,rho)\n  ";
$cases     = substr($cases,0,-3);

$code =~ s/EXTERNAL_DEFS/$externals/;
$code =~ s/CASE_DEFS/$cases/;
$code =~ s/HELPBLOCK/&makehelpblock()/e;
$code =~ s/KEYCHECKBLOCK/&makekeycheckblock/e;
print $code;

exit(0);

sub usage {
  print <<'EOF'
Usage: ./ingestlnk.pl path/to/lnkfile.dat [...]
This should be called from that parent directory of the lnk-data directory
The best way to call it is using "make ingest"
EOF
}

sub makebody {
  # Turn an array of numbers into something that fits into a data staement
  my @a = @{$_[0]};
  my $n_per_line = 5; # How many numbers per line?
  my $str = "";
  my $line = "  ";
  while (@a) {
    $_ = shift(@a);
    if (length($line)+1+length($_)<70) {
      $line = $line . $_ . ",";
    } else {
      $str .= " &\n$line";
      $line = "  $_,";
    }
  }
  $str = $str . " &\n" . substr($line,0,-1) . "/";
  return $str;
}

sub makehelpblock {
  my $helpblock;
  open my $fh, '<', 'lnk_data/lnk-help.txt' or die "Can't open file $!";
  while (<$fh>) {
    chomp $_;
    next if /^#/;
    $helpblock .= "  write(stdo,'(\"$_\")')\n" unless /^\s*$/;
  }
  close $fh;
  return $helpblock;
}

sub makekeycheckblock {
  my @unique;
  my %seen;
  my $str = "";
  foreach my $value (@allkeys) {
    if (! $seen{$value}) {
      push @unique, $value;
      $seen{$value} = 1;
    }
  }
  foreach my $key (@unique) {
    $str .= "  else if (trim(name) .eq. \"$key\") then\n    is_key_or_file=.true.\n";
  }
  chomp($str);
  return $str;
}
