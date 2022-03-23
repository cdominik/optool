#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

use Digest::SHA qw/ sha1_hex /;
use Getopt::Std;
getopts('s');
$dir = 'selftest_optool';

my $tests = [
  {name => 'noscat',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'full-diana-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'lognormal-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -a 0.001 4.9 0.1:1.0 30 -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'diana-plus-ice-mantle',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -na 20 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'div-in-10-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -s -d 3 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat_*.dat"},

  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 10 -na 20 -s 720 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'chop-peak',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -chop 4 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'cde-rayleigh',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -a 0.01 0.1 3.5 10 -l 10 60 30 -s -cde -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'mmf-scat',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool pyr 1 -q -a 10 30 -mmf 0.3 .005 -l 30 3000 20 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'mmf-opac',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr 1 -q -a 10 30 3.5 2 -mmf 0.3 2.2 -nlam 3 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => 'f9c0abb6bd112b08ad1e25edff010173fa8b3369'}
  ];

# Read the hashes
while (<DATA>) {
  next if /^\s*#/;
  next if /^\s*$/;
  ($name,$hashes) = split(/:::/,$_);
  $hashes{$name} = $hashes;
}

foreach $test (@$tests) {
  $name = $test->{name};
  #next unless $test->{name} eq "noscat";
  print "\nRunning test: $name....\n";
  system($test->{prepare}) if $test->{prepare};
  system($test->{cmd});
  $getsha = $test->{getsha};
  if ($getsha eq "accuracy") {
    $glob = $test->{glob};
    @files = glob($glob);
    #print "files are: @files\n"; 
    if ($opt_s) {
      $hashes = &accuracy_hashes(@files);
      $allhashes .= "${name}:::${hashes}\n";
    } else {
      $acc = &find_accuracy($hashes{$name},@files);
      #print "accuracy found is $acc\n";
      if ($acc == 1000) {
        $result = "Test $test->{name} perfect pass\n";
        $r = sprintf "Test %-30s passed",$test->{name};
      } elsif ($acc == 0) {
        $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
        $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
      } else {
        $acc = $1+1 if $acc =~ /([0-9]+):/;
        $result = "Test $test->{name} passed with $acc significant digits\n";
        $r = sprintf("Test %-30s OK to $acc significant digits",$test->{name})
      }
    }
  } else {
    $sha = qx($getsha);
    $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
    if ($sha eq $test->{sha}) {
      $result = "Test $test->{name} passed.\n";
      $r = sprintf "Test %-30s passed",$test->{name};
    } else {
      $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
      $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
    }
  }
  print "$result";
  push @results,$r;
}
if ($opt_s) {
  $date = `date`;
  open SELF,">> $0" or die "Unable to open SELF\n";
  print SELF "\n\n# HASHES set on $date\n",$allhashes,"\n";
} else {
  print "\n";
  print "==========================================\n";
  print "     Summary of optool test results\n";
  print "==========================================\n";
  for (@results) {print "$_\n"}
}

sub accuracy_hashes {
  # Create a series of hashes of the concatenation of all
  # files in @_, where the accuracy of number with exponential format
  # is limited to N digits after the comma.
  # The return value is a string with 10 different hashes,
  # accuracy 10 to 1. The hashes have the format ACC:HASH where ACC
  # is an integer and HASH is the first 10 digits of the sha1 hash of
  # the concatenated files with whitespace compressed and accuracy limited.
  # The different hashes are joined with "::" into a single string.
  my @files = sort @_;
  my $all = "";
  my @all = ();
  my @numbers = ();
  my $file,$content,$fh;
  my $nh = 10;            # Number of chars from sha1 hash to use
  my $debug = 0;
  my $maxacc = 10;        # Maximum accuracy for which to produce a hash
  my $perfect = 1000;     # "accuracy: key for perfect, i.e. unmodified
  my $scatlim = 1e-2;
  
  foreach $file (@files) {
    # print "Doing file $file\n";
    open $fh, '<', $file or die "Can't open file $!";
    $content = do { local $/; <$fh> };

    #$content =~ s/[-+][0-9]\.[0-9]+[eE]-[0-9]+/1.234567E-01/g;
    $content =~ s/[\n\t ]+/ /g; # turn multiple whitespace into single space
    $content =~ s/\s+$//;
    $content =~ s/^\s+//;
    print "c3",$content if $debug;
    push @all,$content;
  }

  $all = join("\n",@all);
  print "all is ",$all if $debug;

  @digits = (1..$maxacc);
  push @digits,$perfect;
  @hashes = ();
  for $i (reverse @digits) {
    $re = "([-+]?[0-9]+[.][0-9]{$i})([0-9]+)([eE][-+][0-9]+)";
    #print "re is $re\n";
    $all =~ s/$re/"$1" . "0" x length($2) . "$3"/ge if $i<$perfect;
    #print $i,":::\n",$all,"\n\n";
    push @hashes,"$i:" . substr(&sha1_hex($all),0,$nh);
  }
  return join("::",@hashes);
}

sub find_accuracy {
  # Find the highest accuracy for which the files are identical
  # to the files used to produce the hashes given
  # in the first argument
  #print "@_";
  my $goalhashes = shift @_;
  my @files = sort(@_);
  my $hashes = &accuracy_hashes(@files);
  $hashes = "::" . $hashes . "::";
  @goalhashes = split('::',$goalhashes);
  for $goal (@goalhashes) {
    return $goal if $hashes =~ m/::${goal}::/;
  }
  return 0;
}


__DATA__

# HASHES set on Wed Mar 23 09:04:46 CET 2022

noscat:::1000:da718bde67::10:da718bde67::9:da718bde67::8:da718bde67::7:da718bde67::6:da718bde67::5:3f4021068e::4:e98093c1d5::3:5523f2812f::2:787172fae1::1:275000a906
quick:::1000:5edb7ab23d::10:5edb7ab23d::9:5edb7ab23d::8:5edb7ab23d::7:5edb7ab23d::6:5edb7ab23d::5:cac77430ea::4:762a325380::3:176d4b54ec::2:d8187069a7::1:a8814aeede
full-diana-lowres:::1000:3ea7cb24f1::10:3ea7cb24f1::9:3ea7cb24f1::8:3ea7cb24f1::7:3ea7cb24f1::6:3ea7cb24f1::5:a246d4a9b3::4:5bc5c1ac58::3:e23f6da565::2:608bbbc793::1:ec7788cb46
lognormal-lowres:::1000:46587e5e53::10:46587e5e53::9:46587e5e53::8:46587e5e53::7:46587e5e53::6:46587e5e53::5:ff91b21175::4:7e60d04fb6::3:4d806f3e8e::2:552ee2ebed::1:401b4da754
diana-plus-ice-mantle:::1000:1dbd6b4e4a::10:1dbd6b4e4a::9:1dbd6b4e4a::8:1dbd6b4e4a::7:1dbd6b4e4a::6:1dbd6b4e4a::5:eff6e0a5e4::4:b65f5e9337::3:88699024f0::2:307539603c::1:6a8f4c4b8c
div-in-10-parts:::1000:15b1bde278::10:15b1bde278::9:15b1bde278::8:15b1bde278::7:15b1bde278::6:15b1bde278::5:8ff3a42df3::4:848b47786a::3:9df87d1355::2:2989369c2d::1:97788937cf
high-angular-resolution:::1000:e84201ec9b::10:e84201ec9b::9:e84201ec9b::8:e84201ec9b::7:e84201ec9b::6:e84201ec9b::5:1e88de0e26::4:fdcfc911f8::3:18e9869eda::2:ef09826527::1:8266231a22
pure-ice-grain:::1000:b06b419a5a::10:b06b419a5a::9:b06b419a5a::8:b06b419a5a::7:b06b419a5a::6:b06b419a5a::5:e970744d9e::4:586e7113d9::3:47f5d924fc::2:38d778f953::1:4446b78e80
chop-peak:::1000:3f8926838a::10:3f8926838a::9:3f8926838a::8:3f8926838a::7:3f8926838a::6:3f8926838a::5:2813723c66::4:9af50f70a8::3:4d1b3e1eef::2:dd23a11051::1:dc7e38cf18
cde-rayleigh:::1000:1461d7c082::10:1461d7c082::9:1461d7c082::8:1461d7c082::7:1461d7c082::6:1461d7c082::5:b2354b40c7::4:a3fedd3b86::3:3328804a35::2:1e44bf7ab9::1:a8f7f57541
mmf-scat:::1000:cd3b6a51b9::10:cd3b6a51b9::9:cd3b6a51b9::8:cd3b6a51b9::7:cd3b6a51b9::6:cd3b6a51b9::5:cec618a88a::4:0cb36b405f::3:d00ffc2c92::2:b933b47486::1:67532246c8
mmf-opac:::1000:ec28210a34::10:ec28210a34::9:ec28210a34::8:ec28210a34::7:ec28210a34::6:ec28210a34::5:56d89c3302::4:1a341c612e::3:405768cc85::2:3a97e5e2a4::1:5e4701113e

