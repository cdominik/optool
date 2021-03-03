#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

use Digest::SHA qw/ sha1_hex /;
use Getopt::Std;
getopts('s');
#$opt_s=1;
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

  {name => 'mmf-scat',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool pyr 1 -q -a 10 30 -mmf 0.3 .005 -l 30 3000 20 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

  {name => 'mmf-opac',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr 1 -q -a 10 30 5 -mmf 0.3 2.2 -nlam 10 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => 'f3825e408341b92082b29b4b019d9ba6791448f4'}
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

# HASHES set on Wed Mar  3 12:27:10 CET 2021

noscat:::1000:bbf283849a::10:bbf283849a::9:bbf283849a::8:bbf283849a::7:bbf283849a::6:bbf283849a::5:c088b8ac5d::4:a2a0d382b9::3:7b77f4ba62::2:f67461f634::1:a112d33dca
quick:::1000:f97c88a0c0::10:f97c88a0c0::9:f97c88a0c0::8:f97c88a0c0::7:f97c88a0c0::6:f97c88a0c0::5:539494a6c9::4:b78d4983f3::3:e49b9a0b64::2:9cae49f013::1:3f15349353
full-diana-lowres:::1000:e867a1b403::10:e867a1b403::9:e867a1b403::8:e867a1b403::7:e867a1b403::6:e867a1b403::5:4f3cec98f0::4:a65f2874c8::3:d866f3221c::2:dd78657084::1:8d6542b969
diana-plus-ice-mantle:::1000:21eebd6077::10:21eebd6077::9:21eebd6077::8:21eebd6077::7:21eebd6077::6:21eebd6077::5:839d046b3b::4:c64bc410f7::3:ff0166a58e::2:8b62600eb6::1:dabc0fa993
div-in-10-parts:::1000:0b8a45279a::10:0b8a45279a::9:0b8a45279a::8:0b8a45279a::7:0b8a45279a::6:0b8a45279a::5:e495379b89::4:30ad35b1d9::3:7538eb7514::2:e062eea233::1:a8ca6f3d02
high-angular-resolution:::1000:42711eaa45::10:42711eaa45::9:42711eaa45::8:42711eaa45::7:42711eaa45::6:42711eaa45::5:023cee39e1::4:9538a4ddb9::3:fce67b9dad::2:25fb7e8f7b::1:34fa47f12e
pure-ice-grain:::1000:2383df5e97::10:2383df5e97::9:2383df5e97::8:2383df5e97::7:2383df5e97::6:2383df5e97::5:b1281e49ae::4:076d1b07dd::3:dcad6d9bf1::2:1d4854116a::1:41ab50f7fe
chop-peak:::1000:52010915c7::10:52010915c7::9:52010915c7::8:52010915c7::7:52010915c7::6:52010915c7::5:24815e8a58::4:824b81d09a::3:aa2bec3994::2:51ad3b64a7::1:7c2cc0b31e
mmf-scat:::1000:a7b6a7b0a0::10:a7b6a7b0a0::9:a7b6a7b0a0::8:a7b6a7b0a0::7:a7b6a7b0a0::6:a7b6a7b0a0::5:af89e82476::4:5972544c62::3:c4066996eb::2:113bc75125::1:aea99d7a70
mmf-opac:::1000:8b50e20bb4::10:8b50e20bb4::9:8b50e20bb4::8:8b50e20bb4::7:8b50e20bb4::6:8b50e20bb4::5:56c5ec03cf::4:880714b108::3:c6a4133664::2:a2f68fbf9e::1:8c8429759a

