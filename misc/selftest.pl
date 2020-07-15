#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.

use Digest::SHA qw/ sha1_hex /;
use Getopt::Std;
getopts('s');
#$opt_s=1;
$dir = 'selftest_optool';

my $tests = [
  {name => 'quick',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat",
   sha => '1000:08c1e5e93e::10:08c1e5e93e::9:08c1e5e93e::8:08c1e5e93e::7:08c1e5e93e::6:08c1e5e93e::5:90f3452d60::4:e4ca4b1256::3:08fc12c7b9::2:d4712e9d2a::1:a448e211f8'},

  {name => 'full-diana-lowres',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat",
   sha => '1000:da886965eb::10:da886965eb::9:da886965eb::8:da886965eb::7:da886965eb::6:da886965eb::5:46e523a4b5::4:5bfdada8c5::3:189ccc4ada::2:26b3ca4e4e::1:f39a6c5ce8'},

  {name => 'diana-plus-ice-mantle',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.87 -c c-z 0.13 -m ice-w 0.2 -p 0.25 -s -nl 30 -na 20 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat",
   sha => '1000:4767575447::10:4767575447::9:4767575447::8:4767575447::7:4767575447::6:4767575447::5:38c6b004b3::4:c08c36c0f0::3:ce8b9ccc27::2:cbbee6c5b3::1:222eba0515'},

  {name => 'div-in-10-parts',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -na 10 -nl 30 -s -d 3 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat_*.dat",
   sha => '1000:166fb6e527::10:166fb6e527::9:166fb6e527::8:166fb6e527::7:166fb6e527::6:166fb6e527::5:06535b616b::4:3f33db8f72::3:2050818458::2:c7a33f2081::1:260fbf285d'},

  {name => 'high-angular-resolution',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -nl 10 -na 20 -s 720 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat",
   sha => '1000:08e101eb54::10:08e101eb54::9:08e101eb54::8:08e101eb54::7:08e101eb54::6:08e101eb54::5:039187374a::4:bfe7c3c3d8::3:a7f234083f::2:bb7b8f7f57::1:4f405c379d'},

  {name => 'pure-ice-grain',
   prepare => "rm -f $dir/dust*",
   cmd  => "./optool -q -c ice-w 1.0 -a 1 3 2.5 15 -l 10 100 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat",
   sha => '1000:37c201f217::10:37c201f217::9:37c201f217::8:37c201f217::7:37c201f217::6:37c201f217::5:7d23baf6cf::4:538f62eb8b::3:406e91c719::2:c92727fd2c::1:8f7afd10b2'},

  {name => 'chop-peak',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd => "./optool -q -na 10 -nl 30 -s -chop 4 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat",
   sha => '1000:1cdd2dcfc3::10:1cdd2dcfc3::9:1cdd2dcfc3::8:1cdd2dcfc3::7:1cdd2dcfc3::6:1cdd2dcfc3::5:d050dcf4ca::4:08325a421f::3:8464dd6cb9::2:051759e079::1:0a1be351e8'},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => '5db1e1afef60ffdfd94b9e4cb7a9626790d7a501'}
  ];

foreach $test (@$tests) {
  #next unless $test->{name} eq debug;
  print "\nRunning test: $test->{name}....\n";
  system($test->{prepare}) if $test->{prepare};
  system($test->{cmd});
  $getsha = $test->{getsha};
  if ($getsha eq "accuracy") {
    $glob = $test->{glob};
    @files = glob($glob);
    #print "files are: @files\n"; 
    if ($opt_s) {
      $hashes = &accuracy_hashes(@files);
      print "$test->{name} :: $hashes\n";
    } else {
      $acc = &find_accuracy($test->{sha},@files);
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
print "\n";
print "==========================================\n";
print "     Summary of optool test results\n";
print "==========================================\n";
for (@results) {print "$_\n"}

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

  foreach $file (@files) {
    # print "Doing file $file\n";
    open $fh, '<', $file or die "Can't open file $!";
    $content = do { local $/; <$fh> };
    print "c1",$content if $debug;
    $content =~ s/^\s*$//mg; $content =~ s/\n+/\n/mg; # remove empty lines
    print "c2",$content if $debug;
    $content =~ s/\s\s+/ /g; # turn multiple whitespace into single space
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
