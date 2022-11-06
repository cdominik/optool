#!/usr/bin/perl

# Run a few standard cases and compare the results to stored hashes.
# When run with -s, the hashes at the end of the file are set.  You may
# want to clean out the old velues by hand.

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

  {name => 'sparse-file',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -s -nl 30 -na 20 -sparse 2.2 -sp 300 -o $dir",
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
   cmd  => "./optool -q -na 10 -nl 30 -d 3 -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa_*.dat"},

  {name => 'large-grain',
   prepare => "rm -f $dir/dustkapscatmat.dat",
   cmd  => "./optool -q -c pyr-mg70 0.696 -c c-z 0.104 -m h2o-w 0.2 -p 0.25 -a 1000 -s -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkapscatmat.dat"},

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
   sha => 'fb1a85451e94e2d4d8a5d05317ebde9b6fc9e48e'}
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

# HASHES set on Sun Nov  6 21:34:14 CET 2022

noscat:::1000:9072a6ed8d::10:9072a6ed8d::9:9072a6ed8d::8:9072a6ed8d::7:9072a6ed8d::6:9072a6ed8d::5:7ae138c623::4:947fa18f7a::3:6a52eda424::2:f10735c019::1:735d6b4bfc
quick:::1000:e36ece4c73::10:e36ece4c73::9:e36ece4c73::8:e36ece4c73::7:e36ece4c73::6:e36ece4c73::5:090fb54eb2::4:42e115764d::3:3e5ebda2cb::2:f110bf4b12::1:322fc08e00
full-diana-lowres:::1000:9ea6b6a191::10:9ea6b6a191::9:9ea6b6a191::8:9ea6b6a191::7:9ea6b6a191::6:9ea6b6a191::5:b0eb48e95f::4:203feabcdb::3:5762bc1600::2:7503525287::1:e8859fe950
sparse-file:::1000:d3317d994c::10:d3317d994c::9:d3317d994c::8:d3317d994c::7:d3317d994c::6:d3317d994c::5:cfcae31444::4:58b5420597::3:36712d26d3::2:4ac5fd0dc6::1:6eb3b3d969
lognormal-lowres:::1000:0657684c37::10:0657684c37::9:0657684c37::8:0657684c37::7:0657684c37::6:0657684c37::5:d95a419365::4:e5839a581e::3:1ebfb49e26::2:feb2f63812::1:62ce857d5f
diana-plus-ice-mantle:::1000:6cf1c4986f::10:6cf1c4986f::9:6cf1c4986f::8:6cf1c4986f::7:6cf1c4986f::6:6cf1c4986f::5:22c6728c33::4:be16f59d9d::3:f7ee4adc9b::2:17eb6c74b4::1:ad3825a6d9
div-in-10-parts:::1000:6a70384b2c::10:6a70384b2c::9:6a70384b2c::8:6a70384b2c::7:6a70384b2c::6:6a70384b2c::5:69f9ac1a41::4:1278f6d160::3:70b628e1c3::2:e0ac13d9b4::1:a5e5c99b3c
large-grain:::1000:a763df5351::10:a763df5351::9:a763df5351::8:a763df5351::7:a763df5351::6:a763df5351::5:4bc1c32b68::4:7ace124dfe::3:45bfb56578::2:b6493bc275::1:ddc07bb95e
high-angular-resolution:::1000:e605872e2c::10:e605872e2c::9:e605872e2c::8:e605872e2c::7:e605872e2c::6:e605872e2c::5:dd19f975a8::4:1e66e9685d::3:3ecf68f980::2:946272ac50::1:51a81b6a09
pure-ice-grain:::1000:1f1dc3247e::10:1f1dc3247e::9:1f1dc3247e::8:1f1dc3247e::7:1f1dc3247e::6:1f1dc3247e::5:4cd021606f::4:2a2008a9cd::3:bf7aa7396f::2:81724ae0fd::1:4446b78e80
chop-peak:::1000:161eba279f::10:161eba279f::9:161eba279f::8:161eba279f::7:161eba279f::6:161eba279f::5:721d991fcb::4:f379c06917::3:66afd01d52::2:8e2703236a::1:14162288d1
cde-rayleigh:::1000:503e063b60::10:503e063b60::9:503e063b60::8:503e063b60::7:503e063b60::6:503e063b60::5:de3c2965c1::4:950b7e9122::3:1a286b7ca5::2:1e44bf7ab9::1:a8f7f57541
mmf-scat:::1000:a919f27210::10:a919f27210::9:a919f27210::8:a919f27210::7:a919f27210::6:a919f27210::5:cacdd11446::4:aa09fa0c4e::3:c6930562d4::2:f34ba2591e::1:9cb68b3ded
mmf-opac:::1000:c1637cf821::10:c1637cf821::9:c1637cf821::8:c1637cf821::7:c1637cf821::6:c1637cf821::5:2f79f583e7::4:7f738c463f::3:e03c264ce5::2:3a97e5e2a4::1:5e4701113e

