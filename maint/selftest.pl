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

  {name => 'material-astrosil',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool astrosil 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-gra',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-gra 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-nano',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-nano 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-org',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-org 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-p',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-p 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-c-z',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool c-z 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch3oh-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch3oh-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch3oh-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch3oh-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch4-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch4-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ch4-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ch4-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co2-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co2-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co2-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co2-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-co2-w',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool co2-w 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-cor-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool cor-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-fe-c',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool fe-c 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-fes',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool fes 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-h2o-a',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool h2o-a 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-h2o-w',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool h2o-w 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-nh3-m',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool nh3-m 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-c-mg00',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-c-mg00 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-c-mg100',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-c-mg100 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-c-mg95',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-c-mg95 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-mg40',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-mg40 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-ol-mg50',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool ol-mg50 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-c-mg96',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-c-mg96 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg100',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg100 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg40',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg40 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg50',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg50 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg60',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg60 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg70',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg70 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg80',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg80 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-pyr-mg95',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool pyr-mg95 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-sic',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool sic 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'material-sio2',
   prepare => "rm -f $dir/dustkappa.dat",
   cmd => "./optool sio2 1 -a 1 -dhs -o $dir",
   getsha => "accuracy",
   glob => "$dir/dustkappa.dat"},

  {name => 'fits-output',
   prepare => "rm -f $dir/dustkappa.fits",
   cmd => "./optool -q -na 10 -nl 30 -s -fits -o $dir",
   getsha => "shasum $dir/dustkappa.fits",
   sha => 'fb1a85451e94e2d4d8a5d05317ebde9b6fc9e48e'}
  ];

die "./optool file not present" if not -e "./optool";
die "./optool is not executable" if not -x "./optool";

# Read the hashes
while (<DATA>) {
  next if /^\s*#/;
  next if /^\s*$/;
  ($name,$hashes) = split(/:::/,$_);
  $hashes{$name} = $hashes;
}

$npassed=0;$nfailed=0;$nskipped=0;
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
        $npassed++;
      } elsif ($acc == 0) {
        $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
        $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
        $nfailed++;
      } else {
        $acc = $1+1 if $acc =~ /([0-9]+):/;
        $result = "Test $test->{name} passed with $acc significant digits\n";
        $r = sprintf("Test %-30s OK to $acc significant digits",$test->{name});
        $npassed++;          
      }
    }
  } else {
    $sha = qx($getsha);
    $sha = $1 if $sha =~ /([0-9a-fA-F]+)/;
    if ($sha eq $test->{sha}) {
      $result = "Test $test->{name} passed.\n";
      $r = sprintf "Test %-30s passed",$test->{name};
      $npassed++;
    } elsif (($name eq 'fits-output') and (`./optool -feature fits` =~/False/)) {
      $result = "Test $test->{name} skipped (fits support not implemented).\n";
      $r = sprintf "Test %-30s skipped",$test->{name};
      $nskipped++;
    } else {
      $result = "Test $test->{name} failed: resulting SHA1 $sha does not match\n";
      $r = sprintf "Test %-30s FAILED %s",$test->{name},$sha;
      $nfailed++;
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
  print "===========================================\n";
  print "     Summary of optool test results\n";
  print "===========================================\n";
  for (@results) {print "$_\n"}
  print "-------------------------------------------\n";
  print "$npassed passed; $nfailed failed; $nskipped skipped\n";
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

# HASHES set on Tue Jun 23 17:47:09 CEST 2026

noscat:::1000:9072a6ed8d::10:9072a6ed8d::9:9072a6ed8d::8:9072a6ed8d::7:9072a6ed8d::6:9072a6ed8d::5:7ae138c623::4:947fa18f7a::3:6a52eda424::2:f10735c019::1:735d6b4bfc
quick:::1000:50b04c7fb1::10:50b04c7fb1::9:50b04c7fb1::8:50b04c7fb1::7:50b04c7fb1::6:50b04c7fb1::5:450608b209::4:05e0b78365::3:d5879f62fd::2:c6b4e0f4f1::1:0e8545e812
full-diana-lowres:::1000:f713d8fee4::10:f713d8fee4::9:f713d8fee4::8:f713d8fee4::7:f713d8fee4::6:f713d8fee4::5:a915addc59::4:2fbf0da698::3:1eb364a5cb::2:470e40be74::1:4c7810ea34
sparse-file:::1000:729a13f00f::10:729a13f00f::9:729a13f00f::8:729a13f00f::7:729a13f00f::6:729a13f00f::5:0cd24fd5f3::4:e557507ba7::3:96be51faa2::2:6ac3255198::1:33192a748e
lognormal-lowres:::1000:c981d20c94::10:c981d20c94::9:c981d20c94::8:c981d20c94::7:c981d20c94::6:c981d20c94::5:9ecb420553::4:3a7838cda7::3:b52d7d6ccf::2:f6b8ee31d8::1:eef6209b43
diana-plus-ice-mantle:::1000:e36055d501::10:e36055d501::9:e36055d501::8:e36055d501::7:e36055d501::6:e36055d501::5:6e285ba351::4:1aa432db5b::3:44c979efe8::2:9fc1bc30bd::1:e2e3ab3095
div-in-10-parts:::1000:543173fbbf::10:543173fbbf::9:543173fbbf::8:543173fbbf::7:543173fbbf::6:543173fbbf::5:dd5cd73b3e::4:7f9480294c::3:70b628e1c3::2:e0ac13d9b4::1:a5e5c99b3c
large-grain:::1000:a462478131::10:a462478131::9:a462478131::8:a462478131::7:a462478131::6:a462478131::5:272ba71ce2::4:a4f01b6a74::3:560c89df5d::2:ce54641700::1:a502b7b366
high-angular-resolution:::1000:348e4d854e::10:348e4d854e::9:348e4d854e::8:348e4d854e::7:348e4d854e::6:348e4d854e::5:3d01de6a76::4:31a7db852f::3:5e94854926::2:208da2f3cc::1:5f79fa7ac0
pure-ice-grain:::1000:c828a10bbb::10:c828a10bbb::9:c828a10bbb::8:c828a10bbb::7:c828a10bbb::6:c828a10bbb::5:078488a6ae::4:39dda3e651::3:8dd2744937::2:f8b37e77b1::1:fe91fc3062
chop-peak:::1000:61cffefdfe::10:61cffefdfe::9:61cffefdfe::8:61cffefdfe::7:61cffefdfe::6:61cffefdfe::5:dcb8aa271a::4:1c2dff56fb::3:80222f3c41::2:74f3951f6c::1:daff212cdf
cde-rayleigh:::1000:503e063b60::10:503e063b60::9:503e063b60::8:503e063b60::7:503e063b60::6:503e063b60::5:de3c2965c1::4:950b7e9122::3:1a286b7ca5::2:1e44bf7ab9::1:a8f7f57541
mmf-scat:::1000:0cdd617d50::10:0cdd617d50::9:0cdd617d50::8:0cdd617d50::7:0cdd617d50::6:0cdd617d50::5:3ca2219b86::4:aa09fa0c4e::3:c6930562d4::2:f34ba2591e::1:9cb68b3ded
mmf-opac:::1000:c1637cf821::10:c1637cf821::9:c1637cf821::8:c1637cf821::7:c1637cf821::6:c1637cf821::5:2f79f583e7::4:7f738c463f::3:e03c264ce5::2:3a97e5e2a4::1:5e4701113e
material-astrosil:::1000:474f4d632f::10:474f4d632f::9:474f4d632f::8:474f4d632f::7:474f4d632f::6:474f4d632f::5:fb8bb5b544::4:157ad05256::3:4277bed275::2:267fe07096::1:2f1ce61d69
material-c-gra:::1000:f602f99aec::10:f602f99aec::9:f602f99aec::8:f602f99aec::7:f602f99aec::6:f602f99aec::5:a0aefae420::4:0a4d229c2a::3:e2eff1c0e3::2:433cee7176::1:13ba3402ff
material-c-nano:::1000:e3eec4947d::10:e3eec4947d::9:e3eec4947d::8:e3eec4947d::7:e3eec4947d::6:e3eec4947d::5:9069bbfbe5::4:6012c16267::3:9a31a94e1e::2:439725a5ee::1:7f7e00cb10
material-c-org:::1000:37c76c4162::10:37c76c4162::9:37c76c4162::8:37c76c4162::7:37c76c4162::6:37c76c4162::5:9b7873962c::4:4c0201b2cd::3:6117cd5361::2:cbc64976d1::1:1b959b90ea
material-c-p:::1000:507e1407cb::10:507e1407cb::9:507e1407cb::8:507e1407cb::7:507e1407cb::6:507e1407cb::5:47676d55c0::4:e05c6e142d::3:8bfa9e817e::2:95171f7a47::1:be27aa0858
material-c-z:::1000:834ab13b56::10:834ab13b56::9:834ab13b56::8:834ab13b56::7:834ab13b56::6:834ab13b56::5:66eb4b81e6::4:16bc401db6::3:e89c663a32::2:1f742af4b9::1:4cde488186
material-ch3oh-a:::1000:08c7bfdcaf::10:08c7bfdcaf::9:08c7bfdcaf::8:08c7bfdcaf::7:08c7bfdcaf::6:08c7bfdcaf::5:2383bb99cb::4:c4fc2ed99e::3:5eb0fb1ab2::2:61d3cb8c4b::1:18833d17fe
material-ch3oh-c:::1000:ba28387c2f::10:ba28387c2f::9:ba28387c2f::8:ba28387c2f::7:ba28387c2f::6:ba28387c2f::5:2b01e7562b::4:29d9f1e16d::3:7eb4009cfd::2:31926938c1::1:3e9317f784
material-ch4-a:::1000:2537b3d2d0::10:2537b3d2d0::9:2537b3d2d0::8:2537b3d2d0::7:2537b3d2d0::6:2537b3d2d0::5:b9f1db450f::4:bb071bc67e::3:c7add95763::2:7ffaf925ae::1:4cccad6316
material-ch4-c:::1000:fad71acaee::10:fad71acaee::9:fad71acaee::8:fad71acaee::7:fad71acaee::6:fad71acaee::5:8cd4050ad9::4:61e9b30a88::3:7b85a5af64::2:f6525bab9b::1:44e7f4bc2c
material-co-a:::1000:44c3c423ae::10:44c3c423ae::9:44c3c423ae::8:44c3c423ae::7:44c3c423ae::6:44c3c423ae::5:d72dce6f9c::4:53c336b02b::3:4e3384c5f6::2:0102deafce::1:77a4b81da9
material-co2-a:::1000:c4b2241ea5::10:c4b2241ea5::9:c4b2241ea5::8:c4b2241ea5::7:c4b2241ea5::6:c4b2241ea5::5:d030d2faf4::4:5f48cf7662::3:1dfddc8108::2:7647ba8974::1:0f46bf90c4
material-co2-c:::1000:a9aab3bd08::10:a9aab3bd08::9:a9aab3bd08::8:a9aab3bd08::7:a9aab3bd08::6:a9aab3bd08::5:f7bc05e601::4:2df4a95e2a::3:3fcca4bbf6::2:c424f447b2::1:55b87457e5
material-co2-w:::1000:20ea8af16b::10:20ea8af16b::9:20ea8af16b::8:20ea8af16b::7:20ea8af16b::6:20ea8af16b::5:f452a67cd6::4:20a600b620::3:7732fa27ec::2:9f26e7a18c::1:2feb4dc9ef
material-cor-c:::1000:ea8846f351::10:ea8846f351::9:ea8846f351::8:ea8846f351::7:ea8846f351::6:ea8846f351::5:7ccff1b0e0::4:d22420fa9f::3:9dd5e44775::2:d0c98e7715::1:e490309d5e
material-fe-c:::1000:a7e0882858::10:a7e0882858::9:a7e0882858::8:a7e0882858::7:a7e0882858::6:a7e0882858::5:9b6422af5f::4:d71a880508::3:580625a496::2:a1bb641af8::1:4f6440d9a2
material-fes:::1000:f1a6e85cd2::10:f1a6e85cd2::9:f1a6e85cd2::8:f1a6e85cd2::7:f1a6e85cd2::6:f1a6e85cd2::5:f320967f7e::4:c3f525dbd4::3:19513c7582::2:c45ad130db::1:81c9233fd2
material-h2o-a:::1000:1cc2f7d8b0::10:1cc2f7d8b0::9:1cc2f7d8b0::8:1cc2f7d8b0::7:1cc2f7d8b0::6:1cc2f7d8b0::5:1edf3a9de0::4:c33c57de74::3:10042a4608::2:78d7a21fe4::1:4e5a9d3e79
material-h2o-w:::1000:c7b3b977aa::10:c7b3b977aa::9:c7b3b977aa::8:c7b3b977aa::7:c7b3b977aa::6:c7b3b977aa::5:b82a5a69a4::4:580db46c9c::3:1d7f5fd355::2:e87ef189f7::1:6f8c6856a2
material-nh3-m:::1000:e52a3e5379::10:e52a3e5379::9:e52a3e5379::8:e52a3e5379::7:e52a3e5379::6:e52a3e5379::5:3ecb5ed601::4:a803ab189c::3:9c830e0b36::2:5bb922cef9::1:38bf3a1da3
material-ol-c-mg00:::1000:03a1440730::10:03a1440730::9:03a1440730::8:03a1440730::7:03a1440730::6:03a1440730::5:74b5f5e652::4:fb035f29da::3:5a07000d15::2:c8d798d98f::1:a6990c325f
material-ol-c-mg100:::1000:cbfa010468::10:cbfa010468::9:cbfa010468::8:cbfa010468::7:cbfa010468::6:cbfa010468::5:fcdc90bc5f::4:0eaf964fcc::3:5ef1013e6a::2:b967df9e82::1:0f8bdf12c0
material-ol-c-mg95:::1000:51872a2b0d::10:51872a2b0d::9:51872a2b0d::8:51872a2b0d::7:51872a2b0d::6:51872a2b0d::5:e8700789c9::4:e5c927adac::3:ae753ef078::2:1f1b5c4954::1:6d850ba6c7
material-ol-mg40:::1000:e100c9fd0f::10:e100c9fd0f::9:e100c9fd0f::8:e100c9fd0f::7:e100c9fd0f::6:e100c9fd0f::5:a6de8e691e::4:70c979f4aa::3:7f9b5126c0::2:f7e3c07fad::1:9a407f0070
material-ol-mg50:::1000:35981058ab::10:35981058ab::9:35981058ab::8:35981058ab::7:35981058ab::6:35981058ab::5:06785d8eee::4:6bd4933507::3:24ba7df9d1::2:b0893e5222::1:56cc06805b
material-pyr-c-mg96:::1000:9a763f6cf1::10:9a763f6cf1::9:9a763f6cf1::8:9a763f6cf1::7:9a763f6cf1::6:9a763f6cf1::5:fd0b5fe2b1::4:7279776456::3:bde8a6c41b::2:491e4f9b60::1:edd737c0c4
material-pyr-mg100:::1000:006f826aa8::10:006f826aa8::9:006f826aa8::8:006f826aa8::7:006f826aa8::6:006f826aa8::5:f298e15f28::4:d0096f5ce1::3:3464273c36::2:df76bb901c::1:2069d7332d
material-pyr-mg40:::1000:7936cf927e::10:7936cf927e::9:7936cf927e::8:7936cf927e::7:7936cf927e::6:7936cf927e::5:980f41b92a::4:952afb9689::3:c3f26eaaa5::2:a34fd5dbc4::1:9c67308a63
material-pyr-mg50:::1000:1110c58f18::10:1110c58f18::9:1110c58f18::8:1110c58f18::7:1110c58f18::6:1110c58f18::5:93464785cb::4:f74e8b4b40::3:53088551fc::2:03d47c29a7::1:d43e6ca640
material-pyr-mg60:::1000:2a779c5745::10:2a779c5745::9:2a779c5745::8:2a779c5745::7:2a779c5745::6:2a779c5745::5:de5cd9768b::4:2fa7f574ce::3:3b7c173a88::2:fce18be067::1:3ffe7f1514
material-pyr-mg70:::1000:ea759b2a2e::10:ea759b2a2e::9:ea759b2a2e::8:ea759b2a2e::7:ea759b2a2e::6:ea759b2a2e::5:8b25cd0203::4:0be6db58ef::3:bbe945cbe3::2:3465013ff7::1:8a443061f1
material-pyr-mg80:::1000:68c16e5dae::10:68c16e5dae::9:68c16e5dae::8:68c16e5dae::7:68c16e5dae::6:68c16e5dae::5:7c6f0907e4::4:9e7a1a6f64::3:53c6202ccc::2:18cc2da261::1:f7852660aa
material-pyr-mg95:::1000:49cff86e1d::10:49cff86e1d::9:49cff86e1d::8:49cff86e1d::7:49cff86e1d::6:49cff86e1d::5:97333a64dd::4:9608bb8be7::3:09052695d2::2:d17cf5066b::1:02f9219d27
material-sic:::1000:6fdd610f3a::10:6fdd610f3a::9:6fdd610f3a::8:6fdd610f3a::7:6fdd610f3a::6:6fdd610f3a::5:b96721b850::4:b8260feae3::3:8cfe457909::2:fcc28b58a1::1:f9380c7083
material-sio2:::1000:3562056501::10:3562056501::9:3562056501::8:3562056501::7:3562056501::6:3562056501::5:58c6efc29d::4:9ccfb5679f::3:195356b4ce::2:fd4f6e986e::1:3d3b4152e2

