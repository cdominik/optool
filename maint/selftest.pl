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
    } elsif (($name eq 'fits-output') and (`./optool -feature fits` =~/False/)) {
      $result = "Test $test->{name} skipped (fits support not implemented).\n";
      $r = sprintf "Test %-30s skipped",$test->{name};
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

# HASHES set on Thu Apr 13 10:01:22 CEST 2023

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
material-astrosil:::1000:c48c0a31c3::10:c48c0a31c3::9:c48c0a31c3::8:c48c0a31c3::7:c48c0a31c3::6:c48c0a31c3::5:cb378303e5::4:51886437ed::3:4277bed275::2:267fe07096::1:2f1ce61d69
material-c-gra:::1000:8990ebb988::10:8990ebb988::9:8990ebb988::8:8990ebb988::7:8990ebb988::6:8990ebb988::5:a0aefae420::4:0a4d229c2a::3:e2eff1c0e3::2:433cee7176::1:13ba3402ff
material-c-nano:::1000:33937fda39::10:33937fda39::9:33937fda39::8:33937fda39::7:33937fda39::6:33937fda39::5:26116efc2d::4:7604f8cf0c::3:d17b2eaf2a::2:439725a5ee::1:7f7e00cb10
material-c-org:::1000:8850afa3cf::10:8850afa3cf::9:8850afa3cf::8:8850afa3cf::7:8850afa3cf::6:8850afa3cf::5:bcde7a5466::4:c3f43059f7::3:45ba61352d::2:cbc64976d1::1:1b959b90ea
material-c-p:::1000:7846f004c0::10:7846f004c0::9:7846f004c0::8:7846f004c0::7:7846f004c0::6:7846f004c0::5:7216a03873::4:eeb3b64b14::3:8bfa9e817e::2:95171f7a47::1:be27aa0858
material-c-z:::1000:7dbbba6f5c::10:7dbbba6f5c::9:7dbbba6f5c::8:7dbbba6f5c::7:7dbbba6f5c::6:7dbbba6f5c::5:ce0098773e::4:f69b6262c1::3:e89c663a32::2:1f742af4b9::1:4cde488186
material-ch3oh-a:::1000:fe04f0ff91::10:fe04f0ff91::9:fe04f0ff91::8:fe04f0ff91::7:fe04f0ff91::6:fe04f0ff91::5:3f94478ff6::4:13132cb958::3:5eb0fb1ab2::2:61d3cb8c4b::1:18833d17fe
material-ch3oh-c:::1000:8158ad8375::10:8158ad8375::9:8158ad8375::8:8158ad8375::7:8158ad8375::6:8158ad8375::5:e11684b2c0::4:92bcfc45f2::3:7eb4009cfd::2:31926938c1::1:3e9317f784
material-ch4-a:::1000:8971109d2e::10:8971109d2e::9:8971109d2e::8:8971109d2e::7:8971109d2e::6:8971109d2e::5:2380c89767::4:97bb0e6d48::3:c7add95763::2:7ffaf925ae::1:4cccad6316
material-ch4-c:::1000:64568d10b5::10:64568d10b5::9:64568d10b5::8:64568d10b5::7:64568d10b5::6:64568d10b5::5:542a84dc0e::4:e45576f6b0::3:7b85a5af64::2:f6525bab9b::1:44e7f4bc2c
material-co-a:::1000:f19a99eaee::10:f19a99eaee::9:f19a99eaee::8:f19a99eaee::7:f19a99eaee::6:f19a99eaee::5:a4e95dca90::4:e47407a208::3:d006a2023e::2:0102deafce::1:77a4b81da9
material-co2-a:::1000:35e40564d5::10:35e40564d5::9:35e40564d5::8:35e40564d5::7:35e40564d5::6:35e40564d5::5:0d0397b504::4:5f48cf7662::3:1dfddc8108::2:7647ba8974::1:0f46bf90c4
material-co2-c:::1000:e8626b1132::10:e8626b1132::9:e8626b1132::8:e8626b1132::7:e8626b1132::6:e8626b1132::5:ab6d7afa97::4:b8bde6bbbd::3:107a2b6e41::2:c424f447b2::1:55b87457e5
material-co2-w:::1000:2e56fbda5f::10:2e56fbda5f::9:2e56fbda5f::8:2e56fbda5f::7:2e56fbda5f::6:2e56fbda5f::5:dceb3aaa79::4:1ae0a93902::3:7732fa27ec::2:9f26e7a18c::1:2feb4dc9ef
material-cor-c:::1000:8c2682bc02::10:8c2682bc02::9:8c2682bc02::8:8c2682bc02::7:8c2682bc02::6:8c2682bc02::5:0f88e1e833::4:8677cc6556::3:9dd5e44775::2:d0c98e7715::1:e490309d5e
material-fe-c:::1000:a7e0882858::10:a7e0882858::9:a7e0882858::8:a7e0882858::7:a7e0882858::6:a7e0882858::5:9b6422af5f::4:d71a880508::3:580625a496::2:a1bb641af8::1:4f6440d9a2
material-fes:::1000:90503cad36::10:90503cad36::9:90503cad36::8:90503cad36::7:90503cad36::6:90503cad36::5:48e460aa6c::4:d083f32f2c::3:19513c7582::2:c45ad130db::1:81c9233fd2
material-h2o-a:::1000:7adf8807fc::10:7adf8807fc::9:7adf8807fc::8:7adf8807fc::7:7adf8807fc::6:7adf8807fc::5:08468d8ffa::4:30cffe9348::3:10042a4608::2:78d7a21fe4::1:4e5a9d3e79
material-h2o-w:::1000:a33d3201f2::10:a33d3201f2::9:a33d3201f2::8:a33d3201f2::7:a33d3201f2::6:a33d3201f2::5:6a87ebf4cc::4:056e6dc68b::3:1d7f5fd355::2:e87ef189f7::1:6f8c6856a2
material-nh3-m:::1000:bf6bec4533::10:bf6bec4533::9:bf6bec4533::8:bf6bec4533::7:bf6bec4533::6:bf6bec4533::5:81062860b5::4:a803ab189c::3:9c830e0b36::2:5bb922cef9::1:38bf3a1da3
material-ol-c-mg00:::1000:0ea06deeca::10:0ea06deeca::9:0ea06deeca::8:0ea06deeca::7:0ea06deeca::6:0ea06deeca::5:c199b499f6::4:c4820b9343::3:5a07000d15::2:c8d798d98f::1:a6990c325f
material-ol-c-mg100:::1000:dca5e53109::10:dca5e53109::9:dca5e53109::8:dca5e53109::7:dca5e53109::6:dca5e53109::5:59f6ac7ae7::4:c3b93a0997::3:207b04d511::2:b967df9e82::1:0f8bdf12c0
material-ol-c-mg95:::1000:fdda20b18f::10:fdda20b18f::9:fdda20b18f::8:fdda20b18f::7:fdda20b18f::6:fdda20b18f::5:5e6fde0da3::4:e5c927adac::3:ae753ef078::2:1f1b5c4954::1:6d850ba6c7
material-ol-mg40:::1000:ce6315e95c::10:ce6315e95c::9:ce6315e95c::8:ce6315e95c::7:ce6315e95c::6:ce6315e95c::5:54a7f67600::4:c2617f12bb::3:7f9b5126c0::2:f7e3c07fad::1:9a407f0070
material-ol-mg50:::1000:4fe368f4f9::10:4fe368f4f9::9:4fe368f4f9::8:4fe368f4f9::7:4fe368f4f9::6:4fe368f4f9::5:80d323f106::4:d328e0cc38::3:24ba7df9d1::2:b0893e5222::1:56cc06805b
material-pyr-c-mg96:::1000:c192d05c4f::10:c192d05c4f::9:c192d05c4f::8:c192d05c4f::7:c192d05c4f::6:c192d05c4f::5:023df5fd60::4:7279776456::3:bde8a6c41b::2:491e4f9b60::1:edd737c0c4
material-pyr-mg100:::1000:601a282cba::10:601a282cba::9:601a282cba::8:601a282cba::7:601a282cba::6:601a282cba::5:e98759d499::4:9067da4dc0::3:3464273c36::2:df76bb901c::1:2069d7332d
material-pyr-mg40:::1000:b50ac5d2b1::10:b50ac5d2b1::9:b50ac5d2b1::8:b50ac5d2b1::7:b50ac5d2b1::6:b50ac5d2b1::5:e101c928fa::4:776717c35d::3:c3f26eaaa5::2:a34fd5dbc4::1:9c67308a63
material-pyr-mg50:::1000:8238849f54::10:8238849f54::9:8238849f54::8:8238849f54::7:8238849f54::6:8238849f54::5:011299422f::4:97d73b228e::3:c87bc600a2::2:03d47c29a7::1:d43e6ca640
material-pyr-mg60:::1000:801bdfe057::10:801bdfe057::9:801bdfe057::8:801bdfe057::7:801bdfe057::6:801bdfe057::5:7062cc7077::4:59a1ae7fd0::3:3b7c173a88::2:fce18be067::1:3ffe7f1514
material-pyr-mg70:::1000:b4f5d88b03::10:b4f5d88b03::9:b4f5d88b03::8:b4f5d88b03::7:b4f5d88b03::6:b4f5d88b03::5:a5f2096db4::4:30715dae62::3:bbe945cbe3::2:3465013ff7::1:8a443061f1
material-pyr-mg80:::1000:f437114526::10:f437114526::9:f437114526::8:f437114526::7:f437114526::6:f437114526::5:00721f0bbb::4:a8a070b9e8::3:53c6202ccc::2:18cc2da261::1:f7852660aa
material-pyr-mg95:::1000:2bd718c8c9::10:2bd718c8c9::9:2bd718c8c9::8:2bd718c8c9::7:2bd718c8c9::6:2bd718c8c9::5:4023b1a32c::4:2f579f72c9::3:09052695d2::2:d17cf5066b::1:02f9219d27
material-sic:::1000:e439e7183c::10:e439e7183c::9:e439e7183c::8:e439e7183c::7:e439e7183c::6:e439e7183c::5:664cf5a36b::4:4dba745736::3:8cfe457909::2:fcc28b58a1::1:f9380c7083
material-sio2:::1000:0ae08d69d6::10:0ae08d69d6::9:0ae08d69d6::8:0ae08d69d6::7:0ae08d69d6::6:0ae08d69d6::5:e16a20184d::4:150c554188::3:8b3ce1b132::2:fd4f6e986e::1:3d3b4152e2

