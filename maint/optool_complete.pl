#!/usr/bin/perl


use Text::ParseWords ();

$line = $ENV{COMP_LINE};
$point = $ENV{COMP_POINT};
$before = substr($line,0,$point);
$after = substr($line,$point);

print STDERR "in here, at point $point\n";
print STDERR "in here, in line: $line\n";
print STDERR "in here, before: $before\n";
print STDERR "in here, after: $after\n";
my @list = Text::ParseWords::shellwords($line);
print STDERR join(' &&& ',@list),"\n";

print "uck\n";
print "oam\n";

