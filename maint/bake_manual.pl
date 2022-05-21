#!/usr/bin/perl

# Program to extract the text from UserGuide.org and bake it into
# a fortran routine so that the manual and parts of it can be displayed
# by optool.

$cmd = "/Applications/Emacs.app/Contents/MacOS/Emacs UserGuide.org --batch -f org-ascii-export-to-ascii --kill";
system($cmd);

$manual = do {
  local $/ = undef;
  open my $fh, "<", "UserGuide.txt" or die "could not open UserGuide.txt: $!";
  <$fh>;
};

@lines = split(/^/,$manual);
$manf90 = "subroutine manual(what)\n  character*(*) what\n  write(*,*)\n  if (what.eq.'all') then\n";
$inopt = 0;
foreach (@lines) {
  last if /^Footnotes/;
  if ($skip and /^$/) {
    $skip = 0;
    next;
  }
  if ($skip) {
    next;
  }
  if (/^[[]/) {
    $skip = 1;
    next;
  }
  if (/^ +`\[?-([a-zA-Z]+)/) {
    # Start of an option
    $inopt = 1;
    if ($1 eq "a") {
      $manf90 .= "  endif\n  if((what.eq.'-a').or.(what.eq.'-amin').or.(what.eq.'-amax').or.(what.eq.'-apow').or.(what.eq.'-na').or.(what.eq.'all')&\n          .or.(what.eq.'-amean').or.(what.eq.'-asig')) then\n";
    } elsif ($1 eq "l") {
      $manf90 .= "  endif\n  if((what.eq.'-l').or.(what.eq.'-lmin').or.(what.eq.'-lmax').or.(what.eq.'-nl').or.(what.eq.'-nlam').or.(what.eq.'all')) then\n";
    } else {
      $manf90 .= "  endif\n  if((what.eq.'-$1').or.(what.eq.'all')) then\n";
    }
  }
  if ($inopt and /^\S/) {
    $inopt = 0;
    $manf90 .= "  endif\n  if(what.eq.'all') then\n";
  }
  if (length($_) >= 110 and /\\quad/) {
    # This is the size distribution table. Fix it for this purpose.
    $_=~ s/\\left//g;
    $_=~ s/\\right//g;
    $_=~ s/\{\\rm m\}/m/g;
    $_=~ s/\\sigma\b/sig/g;
    $_=~ s/\\exp/exp\b/g;
    $_=~ s/\\ln\b/ln/g;
    $_=~ s/\\propto\b/ ~/g;
    $_=~ s/\$//g;
    $_=~ s/\\quad\b/  /g;
    $_=~ s/\\frac\{1\}\{2\}/0.5/g;
    $_=~ s/|//g;
    $_=~ s/[ \t]+$//;
    $_=~ s/powerlaw     /powerlaw/;
  }
    elsif (length($_) >= 110) {
    # This is the material table, do special stuff to make it look ok.
    $a = substr($_,2,58);
    $b = substr($_,98,14);
    $b =~ s/\]\s*$//;
    $_ = "$a $b";
    s/\*/ /g;
    s/{//g;
    s/}/  /g;
    s/\\fbox/     /g;
  }
  s/\\hfill\{\}/     /g;
  chomp;
  s/'/ /g;
  s/`/ /g;
  next if ($last eq "") and ($_ eq "");
  $last = $_;
  $manf90 .= "     write(*,'(\"$_\")')\n";
}
$manf90 .= "  endif\nend subroutine manual\n";

print "$manf90";
