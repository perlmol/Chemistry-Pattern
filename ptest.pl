#!/home/ivan/bin/perl -s

$| = 1;

use Chemistry::File::SMILES;
use Chemistry::Mol;
use blib;
use Chemistry::Pattern;
use strict;
use warnings;
#use diagnostics;
use Data::Dumper;
use Benchmark;

our $debug ||= 0;
our $permute ||= 0;
our $overlap ||= 1;
our $anchor;
my %options = (permute => $permute, overlap => $overlap);

#$Chemistry::Pattern::Atom::DEBUG = $debug;
$Chemistry::Pattern::DEBUG = $debug;

if (@ARGV < 2) {
    die "ptest.pl <pattern> <mol>...\n";
}

my ($patt_str, @mol_strs) = @ARGV;

print "Pattern: $patt_str\n";
print "Options: ", join(" ", %options), "\n";
print "Anchor: ", $anchor || '', "\n";

my $patt = Chemistry::Pattern->parse($patt_str, format => 'smiles');
$patt->options(%options);
#$patt->atoms(1)->test_sub(sub{
    #my ($what, $where) = @_;
    #$where->bonds == 3 ? 1 : 0;
#});
    
for my $mol_str (@mol_strs) {
    print "Mol: $mol_str\n";
    my $mol = Chemistry::Mol->parse($mol_str, format => 'smiles');

    my @ret;
    my $atom;
    $atom = $mol->atoms($anchor) if $anchor;
    while ($patt->match($mol, atom => $atom) ) {
        @ret = $patt->atom_map;
        print "Matched: (@ret)\n";
    }
    print "Matched: ()\n";
}

#timethis(100, \&do_test);

sub do_test {
    for my $mol_str (@mol_strs) {
        my $mol = Chemistry::Mol->parse($mol_str, format => 'smiles');
        my @ret;
        while ($patt->match($mol) ) {
            @ret = $patt->atom_map;
        }
    }
}

