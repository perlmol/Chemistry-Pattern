#!/home/ivan/bin/perl -s

$| = 1;

use Chemistry::Smiles;
use Chemistry::Mol;
use blib;
use Chemistry::Pattern;
use strict;
use warnings;
use Data::Dumper;

our $Debug ||= 0;
our $permute ||= 0;
our $overlap ||= 1;
$Chemistry::Pattern::Atom::Debug = $Debug;
$Chemistry::Pattern::Debug = $Debug;

if (@ARGV < 2) {
    die "ptest.pl <pattern> <mol>...\n";
}

my ($patt_str, @mol_strs) = @ARGV;

print "Pattern: $patt_str\n";

my $mol_parser = new Chemistry::Smiles();
my ($mol, $patt);

$mol_parser->parse($patt_str, $patt = Chemistry::Pattern->new);
$patt->options(permute => $permute, overlap => $overlap);

for my $mol_str (@mol_strs) {
    print "Mol: $mol_str\n";
    $mol_parser->parse($mol_str, $mol = Chemistry::Mol->new);

    my @ret;
    while ($patt->match($mol) ) {
        @ret = $patt->atom_map;
        print "Matched: (@ret)\n";
    }
    print "Matched: ()\n";
}


