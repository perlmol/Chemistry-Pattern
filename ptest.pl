#!/home/ivan/bin/perl

use Chemistry::Smiles;
use Chemistry::Mol;
use blib;
use Chemistry::Pattern;
use strict;
use warnings;
use Data::Dumper;


my $mol_str = $ARGV[0] || "CCCC";
my $patt_str = $ARGV[1] || "CC";

print "Mol: $mol_str\nPattern: $patt_str\n";

my $mol_parser = new Chemistry::Smiles();
my ($mol, $patt);

$mol_parser->parse($mol_str, $mol = Chemistry::Mol->new);
$mol_parser->parse($patt_str, $patt = Chemistry::Pattern->new);

my @ret;

@ret = $patt->match_first($mol);
print "Matched (@ret)\n";
