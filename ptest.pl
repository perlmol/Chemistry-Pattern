#!/home/ivan/bin/perl -s

use Chemistry::Smiles;
use Chemistry::Mol;
use blib;
use Chemistry::Pattern;
use strict;
use warnings;
use Data::Dumper;

our $Debug ||= 0;
$Chemistry::Pattern::Atom::Debug = $Debug;

my $mol_str = $ARGV[0] || "CCCC";
my $patt_str = $ARGV[1] || "CC";

print "Mol: $mol_str\nPattern: $patt_str\n";

my $mol_parser = new Chemistry::Smiles();
my ($mol, $patt);

$mol_parser->parse($mol_str, $mol = Chemistry::Mol->new);
$mol_parser->parse($patt_str, $patt = Chemistry::Pattern->new);

my @ret;

#$patt->atoms(1)->test_sub(sub{1});
@ret = $patt->match_first($mol);
print "Matched: (@ret)\n";

