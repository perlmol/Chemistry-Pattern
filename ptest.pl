#!/home/ivan/bin/perl -s

$| = 1;

use Chemistry::File::SMILES;
use Chemistry::Mol;
use blib;
use Chemistry::Pattern;
use strict;
use warnings;
use Data::Dumper;

our $Debug ||= 0;
our $permute ||= 0;
our $overlap ||= 1;
my %options = (permute => $permute, overlap => $overlap);

$Chemistry::Pattern::Atom::Debug = $Debug;
$Chemistry::Pattern::Debug = $Debug;

if (@ARGV < 2) {
    die "ptest.pl <pattern> <mol>...\n";
}

my ($patt_str, @mol_strs) = @ARGV;

print "Pattern: $patt_str\n";
print "Options: ", join(" ", %options), "\n";

my $patt = Chemistry::Pattern->parse($patt_str, format => 'smiles');
$patt->options(%options);
$patt->atoms(1)->test_sub(sub{
    my ($what, $where) = @_;
    $where->bonds == 3 ? 1 : 0;
});
    

for my $mol_str (@mol_strs) {
    print "Mol: $mol_str\n";
    my $mol = Chemistry::Mol->parse($mol_str, format => 'smiles');


    my @ret;
    while ($patt->match($mol) ) {
        @ret = $patt->atom_map;
        print "Matched: (@ret)\n";
    }
    print "Matched: ()\n";
}


