use Test::More;
use Chemistry::Smiles;
use Chemistry::Mol;
use strict;

my @files;

BEGIN { 
    @files = glob "t/*.pat";
    plan tests => 1 + @files;

    use_ok('Chemistry::Pattern');
};


my $mol_parser = new Chemistry::Smiles();

for my $file (@files) {
    open F, $file or die "couldn't open $file\n";   
    my ($mol_str, $patt_str, $expected_match) = map { /: ([^\n\r]*)/g } <F>;
    
    my ($mol, $patt);
    Chemistry::Atom->reset_id;
    $mol_parser->parse($mol_str, $mol = Chemistry::Mol->new);
    $mol_parser->parse($patt_str, $patt = Chemistry::Pattern->new);

    my @ret = $patt->match_first($mol);
    my $match = "(@ret)";
    is ($match, $expected_match, "$mol_str =~ $patt_str");
}
