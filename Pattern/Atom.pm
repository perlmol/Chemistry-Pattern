package Chemistry::Pattern::Atom;
$VERSION = '0.01';

=head1 NAME

Chemistry::Pattern::Atom

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Atom Attributes

=cut

use 5.006;
use strict;
use Carp;
use base qw(Chemistry::Atom);

=head1 METHODS

=over 4

=cut

our $Debug = 0;

sub test {
    my ($what, $where) = @_;
    if ($what->test_sub) {
         return $what->test_sub->($what, $where);
    } else {
         return $what->symbol eq $where->symbol;
    }
}

Chemistry::Obj::accessor('test_sub');

sub match {
    my $what = shift; # self - the pattern atom to match
    my (%opts) = @_;
    my ($where, $from_where, $from_what,) = 
        (@opts{qw(where from_where from_what )});
    my @ret;
    my $match = 0;
    print "match_atom -- where: $where; what: $what; ",
        "from_where: [@$from_where]; from_what: [@$from_what];\n" if $Debug;
    
    if ($what->attr("painted")) { # ring closure in pattern
        if ($where->attr("painted")) { # ring also closed ok in mol
            print "\tring closed at $where\n" if $Debug; 
            ($match, @ret) = $what->match_bonds(@_); # continue with same params
        } else {
            print "\tring didn't close at $where\n" if $Debug;
        }
    } elsif ($where->attr("painted")) { # ring closure in mol
        print "\tatom $where already visited\n" if $Debug;
    } elsif ($what->test($where)) { ### ATOM TEST
        print "\tatom $where matches $what\n" if $Debug;

        # Now check bonds
        $where->attr("painted", 1);
        $what->attr("painted", 1);
        ($match, @ret) = $what->match_bonds(@_); # continue with same params
        $what->attr("painted", 0);
        $where->attr("painted", 0);

        unshift @ret, $where if $match;
    } else {
        print "\tatom $where does not match $what\n" if $Debug;
    }
    return ($match, @ret);
}

sub match_bonds {
    my $what = shift; # self - the pattern atom to match
    my (%opts) = @_;
    my ($where, $from_where, $from_what) = 
        (@opts{qw(where from_where from_what)});
    my $match = 0;
    my @ret;

    print "match_bonds -- where: $where; what: $what; ",
        "from_where: [@$from_where]; from_what: [@$from_what];\n" if $Debug;

    my ($patt_bond) = grep {! $_->attr("painted")} $what->bonds;
    if (!$patt_bond) { # no more bonds to match?
        print "\tNo more bonds to match at $what\n" if $Debug;
        if (@$from_where) { # go back and finish previous atom
            #pop @$from
            $where = pop @$from_where;
            $what = pop @$from_what;
            ($match, @ret) = $what->match_bonds(
                where=>$where, 
                from_where=>$from_where, from_what=>$from_what
            );
        } else {
            $match = 1; # Finally matched! This is the deepest point
        }
    } else { # Match next bond
        for my $mol_bond (grep {!$_->attr("painted")} $where->bonds) {
            if ($patt_bond->type eq $mol_bond->type) { ### BOND TEST
                print "\tbond $mol_bond matches $patt_bond\n" if $Debug;

                # now check the atom on the other side. First, get atom
                my ($patt_nei) = grep {$_ != $what} $patt_bond->atoms;
                my ($mol_nei)  = grep {$_ != $where} $mol_bond->atoms;
                print "\tChecking neighbor $mol_nei with $patt_nei\n" if $Debug;

                # recursive call to match atom
                $patt_bond->attr("painted", 1);
                $mol_bond->attr("painted", 1);
                push @$from_where, $where;
                push @$from_what, $what;
                ($match, @ret) = $patt_nei->match(
                    where=>$mol_nei, 
                    from_where=>$from_where, from_what=>$from_what
                );
                $patt_bond->attr("painted", 0);
                $mol_bond->attr("painted", 0);

                last if $match;
            } else {
                print "\tbond $mol_bond does not match $patt_bond\n" if $Debug;
            }
        }
    }
    return ($match, @ret);
}
1;

=back

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Bond>, 
L<Math::VectorReal>, L<Chemistry::Tutorial>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2003 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

