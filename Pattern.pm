package Chemistry::Pattern;
$VERSION = '0.01';

=head1 NAME

Chemistry::Pattern - Substructure pattern matching

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

use 5.006;
use strict;
use Carp;
use base qw(Chemistry::Mol);
use Chemistry::Pattern::Atom;
use Chemistry::Pattern::Bond;

=head1 METHODS

=over 4

=item Chemistry::Pattern->new(name => value, ...)

=cut

# so far, inherit and do nothing...
sub new {
    my $class = shift;
    $class->SUPER::new(@_);
}


=item $mol->new_atom(name => value, ...)

Shorthand for $mol->add_atom(Chemistry::Atom->new(name => value, ...));
It has the disadvantage that it doesn't let you create a subclass of 
Chemistry::Atom.

=cut

sub new_atom {
    my $self = shift;
    $self->add_atom(Chemistry::Pattern::Atom->new(@_));
}


=item $mol->new_bond(name => value, ...)

Shorthand for $mol->add_bond(Chemistry::Bond->new(name => value, ...));
It has the disadvantage that it doesn't let you create a subclass of 
Chemistry::Atom.

=cut

sub new_bond {
    my $self = shift;
    $self->add_bond(Chemistry::Pattern::Bond->new(@_));
}

our $Debug = 0;

sub match_all {
    my ($patt, $mol) = @_;
    print "first_match: $mol with $patt\n" if $Debug;
    my @ret;
    for my $atom ($mol->atoms) {
        my @match = $patt->match_local($atom);
        push @ret, \@match if @match;
    }
    @ret;
}

sub match_first {
    my ($patt, $mol) = @_;
    print "first_match: $mol with $patt\n" if $Debug;
    my @ret;
    for my $atom ($mol->atoms) {
        @ret = $patt->match_local($atom);
        return @ret if @ret;
    }
}

sub match_local {
    my ($patt, $atom) = @_;
    print "local_match: $atom with $patt\n" if $Debug;
    my ($match, @ret) = $patt->atoms(1)->match(
        where => $atom,
        #what => $patt->atoms(1),
        from_where => '',
        from_what => '',
    );
    return @ret;
}

1;

=back

=head1 BUGS

Blatant memory leaks. Due to the use of circular references, Perl's current
garbage collector never cleans up molecule, atom, and bond objects. A future
version should address this.

=head1 SEE ALSO

L<Chemistry::Atom>, L<Chemistry::Bond>, L<Chemistry::File>,
L<Chemistry::Tutorial>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2003 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

