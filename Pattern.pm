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

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    $self->reset;
    $self;
}


=item $mol->new_atom(name => value, ...)

Shorthand for $mol->add_atom(Chemistry::Atom->new(name => value, ...));

=cut

sub new_atom {
    my $self = shift;
    $self->add_atom(Chemistry::Pattern::Atom->new(@_));
}


=item $mol->new_bond(name => value, ...)

Shorthand for $mol->add_bond(Chemistry::Bond->new(name => value, ...));

=cut

sub new_bond {
    my $self = shift;
    $self->add_bond(Chemistry::Pattern::Bond->new(@_));
}

our $Debug = 0;

sub reset {
    my ($self, @mols) = @_;
    print "Resetting to (@mols)\n" if $Debug;
    $self->{pending_mols} = [@mols];
    $self->{pending_atoms} = [];
    $self->{already_matched} = {};
}

sub already_matched {
    my ($self, $key) = @_;
    if ($self->{already_matched}{$key}) {
        print "already matched $key\n" if $Debug;
        return 1;
    } else {
        $self->{already_matched}{$key} = 1;
        print "first match of $key\n" if $Debug;
        return 0;
    }
}

sub next_atom {
    my $self = shift;
    my $atom;
    print "next_atom\n" if $Debug;
    if (@{$self->{pending_atoms}}) {
        $atom = shift @{$self->{pending_atoms}};
        print "\tatom $atom\n" if $Debug;
    } elsif (@{$self->{pending_mols}}) {
        my $mol = shift @{$self->{pending_mols}};
        print "\tmol $mol\n" if $Debug;
        $self->map_to($mol);
        $self->{pending_atoms} = [$mol->atoms];
        $atom = shift @{$self->{pending_atoms}};
        print "\tatom $atom\n" if $Debug;
    }
    $atom;
}

Chemistry::Obj::accessor "map_to";

sub match {
    my ($self, $mol) = @_;
    print "match $self $mol\n" if $Debug;
    if (defined($mol) and $self->map_to != $mol) {
        $self->reset($mol);
    }
    my $match = $self->match_next;
    print "returning match: '$match'\n" if $Debug;
    $match;
}

sub match_next {
    my $self = shift;
    my $match;
    print "match_next\n" if $Debug;
    while (my $atom = $self->next_atom) {
        ($match) = $self->match_local($atom);
        if ($match) {
            my $match_key = join " ", sort map {$_->id} $self->atom_map;
            if ($self->already_matched($match_key)) {
                $match = 0;
                next;
            } else {
                last;
            }
        }
    }
    $match;
}

sub match_local {
    my ($patt, $atom) = @_;
    print "local_match: $atom with $patt\n" if $Debug;
    for ($patt->atoms, $patt->bonds) { $_->map_to(undef) }
    $patt->atoms(1)->match(
        where => $atom,
        from_where => [],
        from_what => [],
    );
}

sub atom_map {
    my $self = shift;
    map { $_->map_to } $self->atoms(@_);
}

sub bond_map {
    my $self = shift;
    map { $_->map_to } $self->bonds(@_);
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

