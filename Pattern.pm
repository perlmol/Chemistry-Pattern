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
    $self->{options} = {overlap=>1, permute=>0};
    $self;
}

sub options {
    my $self = shift;
    if (@_ == 1) {
        $self->{options} = {%{$self->{options}}, %{$_[0]}};
    } else {
        $self->{options} = {%{$self->{options}}, @_};
    }
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
    $self->{current_atom} = $self->next_atom;
    $self->match_local_init($self->{current_atom});
}

sub already_matched {
    my ($self, @atoms) = @_;
    my @ids = map {$_->id} @atoms;
    my $unsorted_key = join " ", @ids;
    my $key;
    if ($self->{options}{permute}) {
        $key = $unsorted_key;
    } else {
        $key = join " ", sort @ids;
    }

    if ($self->{already_matched}{$key}) {
        print "already matched $unsorted_key\n" if $Debug;
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

sub atom_map {
    my $self = shift;
    map { $_->map_to } $self->atoms(@_);
}

sub bond_map {
    my $self = shift;
    map { $_->map_to } $self->bonds(@_);
}


# testing a different implementation
sub match {
    my ($self, $mol) = @_;
    print "match $self $mol\n" if $Debug;
    if (defined($mol) and $self->map_to ne $mol) { #TODO
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

    $self->match_local_init($self->{current_atom}) unless $self->{options}{overlap};
    while (1) {
        ($match) = $self->match_local_next;
        if ($match) {
            if ($self->already_matched($self->atom_map)) {
                $match = 0;
                next; # already matched this; try again
            } else {
                unless ($self->{options}{overlap}) {
                    
                    $self->next_atom;
                }
                last; # matched!
            }
        } else { 
            my $atom = $self->next_atom;
            last unless $atom;
            $self->match_local_init($atom);
        }
    }
    $match;
}


sub match_local_init {
    my ($patt, $atom) = @_;
    my $state = 1;
    my @stack = {
        from_where_bond_i => [0], 
        from_what_bond_i => [0],
        from_where => [$atom], 
        from_what => [$patt->atoms(1)],
    };
    my $mol = $patt->map_to;

    print "match_local_init(",$patt->atoms(1),", $atom)\n" if $Debug;
    for ($patt->atoms, $patt->bonds) { $_->map_to(undef) }
    if ($mol and $patt->{options}{overlap}) {
        for ($mol->atoms, $mol->bonds) { $_->del_attr('painted') }
    }
    $patt->{stack} = \@stack;
    $patt->{state} = $state;

}


sub match_local_next {
    my ($patt) = @_;
    my $match = 0;
    my $state = \$patt->{state};
    my $stack = $patt->{stack};

    print "match_local_next\n" if $Debug;

    LOOP: while (1) {
        # initialize variables for this iteration
        my $from_where = $stack->[-1]{from_where};
        my $from_what = $stack->[-1]{from_what};
        my $from_where_bond_i = $stack->[-1]{from_where_bond_i};
        my $from_what_bond_i = $stack->[-1]{from_what_bond_i};
        my $where = $from_where->[-1];
        my $what  = $from_what->[-1];
        my $where_bond_i = $from_where_bond_i->[-1];
        my $what_bond_i = $from_what_bond_i->[-1];

        last unless $what;
        print "    $$state: $what($what_bond_i),$where($where_bond_i)\n" if $Debug >= 2;
        dump_stack($stack) if $Debug >= 3;

        if ($$state == 1) {
            if ($what->map_to) { # ring closure in pattern
                if ($where eq $what->map_to) { # ring also closed ok in mol
                    print "\tring closed at where: $where; what: $what; map: ",
                        $what->map_to, "\n" if $Debug; 
                    #continue to check bonds
                    $$state = 2;
                    next
                } else {
                    print "\tring didn't close at where: $where; what: $what; map: ",
                        $what->map_to, "\n" if $Debug; 
                }
            } elsif ($where->attr("painted")) { # ring closure in mol
                print "\tatom $where already visited\n" if $Debug;
            } elsif ($what->test($where)) { ### ATOM TEST
                print "\tatom $where matches $what\n" if $Debug;
                $where->attr("painted", 1);
                $what->map_to($where);
                $stack->[-1]{paint} = 1;
                # Now check bonds
                $$state = 2;
                next; #XXX
            } 
            print "\tatom $where does not match $what\n" if $Debug;
            # backtrack
            $$state = 3;
            next; #XXX
        } elsif ($$state == 2) {
            ### start of match_bonds
            my ($patt_bond) = ($what->bonds)[$what_bond_i];
            if (!$patt_bond) { # no more bonds to match?
                print "\tNo more bonds to match at $what\n" if $Debug;
                if (@$from_where > 1) { # go back and finish previous atom
                    push @$stack, {
                        from_where_bond_i => 
                        [@$from_where_bond_i[0 .. $#$from_where_bond_i - 2],0], 
                        from_what_bond_i => 
                        [@$from_what_bond_i[0 .. $#$from_what_bond_i - 1]],
                        from_where => 
                        [@$from_where[0 .. $#$from_where - 1]], 
                        from_what => 
                        [@$from_what[0 .. $#$from_what - 1]],
                    };
                    next; #XXX
                } else {
                    $match = 1; # Finally matched! This is the deepest point
                    print "\tFinally matched!\n" if $Debug;
                    $$state = 3; 
                    last;
                } 
            } else { # Match next bond
                ++$from_what_bond_i->[-1], next if $patt_bond->map_to; #XXX
                my $mol_bond = ($where->bonds)[$where_bond_i];
                if (!$mol_bond) {
                    # no more bonds left to try; backtrack
                    print "\tno more bonds left to try at $where; backtracking\n" if $Debug;
                    $$state = 3;
                    next; #XXX
                } else {
                    ++$from_where_bond_i->[-1], next if $mol_bond->attr("painted"); #XXX
                    if ($patt_bond->type eq $mol_bond->type) { ### BOND TEST
                        print "\tbond $mol_bond matches $patt_bond\n" if $Debug;

                        # now check the atom on the other side. First, get atom
                        my ($patt_nei) = grep {$_ ne $what} $patt_bond->atoms;
                        my ($mol_nei)  = grep {$_ ne $where} $mol_bond->atoms;
                        print "\tChecking neighbor $mol_nei with $patt_nei\n" if $Debug;

                        # recursive call to match atom
                        $patt_bond->map_to($mol_bond);
                        $mol_bond->attr("painted", 1);
                        $where->attr("painted", 1);
                        $what->map_to($where);

                        push @$stack, {
                            from_where => [@$from_where, $mol_nei], 
                            from_what => [@$from_what, $patt_nei],
                            from_where_bond_i => [@$from_where_bond_i, 0], 
                            from_what_bond_i => [@$from_what_bond_i, 0],
                        };
                        $$state = 1;
                        next; #XXX
                    } else {
                        print "\tbond $mol_bond does not match $patt_bond\n" if $Debug;
                        ++$from_where_bond_i->[-1], next; # try next bond
                    }
                }
            }
        } elsif ($$state == 3) { # Backtracking mode
            my $fwhatb = 0;
            do {
                last LOOP unless @$stack > 1;
                print "\tpopping stack\n" if $Debug;
                pop @$stack;
                my $sf = $stack->[-1];
                if ($sf->{paint}) {
                    $where->del_attr("painted");
                    $what->map_to(undef);
                }
                ($from_what, $from_where, $from_what_bond_i, 
                $from_where_bond_i) = @$sf{qw(from_what from_where 
                    from_what_bond_i from_where_bond_i)};
                $where = $from_where->[-1];
                $what  = $from_what->[-1];
                ($from_where->[-1]->bonds)[$from_where_bond_i->[-1]]->del_attr("painted");
                $fwhatb = ($from_what->[-1]->bonds)[$from_what_bond_i->[-1]];
            } until ($fwhatb);
            $fwhatb->map_to(undef);
            ++$from_where_bond_i->[-1];
            $$state = 2;
            next;
        } else {
            die "wrong state $$state";
        }
    } 
    print "\treturning from match_local: '$match'\n" if $Debug;
    $match;
}

sub dump_stack {
    my ($stack) = @_;
    my $i = 0;
    for my $sf (@$stack) {
        print "\t", $i++, "\n";
        for my $key (sort keys %$sf) {
            print "\t\t$key: @{$sf->{$key}}\n";
        }
    }
}

1;

=back

=head1 BUGS

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

