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
Chemistry::Obj::accessor('map_to');

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

