package Chemistry::Pattern::Bond;
$VERSION = '0.01';

=head1 NAME

Chemistry::Pattern::Bond

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

use 5.006;
use strict;
use base qw(Chemistry::Bond);

=head1 METHODS

=over 4

=cut


Chemistry::Obj::accessor('map_to');
Chemistry::Obj::accessor('test_sub');

sub test {
    my ($what, $where) = @_;
    if ($what->test_sub) {
         return $what->test_sub->($what, $where);
    } else {
         return $what->order eq $where->order;
    }
}


1;

=back

=head1 SEE ALSO

L<Chemistry::Mol>, L<Chemistry::Atom>, L<Chemistry::Tutorial>

=head1 AUTHOR

Ivan Tubert E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2003 Ivan Tubert. All rights reserved. This program is free
software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

