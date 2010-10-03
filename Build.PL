use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::ANOVA',
    license             => 'perl',
    dist_author         => 'Roderick Garton <rgarton@cpan.org>',
    dist_version_from   => 'lib/Statistics/ANOVA.pm',
    build_requires => {
      'Algorithm::Combinatorics' => '0.23',
      'Math::Cephes' => '0.43',
      'Scalar::Util' => '1.18',
      'Statistics::Lite' => '3.2',
      'Statistics::DependantTTest' => '0.03',
      'Statistics::Descriptive' => '2.6',         
      'Statistics::TTest' => '1.1',
      'Test::More' => '0.0',
    },
    add_to_cleanup      => [ 'Statistics-ANOVA-0.07-*' ],
);

$builder->create_build_script();
