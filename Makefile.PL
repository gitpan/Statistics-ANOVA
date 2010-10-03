use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Statistics::ANOVA',
    AUTHOR              => 'Roderick Garton <rgarton@cpan.org>',
    VERSION_FROM        => 'lib/Statistics/ANOVA.pm',
    ABSTRACT_FROM       => 'lib/Statistics/ANOVA.pm',
    PL_FILES            => {},
    PREREQ_PM => {
        'Algorithm::Combinatorics' => '0.23',
	'Math::Cephes' => '0.43',
        'Scalar::Util' => '1.18',
	'Statistics::Lite' => '3.2',
        'Statistics::DependantTTest' => '0.03',
	'Statistics::Descriptive' => '2.6',         
	'Statistics::TTest' => '1.1',
	'Test::More' => '0.0',
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Statistics-ANOVA-0.07-*' },
);
