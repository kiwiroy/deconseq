#!/usr/bin/perl

#===============================================================================
#   Author: Robert SCHMIEDER, Computational Science Research Center @ SDSU, CA
#
#   File: splitFasta.pl
#   Date: 2013-01-19
#   Version: 0.1
#
#   Usage:
#      perl splitFasta.pl [options]
#
#      Try 'perl splitFasta.pl -h' for more information.
#
#    Purpose: Splits a FASTA file into smaller chunks
#
#===============================================================================

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);

$| = 1; # Do not buffer output

my $VERSION     = '0.1';
my $DEBUG       = 0;

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print "splitFasta-".$VERSION."\n"; exit; },
            'i=s',
            's=f',
            'n=i'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

splitFasta

=head1 VERSION

splitFasta-0.1

=head1 SYNOPSIS

perl splitFasta.pl [-h] [-help] [-version] [-man] [-verbose] [options]

Examples:

perl splitFasta.pl -verbose -i file.fasta -s 2     #chunks of 2MB

perl splitFasta.pl -verbose -i file.fasta -n 10    #10 chunks

=head1 DESCRIPTION

Splits a FASTA file into smaller chunks.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<-i> <file>

Input file in FASTA format.

=item B<-s> <float>

Size of the FASTA file chunks in MB (MegaByte). Cannot be used in combination with -n.

=item B<-n> <integer>

Number of chunks the FASTA file should be separate into. Cannot be used in compination with -s and has to be greater than zero.

=back

=head1 AUTHOR

Robert SCHMIEDER, C<< <rschmieder_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <rschmieder_at_gmail_dot_com> >>.

=head1 COPYRIGHT

Copyright (C) 2013 Robert SCHMIEDER

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

#
################################################################################
## DATA AND PARAMETER CHECKING
################################################################################
#

my $file;

#Check if input file exists and check if file format is correct
if(exists $params{i}) {
    print "Checking input file\n" if(exists $params{verbose});
    if(-e $params{i}) {
        $file = $params{i};
        #check for file format
        my $format = &checkFileFormat($file);
        unless($format eq 'fasta') {
            &printError('input file for -i is in '.uc($format).' format not in FASTA format');
        }
    } else {
        &printError("could not find input file \"".$params{i}."\"");
    }
} else {
    &printError("you did not specify an input file containing the query sequences");
}

#check if split size and num are defined
if(exists $params{s} && exists $params{n}) {
    &printError("The options -s and -n cannot be used together");
}
if(!exists $params{s} && !exists $params{n}) {
    &printError("Either option -s or option -n has to be specified");
}

#Check for split size
if(exists $params{s}) {
    unless(looks_like_number($params{s}) && $params{s} > 0) {
        &printError("The split size has to be a float or integer number greater than zero");
    }
}

#Check for split num
if(exists $params{n}) {
    unless($params{n} =~ /^\d+$/ || $params{n} > 0) {
        &printError("The -n value has to be a number greater than zero");
    }
}


#
################################################################################
## DATA PROCESSING
################################################################################
#

my $seqnum = `grep -c '^>' $file`;

#status bar stuff
my ($numlines,$progress,$counter,$part);
if(exists $params{verbose}) {
    print STDERR "Estimate size of input data for status report (this might take a while for large files)\n";
    print STDERR "\tdone\n";
    $progress = 0;
    $counter = 1;
    $part = int($seqnum/100);
}
print STDERR "Split file into ".(exists $params{n} ? $params{n}.' ' : '')."chunks".(exists $params{s} ? " of max ".$params{s}." MB" : '')."\n" if(exists $params{verbose});
print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});

#read fasta file
open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < ".$file." |") or die "ERROR: Could not open file $file: $! \n";

#check if split by number or size
my $splitnum;
chomp($seqnum);
if(exists $params{n}) {
    $splitnum = int($seqnum/$params{n});
} else {
    my $filesize = -s "$file";
    $splitnum = int($seqnum/int(0.999999+($filesize/($params{s}*1000000))));
}
$splitnum += 1;
my $seqcounter = 0;
my $filecounter = 1;
my $tmpname = $params{i}.'_c'.$filecounter++.'.fasta';
open(OUT,">$tmpname") or die "ERROR: could not write to file $tmpname: $! \n";
open(IN,"<$file") or die "ERROR: could not open file $file: $! \n";
while(<IN>) {
    if(/^>/) {
        $seqcounter++;        
        #progress bar stuff
        $counter++;
        if($counter > $part) {
            $counter = 1;
            $progress++;
            $progress = 99 if($progress > 99);
            print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
        }
    }
    if($seqcounter > $splitnum) {
        close(OUT);
        $tmpname = $params{i}.'_c'.$filecounter++.'.fasta';
        open(OUT,">$tmpname") or die "ERROR: could not write to file $tmpname: $! \n";
        $seqcounter = 1;
    }
    print OUT $_;
    
}
close(IN);
close(OUT);

print STDERR "\r\tdone          \n" if(exists $params{verbose});


#
################################################################################
## MISC FUNCTIONS
################################################################################
#

sub printError {
    my $msg = shift;
    print STDERR "ERROR: ".$msg.".\n\nTry \'perl splitFasta.pl -h\' for more information.\nExit program.\n";
    exit(0);
}

sub checkFileFormat {
    my $file = shift;

    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    while (<FILE>) {
        chomp();
        next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && (/^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx\*-]+/ || /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/)) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && (/^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx\*-]+/ || /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/)) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        }
    }
    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}

