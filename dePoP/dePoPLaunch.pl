#!/usr/bin/perl
use 5.10.0;
use strict;
use warnings;

use File::Basename;
use File::Spec;
use File::Spec::Functions;
use File::Path qw(make_path);
use List::Util qw(max);

#---OPTIONS---
our $path_sep = File::Spec->catfile('', '');

our $help;
our $threads_number = 2;
our $config_file; # = "test.txt";                                               #debugging
our $gatk_location = catfile(catdir(dirname(__FILE__), "bin"), "GenomeAnalysisTK.jar"); 
our $bt2script_location = catfile(catdir(dirname(__FILE__), "bin"), "bt2alignment.sh");
our $depooler_location = catfile(catdir(dirname(__FILE__), "bin"), "s-DePooler-0.2.3.jar"); 
our ($output_folder, $prefix, $opt_file) = ('run');
our ($scheme, $fq_file, $dmp_folder, $algn_folder, $bam_file, $snp_folder, $vcf_file);
our ($demult_seq, @trim5_seq, @trim3_seq);
our ($bt2index, $reference, $regions);
our ($trimert, $trimqual, $trimlen) = (0.2, 25, 40);
our ($organism_ploidy) =(2);
our ($time_limit, $comb_limit, $tries, $success, $noIndels) = (2000);
our ($seqerr, $confidence, $mixing_precision);
our %pools_hash;
our $monitor;

our %OPT = (
    HELP1           => 'help', 
    HELP2           => 'h',
    THREADS         => 'n',
    CFG_FILE        => 'cfg',
    GATK            => 'gatk',
    SCHEME          => 'scheme',
    DEMULT_FOLDER   => 'dmpf',
    FQ_FILE         => 'fq',
    BAM_FILE        => 'bam',
    VCF_FILE        => 'vcf',
    PREFIX          => 'prefix',
    OUTPUT          => 'output',
    DEMULT          => 'demult',
    TR_SEQ5         => 'trseq5',
    TR_SEQ3         => 'trseq3',
    TR_ERROR        => 'trer',
    TR_QUAL         => 'trqual',
    READS_LEN       => 'readlen',
    BT2_INDEX       => 'bt2index',
    REFERENCE       => 'reference',
    REGIONS         => 'regions',
    ORGANISM_PLOIDY => 'org_ploidy',
    GRAFIC_MONITOR  => 'gm',
    TIME_LIMIT      => 'time_limit',
    COMBIN_LIMIT    => 'cmb_limit',
    TRIES           => 'tries',
    SUCCESS         => 'success',
    NO_INDELS       => 'noIndels',
    SEQERR          => 'seqerr',
    CONFIDENCE      => 'conf',
    MIX_PREC        => 'mp',
    );

#---LAUNCH---
&parse_options(@ARGV);
(&print_help() and die)                                         if ($help);
(&parse_options(&read_config_file()) and &parse_options(@ARGV)) if ($config_file);

die "No pooling scheme."    if (not defined($scheme) or not -e $scheme);
%pools_hash = &get_pools_sample_hash();
$prefix = $prefix ? "${prefix}\_" : '';

&create_output_folder();
&define_output_options_file_name();

if      ($vcf_file)     {&do_depooling();}
elsif   ($bam_file)     {&do_snpcalling();}
elsif   ($dmp_folder)   {&do_mapping();}
else                    {&do_demultiplexing();}

#---PROCEDURES---
sub do_demultiplexing {
    die "No fastq source file"              if not defined($fq_file)    or not -e $fq_file;
    die "No demultiplexing sequences file"  if not defined($demult_seq) or not -e $demult_seq;
    
    &create_demultiplexing_folder();
    &write_demultiplexing_options();
    say "***Demultiplexing***";
    
    my $exit_code = system "cutadapt -g file:$demult_seq -o ${dmp_folder}${path_sep}\{name}.fq $fq_file \\
    >> ${dmp_folder}${path_sep}dmp-report.txt";
    die if $exit_code
    
    &do_reads_trimming();
}

sub do_reads_trimming {
    say "***Reads Trimming***";

    foreach my $poolId (keys %pools_hash) {
        my $raw_fq = &define_raw_fastq_file_name($poolId);
        die "No fastq-file for the pool ${raw_fq}. Expected ${raw_fq}."    if not -e "${raw_fq}";
        my $tr_fq = &define_trimmed_fastq_file_name($poolId);
        
        my $command;
        foreach (@trim5_seq) {
            my $source = $command ? '-' : $raw_fq;
            my $part = fileparse($_, qr/\.[^.]*/);
            my $report_file = &define_trimming_report_file_name($poolId, $part);
            $command .= "cutadapt -g file:$_ -e ${trimert} --discard-untrimmed ${source} 2> ${report_file} | ";
        }
        foreach (@trim3_seq) {
            my $source = $command ? '-' : $raw_fq;
            my $part = fileparse($_, qr/\.[^.]*/);
            my $report_file = &define_trimming_report_file_name($poolId, $part);
            $command .= "cutadapt -a file:$_ -e ${trimert} --discard-untrimmed ${source} 2> ${report_file} | ";
        }
        {
            my $source = $command ? '-' : $raw_fq;
            my $report_file = &define_length_trimming_report_file_name($poolId);
            $command .= "cutadapt -q ${trimqual},${trimqual} -m ${trimlen} -o ${tr_fq} ${source} > ${report_file}";
        }
        
        say $poolId;
        my $exit_code = system $command;
        die if $exit_code;
    }

    &do_mapping();
}

sub do_mapping {
    die "No bowtie2 index"  if not defined($bt2index) or not -e "${bt2index}.1.bt2";
    
    &create_alignment_folder();
    &write_mapping_options();
    say "***Reads Mapping***";

    foreach my $poolId (keys %pools_hash) {
        say $poolId;
        my $tr_fq = &define_trimmed_fastq_file_name($poolId);
        die "No fastq-file for the pool ${poolId}. Expected ${tr_fq}." if not -e $tr_fq;
        my $pool_bam_file = &define_pool_bam_file_name($poolId);
        my $exit_code = system "${bt2script_location} ${pool_bam_file} ${bt2index} ${tr_fq} ${poolId} ${threads_number}";
        die if $exit_code;
    }

    &do_alignments_merging();
}

sub do_alignments_merging {
    say "***Alignments Merging***";
    
    $bam_file = &define_merged_bam_file_name();

    my $merging_command = "samtools merge -f ${bam_file}";
    foreach my $poolId (keys %pools_hash) {
        my $pool_bam_file = &define_pool_bam_file_name($poolId);
        die "No bam-file for the pool ${poolId}. Expected ${pool_bam_file}." if not -e $pool_bam_file;
        $merging_command .= " ${pool_bam_file}";
    }
    
    my $exit_code = system $merging_command;
    die if $exit_code;

    $exit_code = system "samtools index ${bam_file}";
    die if $exit_code;
    
    &do_snpcalling();
}

sub do_snpcalling {
    die "No reference"          if not defined($reference)  or not -e "${reference}";
    die "No regions"            if not defined($regions)    or not -e "${regions}";
    die "No bam-file"           if not defined($bam_file)   or not -e "${bam_file}";
    
    &create_snp_calling_folder();
    &write_snpcalling_options();
    say "***SNP-calling***";
    
    $vcf_file = &define_vcf_file_name();
    my $max_ploidy = (max values %pools_hash) * $organism_ploidy;
    my $snp_calling_report = &define_snp_calling_report_file_name();

    my $snp_calling_command = "java -jar ${gatk_location} -T HaplotypeCaller -gt_mode DISCOVERY "
        . "-R ${reference} -L ${regions} -I ${bam_file} --sample_ploidy ${max_ploidy} "
        . "-o ${vcf_file} -nct ${threads_number} 2> ${snp_calling_report}";

    my $exit_code = system $snp_calling_command;
    die if $exit_code;

    &do_depooling() if not $exit_code;
}

sub do_depooling {
    die "No vcf file"           if not defined($vcf_file) or not -e "${vcf_file}";

    &write_depooling_options();
    say "***De-pooling***";

    my $output_file = define_depooling_output_file_name();

    my $command = "java -jar ${depooler_location} -sch ${scheme} -vcf ${vcf_file} -n ${threads_number}";
    $command .= $monitor ? " -g" : (" > " . $output_file);
    $command .= " -timeLimit ${time_limit}"             if $time_limit;
    $command .= " -combLimit ${comb_limit}"             if $comb_limit;
    $command .= " -tries ${tries}"                      if $tries;
    $command .= " -success ${success}"                  if $success;
    $command .= " -noIndels"                            if $noIndels;
    $command .= " -seqerr ${seqerr}"                    if $seqerr;
    $command .= " -confidence ${confidence}"            if $confidence;
    $command .= " -mixingPrecision ${mixing_precision}" if $mixing_precision;

    my $exit_code = system $command;
    die if $exit_code;

    say "Finished successfully" ;
}

#---OPTIONS FUNCTIONS---
sub parse_options {
    use Getopt::Long qw(GetOptionsFromArray);
    
    my $opt_parsing = GetOptionsFromArray(\@_,
        "$OPT{HELP1}|$OPT{HELP2}"   => \$help,
        "$OPT{THREADS}=i"           => \$threads_number,
        "$OPT{CFG_FILE}=s"          => \$config_file,
        "$OPT{GATK}=s"              => \$gatk_location,
        "$OPT{SCHEME}=s"            => \$scheme,
        "$OPT{FQ_FILE}=s"           => \$fq_file,
        "$OPT{DEMULT_FOLDER}=s"     => \$dmp_folder,
        "$OPT{BAM_FILE}=s"          => \$algn_folder,
        "$OPT{VCF_FILE}=s"          => \$vcf_file,
        "$OPT{PREFIX}=s"            => \$prefix,
        "$OPT{OUTPUT}=s"            => \$output_folder,
        "$OPT{DEMULT}=s"            => \$demult_seq,
        "$OPT{TR_SEQ5}=s"           => \@trim5_seq,
        "$OPT{TR_SEQ3}=s"           => \@trim3_seq,
        "$OPT{TR_ERROR}=f"          => \$trimert,
        "$OPT{TR_QUAL}=i"           => \$trimqual, 
        "$OPT{READS_LEN}=i"         => \$trimlen,
        "$OPT{BT2_INDEX}=s"         => \$bt2index,
        "$OPT{REFERENCE}=s"         => \$reference, 
        "$OPT{REGIONS}=s"           => \$regions,
        "$OPT{ORGANISM_PLOIDY}=i"   => \$organism_ploidy,
        "$OPT{GRAFIC_MONITOR}"      => \$monitor,
        "$OPT{TIME_LIMIT}=i"        => \$time_limit,
        "$OPT{COMBIN_LIMIT}=i"      => \$comb_limit,
        "$OPT{TRIES}=i"             => \$tries,
        "$OPT{SUCCESS}=i"           => \$success,
        "$OPT{NO_INDELS}"           => \$noIndels,
        "$OPT{SEQERR}=i"            => \$seqerr,
        "$OPT{CONFIDENCE}=i"        => \$confidence,
        "$OPT{MIX_PREC}=i"          => \$mixing_precision,
    );
    die "Unknown option"  if (not $opt_parsing);
    @trim5_seq = split(/,/,join(',',@trim5_seq));
    @trim3_seq = split(/,/,join(',',@trim3_seq));
}

sub write_demultiplexing_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{FQ_FILE},          $fq_file                if $fq_file;
    printf OPTF "-%s\t%s\n",    $OPT{DEMULT},           $demult_seq;
    printf OPTF "-%s\t%s\n",    $OPT{DEMULT_FOLDER},    $dmp_folder;
    printf OPTF "-%s\t%s\n",    $OPT{TR_SEQ5},          join(',',@trim5_seq)    if @trim5_seq;
    printf OPTF "-%s\t%s\n",    $OPT{TR_SEQ3},          join(',',@trim3_seq)    if @trim3_seq;
    printf OPTF "-%s\t%s\n",    $OPT{TR_ERROR},         $trimert;
    printf OPTF "-%s\t%s\n",    $OPT{TR_QUAL},          $trimqual;
    printf OPTF "-%s\t%s\n",    $OPT{READS_LEN},        $trimlen;
    close OPTF;
}

sub write_mapping_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{BT2_INDEX},        $bt2index;
    close OPTF;
}

sub write_snpcalling_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{BAM_FILE},         $algn_folder;
    printf OPTF "-%s\t%s\n",    $OPT{REFERENCE},        $reference;
    printf OPTF "-%s\t%s\n",    $OPT{REGIONS},          $regions;
    printf OPTF "-%s\t%s\n",    $OPT{ORGANISM_PLOIDY},  $organism_ploidy;
    printf OPTF "-%s\t%s\n",    $OPT{GATK},             $gatk_location;
    close OPTF;
}

sub write_depooling_options {
    open(OPTF, ">>", $opt_file) or die "Cannot write options file .";
    printf OPTF "-%s\t%s\n",    $OPT{VCF_FILE},         $vcf_file;
    printf OPTF "-%s\t%s\n",    $OPT{SCHEME},           $scheme;
    printf OPTF "-%s\t%s\n",    $OPT{TIME_LIMIT},       $time_limit        if $time_limit;
    printf OPTF "-%s\t%s\n",    $OPT{COMBIN_LIMIT},     $comb_limit        if $comb_limit;
    printf OPTF "-%s\t%s\n",    $OPT{TRIES},            $tries             if $tries;
    printf OPTF "-%s\t%s\n",    $OPT{SUCCESS},          $success           if $success;
    printf OPTF "-%s\n",        $OPT{NO_INDELS},                           if $noIndels;
    printf OPTF "-%s\t%s\n",    $OPT{SEQERR},           $seqerr            if $seqerr;
    printf OPTF "-%s\t%s\n",    $OPT{CONFIDENCE},       $confidence        if $confidence;
    printf OPTF "-%s\t%s\n",    $OPT{MIX_PREC},         $mixing_precision  if $mixing_precision;
    printf OPTF "-%s\n",        $OPT{GRAFIC_MONITOR},   $monitor           if $monitor;

    printf OPTF "-%s\t%s\n",    $OPT{PREFIX},           $prefix             if $prefix;
    close OPTF;
}
#---FUNCTIONS---
sub print_help {
    say "i'll help you";
    die;
}

sub read_config_file {
    say "Configuration file is '$config_file'";
    my @array;
    if (open CONFIG, "<", $config_file) {
        while (<CONFIG>) {
            my @tokens = split( /\s+/, $_ );
            if (@tokens) {
                $tokens[0] =~ s/(^\w)/-$1/;
#                say $tokens[0];
                push (@array, @tokens);
            }
        }
        close CONFIG;
        @array;
    } else {
        die "Cannot read contig file .";
    };
}

sub get_pools_sample_hash {
    my %pools_hash = ();
    if (open FILE, "<", $scheme) {
        while (<FILE>) {
            if (! s/^\s+//) {
                my @ar = (split(/\s+/, $_));
                my $count = 0;
                for(my $i = 1; $i <= $#ar; $i++) {
                    ++$count if $ar[$i] != 0;
                }
                $pools_hash{$ar[0]} = $count;
            }
        }
        %pools_hash;
    } else {
        die "Cannot read scheme file.";
    }
}

#---NAMES---
sub create_output_folder {
    my $folder = canonpath($output_folder);
    my $index = 0;
    while(-d $folder) {
        $folder = canonpath($output_folder) . '_' . ++$index;
    }
    die "Cannot create output folder $folder"  if not make_path ($folder);
    $output_folder = $folder;
}
    
sub create_folder {
    my $name = $_[0];
    my $folder = catdir(canonpath($output_folder), $prefix . $name);
    if (not -d $folder) {
        my $s = make_path ($folder);
        die "Cannot create folder $folder"  if not $s;
    } else {
        die "Cannot write to $dmp_folder"   if ((stat($folder))[2] & 060) != 48; 
    }
    $folder;
}

sub define_output_options_file_name {
    $opt_file = catfile($output_folder, $prefix . "cfg.txt");
}
sub create_demultiplexing_folder {
    $dmp_folder = &create_folder("dmp-reads");
}
sub define_raw_fastq_file_name {
    my $poolId = $_[0];
    catfile($dmp_folder, $poolId . '.fq');
}
sub define_trimmed_fastq_file_name {
    my $poolId = $_[0];
    catfile($dmp_folder, $poolId . '_tr.fq');
}
sub define_trimming_report_file_name {
    my $poolId = $_[0];
    my $part = $_[1];
    catfile($dmp_folder, $poolId . '_' . $part .'_tr-rep.txt');
}
sub define_length_trimming_report_file_name {
    my $poolId = $_[0];
    catfile($dmp_folder, $poolId . '_filt-rep.txt');
}
sub create_alignment_folder {
    $algn_folder = &create_folder("alignments");
}
sub define_pool_bam_file_name {
    my $poolId = $_[0];
    catfile($algn_folder, ${prefix} . $poolId . '.bam');
}
sub define_merged_bam_file_name {
    catfile($algn_folder, ${prefix} . 'merged.bam');
}
sub create_snp_calling_folder {
    $snp_folder = &create_folder("snp-calling");
}
sub define_vcf_file_name {
    catfile($snp_folder, $prefix . 'snp.vcf');
}
sub define_snp_calling_report_file_name {
    catfile($snp_folder, $prefix . 'snp-calling_report.txt');
}
sub define_depooling_output_file_name {
    catfile($output_folder, $prefix . 'depooling.txt');
}

