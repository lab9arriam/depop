## dePoP - determination of polymorphism carriers from overlapping DNA pools

### Overview

dePoP is a pipeline designed for analysis of NGS-sequences of pooled samples. It identifies the carriers of rare Nucleotide Variants (NV) using sequence reads of overlapping pools, a process we called de-pooling. dePoP automatizes de-multiplexing, trimming, mapping, snp-calling of raw NGS-reads and performs de-pooling using s-dePooler - a novel java application.

dePoP is adapted for work in Linux OS.

### Installation

The pipeline is Perl-based (Perl v.5) it requires no installation, but the following tools must be present in the $PATH variable:

* Bowtie2 ver 2.2 or newer
* Cutadapt ver 1.8 or newer
* samtools and bcftools ver 1.3 or newer.

The following tools are included in repository and alternatively can be installed separately and included in pipeline by defined options:

* Genome Analysis Toolkit (GATK) .
* s-dePooler  

These two tools require Java8 to be installed in the system.
