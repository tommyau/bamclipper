BAMClipper
==========
Remove gene-specific primer sequences from SAM/BAM alignments of PCR amplicons by soft-clipping

[Download latest version in a ZIP package](https://github.com/tommyau/bamclipper/zipball/master)

### Dependencies, as tested on 64-bit CentOS 5.5
* [SAMtools](http://www.htslib.org/download/) (version 1.3 tested)
* [GNU Parallel](http://www.gnu.org/software/parallel/) (version 20130522 tested)

### Usage
`bamclipper.sh` soft-clips gene-specific primers from BAM alignment file based on *genomic coordinates* of primer pairs in BEDPE format.
>./bamclipper.sh -b _BAM_ -p _BEDPE_ [-n _NTHREAD_] [-s _SAMTOOLS_] [-g _GNUPARALLEL_] [-u _UPSTREAM_] [-d _DOWNSTREAM_]

Given a BAM file called **_NAME_.bam**, a new BAM file (**_NAME_.primerclipped.bam**) and its associated index (**_NAME_.primerclipped.bam.bai**) will be generated in the current working directory.

_Notes_: For the sake of performance and simplicity, soft-clipping is performed solely based on genomic coordinates without involving the underlying sequence. Reference sequence names and coordinates of BAM and BEDPE are assumed to be derived from identical reference sequences (e.g. hg19).

*Required arguments*
- **-b** _FILE_: indexed BAM alignment file
- **-p** _FILE_: [BEDPE](http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) file of primer pair locations

*Options*
- **-n** _INT_: number of threads for clipprimer.pl (the workhorse Perl script of BAMClipper) and two samtools sort instances [1]
- **-s** _FILE_: path to samtools executable [samtools]
- **-g** _FILE_: path to gnu parallel executable [parallel]
- **-u** _INT_: number of nucleotide upstream to 5' most nucleotide of primer (in addition to 5' most nucleotide of primer) for assigning alignments to primers based on the alignment starting position. [1]
- **-d** _INT_: number of nucleotide downstream to 5' most nucleotide of primer (in addition to 5' most nucleotide of primer) for assigning alignments to primers based on the alignment starting position. [5]

### Example
```bash
# Clip primers by BAMClipper
>./bamclipper.sh -b sample1.bam -p trusight_myeloid.bedpe -n 12
# done!
# sample1.primerclipped.bam and its index sample1.primerclipped.bam.bai will be generated.

# Show an example line of primer pair BEDPE file
>head -1 trusight_myeloid.bedpe
chr1	36931667	36931695	chr1	36931911	36931937
```

Citation
--------
Au CH, Ho DN, Kwong A, Chan TL and Ma ESK, 2016. _(submitted)_
