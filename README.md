BAMClipper
==========
Remove gene-specific primer sequences from SAM/BAM alignments of PCR amplicons by soft-clipping

[Download latest version in a ZIP package](https://github.com/tommyau/bamclipper/zipball/master)

### Dependencies, as tested on 64-bit CentOS 5.5
* [SAMtools](http://www.htslib.org/download/) (at least version 1.3.1)
* [GNU Parallel](http://www.gnu.org/software/parallel/) (at least version 20130522)

### Usage
`bamclipper.sh` soft-clips gene-specific primers from BAM alignment file based on *genomic coordinates* of primer pairs in BEDPE format.
>./bamclipper.sh -b _BAM_ -p _BEDPE_ [-n _NTHREAD_] [-s _SAMTOOLS_] [-g _GNUPARALLEL_] [-u _UPSTREAM_] [-d _DOWNSTREAM_]

Given a BAM file called **_NAME_.bam**, a new BAM file (**_NAME_.primerclipped.bam**) and its associated index (**_NAME_.primerclipped.bam.bai**) will be generated in the current working directory.

_Notes_: For the sake of performance and simplicity, soft-clipping is performed solely based on genomic coordinates without involving the underlying sequence. Reference sequence names and coordinates of BAM and BEDPE are assumed to be derived from identical reference sequences (e.g. hg19).

*Required arguments*
- **-b** _FILE_: indexed BAM alignment file
- **-p** _FILE_: [BEDPE](http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) file of primer pair locations

*Options*
- **-n** _INT_: number of threads for clipprimer.pl (the workhorse Perl script of BAMClipper) and samtools sort [1]
- **-s** _FILE_: path to samtools executable [samtools]
- **-g** _FILE_: path to gnu parallel executable [parallel]
- **-u** _INT_: number of nucleotide upstream to 5' most nucleotide of primer (in addition to 5' most nucleotide of primer) for assigning alignments to primers based on the alignment starting position. [1]
- **-d** _INT_: number of nucleotide downstream to 5' most nucleotide of primer (in addition to 5' most nucleotide of primer) for assigning alignments to primers based on the alignment starting position. [5]

### Example using demo data
```bash
# Clip primers by BAMClipper
>./bamclipper.sh -b examples/SRR2075598.bam -p examples/trusight_myeloid.bedpe -n 4
# done!
# SRR2075598.primerclipped.bam and its index SRR2075598.primerclipped.bam.bai are generated.

# the new SRR2075598.primerclipped.bam should be identical to the provided example (compare checksum of alignments and ignore headers)
>samtools view SRR2075598.primerclipped.bam | md5sum
6a431457fd6e892646c17d1c3029c24e  -
>samtools view examples/SRR2075598.primerclipped.bam | md5sum
6a431457fd6e892646c17d1c3029c24e  -

# An example line of primer pair BEDPE file (an amplicon targeting ASXL1)
>grep 31022896 examples/trusight_myeloid.bedpe
chr20   31022896        31022921        chr20   31023096        31023123
```
Details of demo data:
- **examples/SRR2075598.bam**: original BWA-MEM v0.7.7 alignments of a TruSight Myeloid panel dataset (first 100 read pairs of [SRR2075598](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR2075598)) (reference: hg19)
- **examples/SRR2075598.primerclipped.bam**: primer-clipped alignments
- **examples/trusight_myeloid.bedpe**: primer locations derived from manufacturer's TruSight-Myeloid-Manifest.txt

Citation
--------
Au CH, Ho DN, Kwong A, Chan TL and Ma ESK, 2017. [BAMClipper: removing primers from alignments to minimize false-negative mutations in amplicon next-generation sequencing](http://www.nature.com/articles/s41598-017-01703-6). _Scientific Reports_ 7:1567  (doi:10.1038/s41598-017-01703-6)