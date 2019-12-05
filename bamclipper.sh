#!/usr/bin/env bash
# bamclipper.sh
VERSION=1.1.2
NTHREAD=1
SAMTOOLS="samtools"
PARALLEL="parallel"
UPSTREAM=1
DOWNSTREAM=5
OUTPUT=$PWD
SAMTOOLS_VERSION_REQUIRED=1.3.1
PARALLEL_VERSION_REQUIRED=20130522


# show usage if no option is provided
if [[ "$#" -eq 0 ]]; then
    echo >&2 "Program: BAMClipper (Remove primer sequence from BAM alignments by soft-clipping)"
    echo >&2 "Version: $VERSION"
    echo >&2
    echo >&2 "File mode"
    echo >&2 "Usage: $0 -b BAM -p BEDPE [-n NTHREAD] [-s SAMTOOLS] [-g GNUPARALLEL] [-u UPSTREAM] [-d DOWNSTREAM]"
    echo >&2 "Required arguments:"
    echo >&2 "    -b FILE    indexed BAM alignment file"
    echo >&2 "    -p FILE    BEDPE file of primer pair locations"
    echo >&2
    echo >&2 "Pipe mode"
    echo >&2 "Usage: bwa mem ref.fasta r1.fastq r2.fastq | $0 -i -p BEDPE [-n NTHREAD] [-s SAMTOOLS] [-g GNUPARALLEL] [-u UPSTREAM] [-d DOWNSTREAM] | ..."
    echo >&2 "Required arguments:"
    echo >&2 "    -i         read SAM alignment from STDIN"
    echo >&2 "    -p FILE    BEDPE file of primer pair locations"
    echo >&2
    echo >&2 "Options for either mode:"
    echo >&2 "    -n INT     number of threads for clipprimer.pl and samtools sort [$NTHREAD]"
    echo >&2 "    -s FILE    path to samtools executable [$SAMTOOLS]"
    echo >&2 "    -g FILE    path to gnu parallel executable [$PARALLEL]"
    echo >&2 "    -u INT     number of nucleotide upstream to 5' most nucleotide of primer [$UPSTREAM]"
    echo >&2 "    -d INT     number of nucleotide downstream to 5' most nucleotide of primer [$DOWNSTREAM]"
    echo >&2 "    -o DIR     path to write output [$OUTPUT]"
    exit 1
fi

function error {
    echo >&2 "ERROR: $1"
    exit 1
}

while getopts ":ib:p:n:s:g:u:d:o:" o; do
    case "${o}" in
	i)
	    PIPE=1
	    ;;
	b)
	    BAM=${OPTARG}
	    [[ ! -f "$BAM" ]] && error "BAM file is not found ($BAM)"
	    [[ ! -f "$BAM.bai" ]] && error "BAM file index $BAM.bai is not found. (Tips: $SAMTOOLS index $BAM)"
	    BAMbn=$(basename "$BAM")
	    ;;
	p)
	    BEDPE=${OPTARG}
	    [[ -f "$BEDPE" ]] || error "BEDPE file is not found ($BEDPE)"
	    ;;
	n)
	    NTHREAD=${OPTARG}
	    [[ "$NTHREAD" -ge 1 ]] || error "NTHREAD requires non-zero integer"
	    ;;
	s)
	    SAMTOOLS=${OPTARG}
	    ;;
	g)
	    PARALLEL=${OPTARG}
	    ;;
	u)
	    UPSTREAM=${OPTARG}
	    [[ "$UPSTREAM" =~ ^[0-9]+$ ]] || error "UPSTREAM requires non-negative integer"
	    [[ "$UPSTREAM" -ge 0 ]] || error "UPSTREAM requires non-negative integer"
	    ;;
	d)
	    DOWNSTREAM=${OPTARG}
	    [[ "$DOWNSTREAM" =~ ^[0-9]+$ ]] || error "DOWNSTREAM requires non-negative integer"
	    [[ "$DOWNSTREAM" -ge 0 ]] || error "DOWNSTREAM requires non-negative integer"
	    ;;
	o)
	    OUTPUT=${OPTARG}
	    [[ ! -d "$OUTPUT" ]] && error "OUTPUT is not a valid directory ($OUTPUT)."
	    ;;
	*)
	    error "Invalid option: -$OPTARG"
	    ;;
    esac
done
shift $((OPTIND-1))

# assert (either BAM file or PIPE mode, but not both) and BEDPE are defined
([[ -z "$BAM" ]] && [[ -z "$PIPE" ]]) && error "BAM file (-b) or Pipe mode (-i) is not defined"
([[ ! -z "$BAM" ]] && [[ "$PIPE" == 1 ]]) && error "File mode (-b) and Pipe mode (-i) cannot be defined simultaneously"
[[ -z "$BEDPE" ]] && error "BEDPE file (-p) is not defined"

# Absolute path for output
OUTPUT=$(readlink -e $OUTPUT)

# check parallel & version
"$PARALLEL" --minversion $PARALLEL_VERSION_REQUIRED >/dev/null 2>&1 || error "GNU Parallel (provided path: $PARALLEL) is not running properly. Please check the path and/or version (at least $PARALLEL_VERSION_REQUIRED)."

# run bamclipper
SCRIPT_PATH="$(readlink -f $0)"
SCRIPT_DIR="$(dirname $SCRIPT_PATH)"
if [ "$PIPE" == 1 ]; then
    # pipe mode
    cat /dev/stdin | "$SCRIPT_DIR"/injectseparator.pl | "$PARALLEL" -j "$NTHREAD" --keep-order --remove-rec-sep --pipe --remove-rec-sep --recend '__\n' --block 10M "$SCRIPT_DIR/clipprimer.pl --in $BEDPE --upstream $UPSTREAM --downstream $DOWNSTREAM"
else
    # file mode
    # check samtools & version
    function version { echo "$@" | cut -f1 -d"+" | awk -F. '{ printf("%03d%03d%03d\n", $1,$2,$3); }'; }
    "$SAMTOOLS" --version-only >/dev/null 2>&1 || error "SAMtools (provided path: $SAMTOOLS) is not running properly. Please check the path and/or version (at least $SAMTOOLS_VERSION_REQUIRED)"
    SAMTOOLS_VERSION=`"$SAMTOOLS" --version-only`
    [ "$(version "$SAMTOOLS_VERSION")" -lt "$(version "$SAMTOOLS_VERSION_REQUIRED")" ] && error "SAMtools version ($SAMTOOLS_VERSION) is not supported (supported version: at least $SAMTOOLS_VERSION_REQUIRED)."

    "$SAMTOOLS" collate -O --output-fmt SAM "$BAM" "${BAMbn}.sort1" | "$SCRIPT_DIR"/injectseparator.pl | "$PARALLEL" -j "$NTHREAD" --keep-order --remove-rec-sep --pipe --remove-rec-sep --recend '__\n' --block 10M "$SCRIPT_DIR/clipprimer.pl --in $BEDPE --upstream $UPSTREAM --downstream $DOWNSTREAM" | "$SAMTOOLS" sort -T "${BAMbn}.sort2" -@ "$NTHREAD" > "${OUTPUT}/${BAMbn%.bam}.primerclipped.bam" && "$SAMTOOLS" index "${OUTPUT}/${BAMbn%.bam}.primerclipped.bam"
fi
