#!/usr/bin/env bash
# bamclipper.sh
VERSION=1.0.1
NTHREAD=1
SAMTOOLS="samtools"
PARALLEL="parallel"
UPSTREAM=1
DOWNSTREAM=5
SAMTOOLS_VERSION_REQUIRED=1.3.1
PARALLEL_VERSION_REQUIRED=20130522

function usage {
    echo >&2 "Program: BAMClipper (Remove primer sequence from BAM alignments by soft-clipping)"
    echo >&2 "Version: $VERSION"
    echo >&2
    echo >&2 "Usage: $0 -b BAM -p BEDPE [-n NTHREAD] [-s SAMTOOLS] [-g GNUPARALLEL] [-u UPSTREAM] [-d DOWNSTREAM]"
    echo >&2
    echo >&2 "Required arguments:"
    echo >&2 "    -b FILE    indexed BAM alignment file"
    echo >&2 "    -p FILE    BEDPE file of primer pair locations"
    echo >&2
    echo >&2 "Options:"
    echo >&2 "    -n INT     number of threads for clipprimer.pl and samtools sort [$NTHREAD]"
    echo >&2 "    -s FILE    path to samtools executable [$SAMTOOLS]"
    echo >&2 "    -g FILE    path to gnu parallel executable [$PARALLEL]"
    echo >&2 "    -u INT     number of nucleotide upstream to 5' most nucleotide of primer [$UPSTREAM]"
    echo >&2 "    -d INT     number of nucleotide downstream to 5' most nucleotide of primer [$DOWNSTREAM]"
    exit 1
}


while getopts ":b:p:n::s::g::u::d::" o; do
    case "${o}" in
	b)
	    BAM=${OPTARG}
	    BAMbn=$(basename "$BAM")
	    if [[ ! -f "$BAM" ]]; then
		echo >&2 "ERROR: BAM file $BAM cannot be found."
		exit 1
	    fi 
	    if [[ ! -f "$BAM.bai" ]]; then
		echo >&2 "ERROR: BAM file index $BAM.bai cannot be found."
		echo >&2 "Tips: $SAMTOOLS index $BAM"
		exit 1
	    fi 
	    ;;
	p)
	    BEDPE=${OPTARG}
	    [[ -f "$BEDPE" ]] || usage
	    ;;
	n)
	    NTHREAD=${OPTARG}
	    [[ "$NTHREAD" -ge 1 ]] || usage
	    ;;
	s)
	    SAMTOOLS=${OPTARG}
	    ;;
	g)
	    PARALLEL=${OPTARG}
	    ;;
	u)
	    UPSTREAM=${OPTARG}
	    [[ "$UPSTREAM" -ge 0 ]] || usage
	    ;;
	d)
	    DOWNSTREAM=${OPTARG}
	    [[ "$DOWNSTREAM" -ge 0 ]] || usage
	    ;;
	*)
	    usage
	    ;;
    esac
done
shift $((OPTIND-1))

# assert: BAM and BEDPE are defined
if [ -z "$BAM" ] || [ -z "$BEDPE" ]; then
    usage
fi

# check parallel & version
"$PARALLEL" --minversion $PARALLEL_VERSION_REQUIRED >/dev/null 2>&1 || { echo >&2 "ERROR: GNU Parallel (provided path: $PARALLEL) is not running properly. Please check the path and/or version (at least $PARALLEL_VERSION_REQUIRED)."; exit 1; }

# check samtools & version
function version { echo "$@" | cut -f1 -d"+" | awk -F. '{ printf("%03d%03d%03d\n", $1,$2,$3); }'; }
"$SAMTOOLS" --version-only >/dev/null 2>&1 || { echo >&2 "ERROR: SAMtools (provided path: $SAMTOOLS) is not running properly. Please check the path and/or version (at least $SAMTOOLS_VERSION_REQUIRED)"; exit 1; }
SAMTOOLS_VERSION=`"$SAMTOOLS" --version-only`
if [ "$(version "$SAMTOOLS_VERSION")" -lt "$(version "$SAMTOOLS_VERSION_REQUIRED")" ]; then
     echo >&2 "ERROR: SAMtools version ($SAMTOOLS_VERSION) is not supported (supported version: at least $SAMTOOLS_VERSION_REQUIRED)."
     exit 1
fi

# run bamclipper
SCRIPT_PATH="$(readlink -f $0)"
SCRIPT_DIR="$(dirname $SCRIPT_PATH)"
"$SAMTOOLS" collate -O --output-fmt SAM "$BAM" "${BAMbn}.sort1" | "$SCRIPT_DIR"/injectseparator.pl | "$PARALLEL" -j "$NTHREAD" --keep-order --remove-rec-sep --pipe --remove-rec-sep --recend '__\n' --block 1m "$SCRIPT_DIR/clipprimer.pl --in $BEDPE --upstream $UPSTREAM --downstream $DOWNSTREAM" | "$SAMTOOLS" sort -T "${BAMbn}.sort2" -@ "$NTHREAD" > "${BAMbn%.bam}.primerclipped.bam" && "$SAMTOOLS" index "${BAMbn%.bam}.primerclipped.bam"
