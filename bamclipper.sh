#!/usr/bin/env bash
# bamclipper.sh

usage() { echo "Usage: $0 -b BAM -p BEDPE [-n NTHREAD] [-s SAMTOOLS] [-g GNUPARALLEL] [-u UPSTREAM] [-d DOWNSTREAM]" 1>&2; exit 1; }

NTHREAD=1
SAMTOOLS="samtools"
PARALLEL="parallel"
UPSTREAM=1
DOWNSTREAM=5

while getopts ":b:p:n::s::g::u::d::" o; do
    case "${o}" in
        b)
	    BAM=${OPTARG}
	    BAMbn="$(basename $BAM)"
            [[ -f "$BAM" && -f "$BAM.bai" ]] || usage
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

if [ -z "$BAM" ] || [ -z "$BEDPE" ]; then
    usage
fi

SCRIPT_PATH="$(readlink -f $0)"
SCRIPT_DIR="$(dirname $SCRIPT_PATH)"

"$SAMTOOLS" collate -O --output-fmt SAM $BAM ${BAMbn}.sort1 | "$SCRIPT_DIR"/injectseparator.pl | "$PARALLEL" -j "$NTHREAD" --keep-order --remove-rec-sep --pipe --remove-rec-sep --recend '__\n' --block 1m "$SCRIPT_DIR/clipprimer.pl --in $BEDPE --upstream $UPSTREAM --downstream $DOWNSTREAM" | "$SAMTOOLS" sort -T ${BAMbn}.sort2 -@ "$NTHREAD" > ${BAMbn%.bam}.primerclipped.bam && "$SAMTOOLS" index ${BAMbn%.bam}.primerclipped.bam
