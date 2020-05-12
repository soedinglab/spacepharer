#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check number of input variables
#[ "$#" -ne 5 ] && echo "Please provide <queryDB> <targetDB> <controltargetDB> <outputDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -f "$3.dbtype" ] && echo "$3.dbtype not found!" && exit 1;

QUERY="$1"
TARGET="$2"
CONTROLTARGET="$3"
OUTPUT="$4"
TMP_PATH="$5"

#search in target db
if notExists "${TMP_PATH}/result.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
fi

if notExists "${TMP_PATH}/aggregate.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitperset "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${TMP_PATH}/aggregate_merged.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_merged" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

if notExists "${TMP_PATH}/cScore.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" combinepvalperset "${QUERY}" "${TARGET}" "${TMP_PATH}/aggregate_merged" "${TMP_PATH}/cScore" "${TMP_PATH}" ${COMBINEPVALPERSET_PAR} \
        || fail "combinepvalperset failed"
fi

#search in control target db
if notExists "${TMP_PATH}/result_rev.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERY}" "${CONTROLTARGET}" "${TMP_PATH}/result_rev" "${TMP_PATH}/search_rev" ${SEARCH_PAR} \
        || fail "search failed"
fi

if notExists "${TMP_PATH}/aggregate_rev.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitperset "${QUERY}" "${CONTROLTARGET}" "${TMP_PATH}/result_rev" "${TMP_PATH}/aggregate_rev" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${TMP_PATH}/aggregate_merged_rev.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate_rev" "${TMP_PATH}/aggregate_merged_rev" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

if notExists "${TMP_PATH}/cScore_rev.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" combinepvalperset "${QUERY}" "${CONTROLTARGET}" "${TMP_PATH}/aggregate_merged_rev" "${TMP_PATH}/cScore_rev" "${TMP_PATH}" ${COMBINEPVALPERSET_PAR} \
        || fail "combinepvalperset failed"
fi

#fdr analysis
if notExists "${TMP_PATH}/match.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" filtermatchbyfdr "${TMP_PATH}/cScore" "${TMP_PATH}/cScore_rev" "${TMP_PATH}/match" ${FILTERMATCHBYFDR_PAR} \
        || fail "filtermatchbyfdr failed"
fi

#nucleotide alignment
if notExists "${TMP_PATH}/aggregate_truncated.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" truncatebesthits "${QUERY}" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_truncated" ${THREADS_PAR} \
        || fail "truncatebesthits failed"
fi

if notExists "${TMP_PATH}/aggregate_offset.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" offsetalignment "${QUERY}_nucl" "${QUERY}" "${TARGET}_nucl" "${TARGET}" "${TMP_PATH}/aggregate_truncated" "${TMP_PATH}/aggregate_offset" --search-type 4 ${THREADS_PAR} \
        || fail "offsetalignment failed"
fi

if notExists "${TMP_PATH}/aln.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" convertalis "${QUERY}_nucl" "${TARGET}_nucl" "${TMP_PATH}/aggregate_offset" "${TMP_PATH}/aln" \
        --db-output --search-type 3 --format-output "tsetid,query,qset,target,evalue,qstart,qend,qlen,tstart,tend,qaln,taln" ${THREADS_PAR} \
        || fail "convertalis failed"
fi

if notExists "${TMP_PATH}/aln_merge.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_nucl_set_to_contig" "${TMP_PATH}/aln" "${TMP_PATH}/aln_merge" ${THREADS_PAR} \
        || fail "mergeresultsbyset failed"
fi

if [ -n "$REPORT_PAM" ]; then
    if notExists "${TMP_PATH}/aln_merge_pam.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" findpam "${TARGET}_nucl" "${TMP_PATH}/aln_merge" "${TMP_PATH}/aln_merge_pam" ${THREADS_PAR} \
            || fail "findpam failed"
    fi

    if notExists "${OUTPUT}.tsv"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" summarizeresults "${TMP_PATH}/match" "${TMP_PATH}/aln_merge_pam" "${OUTPUT}.tsv" ${SUMMARIZERESULTS_PAR} \
            || fail "summarizeresults failed"
    fi

else
    if notExists "${OUTPUT}.tsv"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" summarizeresults "${TMP_PATH}/match" "${TMP_PATH}/aln_merge" "${OUTPUT}.tsv" ${SUMMARIZERESULTS_PAR} \
            || fail "summarizeresults failed"
    fi
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir "${TMP_PATH}/search"
    "$MMSEQS" rmdb "${TMP_PATH}/result"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_merged"
    "$MMSEQS" rmdb "${TMP_PATH}/cScore"
    "$MMSEQS" rmdb "${TMP_PATH}/result_rev"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_rev"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_merged_rev"
    "$MMSEQS" rmdb "${TMP_PATH}/cScore_rev"
    "$MMSEQS" rmdb "${TMP_PATH}/match"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_truncated"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_offset"
    "$MMSEQS" rmdb "${TMP_PATH}/aln"
    "$MMSEQS" rmdb "${TMP_PATH}/aln_merge"
    if [ -n "$REPORT_PAM" ]; then 
        "$MMSEQS" rmdb "${TMP_PATH}/aln_merge_pam"
    fi
    "$MMSEQS" rmdb "${TMP_PATH}/output"
    rm -f "${TMP_PATH}/output.tsv"
    rm -f "${TMP_PATH}/predictmatch.sh"
fi
