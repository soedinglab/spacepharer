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
[ "$#" -ne 5 ] && echo "Please provide <queryDB> <targetDB> <controltargetDB> <outputDB> <tmpDir>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[ ! -f "$3.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
# TO DO??? add check if $3.dbtype already exists before entire workfolw ???

QUERY="$1"
TARGET="$2"
CONTROLTARGET="$3"
OUTPUT="$4"
TMP_PATH="$5"

if notExists "${TMP_PATH}/result.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" search "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/search" ${SEARCH_PAR} \
        || fail "search failed"
fi

if notExists "${TMP_PATH}/aggregate.index"; then
    # aggregation: take for each target set the best hit
    # shellcheck disable=SC2086
    "${MMSEQS}" besthitperset "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/aggregate" ${BESTHITBYSET_PAR} \
        || fail "aggregate best hit failed"
fi

if notExists "${TMP_PATH}/aggregate_merged.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_merged" ${THREADS_PAR} \
        || fail "mergesetresults failed"
fi

if notExists "${OUTPUT}.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" combinepvalperset "${QUERY}" "${TARGET}" "${TMP_PATH}/aggregate_merged" "${OUTPUT}" ${COMBINEPVAL_PAR} \
        || fail "combinepvalperset failed"
fi

# if notExists "${OUTPUT}.index"; then
#     # shellcheck disable=SC2086
#     "${MMSEQS}" mergeresultsbyset "${QUERY}_set_to_member" "${TMP_PATH}/aggregate" "${OUTPUT}" ${THREADS_PAR} \
#         || fail "mergesetresults failed"
# fi
#if flag --output-alignment
# if notExists "${TMP_PATH}/aggregate_offset.index"; then
#     # shellcheck disable=SC2086
#     "${MMSEQS}" offsetalignment "${QUERY}_nucl" "${QUERY}" "${TARGET}_nucl" "${TARGET}" "${TMP_PATH}/aggregate" "${TMP_PATH}/aggregate_offset" ${OFFSETALIGNMENT_PAR} \
#         || fail "offsetalignment failed"
# fi

# if notExists "${TMP_PATH}/aln.index"; then
#     # shellcheck disable=SC2086
#     "${MMSEQS}" convertalis "${QUERY}_nucl" "${TARGET}_nucl" "${TMP_PATH}/aggregate_offset" "${TMP_PATH}/aln" "--db-output" "--search-type" "3" ${THREADS_PAR} \
#         || fail "convertalis failed"
# fi

# if notExists "${TMP_PATH}/aln_merge.index"; then
#     # shellcheck disable=SC2086
#     "${MMSEQS}" mergeresultsbyset "${QUERY}_nucl_set_to_contig" "${TMP_PATH}/aln" "${TMP_PATH}/aln_merge" ${THREADS_PAR} \
#         || fail "mergeresultsbyset failed"
# fi



if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir "${TMP_PATH}/search"
    "$MMSEQS" rmdb "${TMP_PATH}/result"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate"
    "$MMSEQS" rmdb "${TMP_PATH}/aggregate_merged"
    rm -f "${TMP_PATH}/multihitsearch.sh"
fi



# #!/bin/bash -e

# # predict exons workflow script
# fail() {
#     echo "Error: $1"
#     exit 1
# }

# notExists() {
#     [ ! -f "$1" ]
# }

# abspath() {
#     if [ -d "$1" ]; then
#         (cd "$1"; pwd)
#     elif [ -f "$1" ]; then
#         if [[ "$1" == */* ]]; then
#             echo "$(cd "${1%/*}"; pwd)/${1##*/}"
#         else
#             echo "$(pwd)/$1"
#         fi
#     elif [ -d "$(dirname "$1")" ]; then
#             echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
#     fi
# }

# # check number of input variables
# [ "$#" -ne 4 ] && echo "Please provide <contigsDB> <proteinsDB> <predictmatchBaseName> <tmpDir>" && exit 1;
# # check if files exist
# [ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
# [ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
# [ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

# INPUT_CONTIGS="$(abspath "$1")"
# INPUT_TARGET_PROTEINS="$(abspath "$2")"
# TMP_PATH="$(abspath "$4")"

# # extract coding fragments from input contigs (result in DNA)
# if notExists "${TMP_PATH}/nucl_6f.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" extractorfs "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" ${EXTRACTORFS_PAR} \
#         || fail "extractorfs step died"
# fi

# # write extracted orfs locations on contig in alignment format
# if notExists "${TMP_PATH}/nucl_6f_orf_aligned_to_contig.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" orftocontig "${INPUT_CONTIGS}" "${TMP_PATH}/nucl_6f" "${TMP_PATH}/nucl_6f_orf_aligned_to_contig" ${THREADS_PAR} \
#         || fail "orftocontig step died"
# fi

# # translate each coding fragment (result in AA)
# if notExists "${TMP_PATH}/aa_6f.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" translatenucs "${TMP_PATH}/nucl_6f" "${TMP_PATH}/aa_6f" ${TRANSLATENUCS_PAR} \
#         || fail "translatenucs step died"
# fi

# # when running in null mode (to assess evalues), we reverse the AA fragments:
# AA_FRAGS="${TMP_PATH}/aa_6f"
# if [ -n "$REVERSE_FRAGMENTS" ]; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" reverseseq "${AA_FRAGS}" "${AA_FRAGS}_reverse" ${THREADS_PAR} \
#         || fail "reverseseq step died"
#     AA_FRAGS="${AA_FRAGS}_reverse"
#     echo "Will base search on ${AA_FRAGS}"
# fi

# # search with each aa fragment against a target DB (result has queries as implicit keys)
# if notExists "${TMP_PATH}/search_res.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" search "${AA_FRAGS}" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/tmp_search" ${SEARCH_PAR} \
#         || fail "search step died"
# fi

# # swap results (result has targets as implicit keys)
# if notExists "${TMP_PATH}/search_res_swap.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" swapresults "${AA_FRAGS}" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/search_res" "${TMP_PATH}/search_res_swap" ${SWAPRESULT_PAR} \
#         || fail "swap step died"
# fi

# # join contig information to swapped results (result has additional info about the origin of the AA fragments)
# if notExists "${TMP_PATH}/search_res_swap_w_contig_info.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" filterdb "${TMP_PATH}/search_res_swap" "${TMP_PATH}/search_res_swap_w_contig_info" --join-db "${TMP_PATH}/nucl_6f_orf_aligned_to_contig" --filter-column 1 ${THREADS_PAR} \
#         || fail "filterdb (to join contig info) step died"
# fi

# # sort joined swapped results by contig id
# if notExists "${TMP_PATH}/search_res_swap_w_contig_info_sorted.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" filterdb "${TMP_PATH}/search_res_swap_w_contig_info" "${TMP_PATH}/search_res_swap_w_contig_info_sorted" --sort-entries 1 --filter-column 11 ${THREADS_PAR} \
#         || fail "filterdb (to sort by contig) step died"
# fi

# # for each target, with respect to each contig and each strand, find the optimal set of exons
# if notExists "${TMP_PATH}/dp_protein_contig_strand_map.dbtype"; then
#     # shellcheck disable=SC2086
#     "$MMSEQS" collectoptimalset "${TMP_PATH}/search_res_swap_w_contig_info_sorted" "${INPUT_TARGET_PROTEINS}" "${TMP_PATH}/" ${COLLECTOPTIMALSET_PAR} \
#         || fail "collectoptimalset step died"
# fi

# # post processing
# mv -f "${TMP_PATH}/dp_protein_contig_strand_map" "$3_dp_protein_contig_strand_map" || fail "Could not move result to $3_dp_protein_contig_strand_map"
# mv -f "${TMP_PATH}/dp_protein_contig_strand_map.index" "$3_dp_protein_contig_strand_map.index" || fail "Could not move result to $3_dp_protein_contig_strand_map.index"
# mv -f "${TMP_PATH}/dp_protein_contig_strand_map.dbtype" "$3_dp_protein_contig_strand_map.dbtype" || fail "Could not move result to $3_dp_protein_contig_strand_map.dbtype"

# mv -f "${TMP_PATH}/dp_optimal_exon_sets" "$3_dp_optimal_exon_sets" || fail "Could not move result to $3_dp_optimal_exon_sets"
# mv -f "${TMP_PATH}/dp_optimal_exon_sets.index" "$3_dp_optimal_exon_sets.index" || fail "Could not move result to $3_dp_optimal_exon_sets.index"
# mv -f "${TMP_PATH}/dp_optimal_exon_sets.dbtype" "$3_dp_optimal_exon_sets.dbtype" || fail "Could not move result to $3_dp_optimal_exon_sets.dbtype"


# if [ -n "$REMOVE_TMP" ]; then
#     echo "Removing temporary files from ${TMP_PATH}"
#     rm -f "${TMP_PATH}"/nucl_6f*
#     rm -f "${TMP_PATH}"/aa_6f*
#     rm -f "${TMP_PATH}"/search_res*
#     rm -r "${TMP_PATH}"/tmp_search/
# fi

