#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;


#check if already created db
if notExists "${OUTDB}"; then
    if notExists "$1.dbtype"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "$@" "${OUTDB}" ${CREATEDB_PAR} \
            || fail "createdb failed"
    else
        cp "$1" "${OUTDB}"
        cp "$1.index" "${OUTDB}.index"
        cp "$1.lookup" "${OUTDB}.lookup"
        cp "$1.source" "${OUTDB}.source"
        cp "$1.dbtype" "${OUTDB}.dbtype"

        cp "$1_h" "${OUTDB}_h"
        cp "$1_h.index" "${OUTDB}_h.index"
        cp "$1_h.dbtype" "${OUTDB}_h.dbtype"
    fi

fi


# if notExists "${OUTDB}"; then
#     # shellcheck disable=SC2086
#     "${MMSEQS}" createdb "$@" "${OUTDB}" ${CREATEDB_PAR} \
#         || fail "createdb failed"
# fi


if [ "$("${MMSEQS}" dbtype "${OUTDB}")" = "Nucleotide" ]; then
    mv -f "${OUTDB}" "${OUTDB}_nucl"
    mv -f "${OUTDB}.index" "${OUTDB}_nucl.index"
    mv -f "${OUTDB}.lookup" "${OUTDB}_nucl.lookup"
    mv -f "${OUTDB}.source" "${OUTDB}_nucl.source"
    mv -f "${OUTDB}.dbtype" "${OUTDB}_nucl.dbtype"

    mv -f "${OUTDB}_h" "${OUTDB}_nucl_h"
    mv -f "${OUTDB}_h.index" "${OUTDB}_nucl_h.index"
    mv -f "${OUTDB}_h.dbtype" "${OUTDB}_nucl_h.dbtype"

    if notExists "${OUTDB}_nucl_contig_to_set.index"; then
        awk '{ print $1"\t"$3; }' "${OUTDB}_nucl.lookup" | sort -k1,1n -k2,2n > "${OUTDB}_nucl_contig_to_set.tsv"
        "${MMSEQS}" tsv2db "${OUTDB}_nucl_contig_to_set.tsv" "${OUTDB}_nucl_contig_to_set" "--output-dbtype" "5"\
            || fail "tsv2db failed"
    fi

    if notExists "${OUTDB}_nucl_set_to_contig.index"; then
        awk '{ print $3"\t"$1; }' "${OUTDB}_nucl.lookup" | sort -k1,1n -k2,2n > "${OUTDB}_nucl_set_to_contig.tsv"
        "${MMSEQS}" tsv2db "${OUTDB}_nucl_set_to_contig.tsv" "${OUTDB}_nucl_set_to_contig" "--output-dbtype" "5"\
            || fail "tsv2db failed"
    fi

    if notExists "${OUTDB}_nucl_orf.index"; then
        # shellcheck disable=SC2086
        if [ -n "${EXTRACTORFS_SPACER}" ]; then
            echo "Extract ORF using spacer parameters..."
            #"--reverse-frames" " "
            "${MMSEQS}" extractorfs "${OUTDB}_nucl" "${OUTDB}_nucl_orf" "--min-length" "9" "--orf-start-mode" "1" ${THREADS_PAR}\
                || fail "extractorfs failed"     
        else   
            "${MMSEQS}" extractorfs "${OUTDB}_nucl" "${OUTDB}_nucl_orf" "--orf-start-mode" "0" ${THREADS_PAR}\
                || fail "extractorfs failed"
        fi
    fi

    if notExists "${OUTDB}.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" translatenucs "${OUTDB}_nucl_orf" "${OUTDB}" ${TRANSLATENUCS_PAR} \
            || fail "translatenucs failed"
    fi

    # write extracted orfs locations on contig in alignment format
    if notExists "${OUTDB}_nucl_orf_aligned_to_contig.index"; then
        "${MMSEQS}" orftocontig "${OUTDB}_nucl" "${OUTDB}_nucl_orf" "${OUTDB}_nucl_orf_aligned_to_contig" \
            || fail "orftocontig step died"
    fi

    if notExists "${OUTDB}_nucl_orf_to_contig.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" filterdb "${OUTDB}_nucl_orf_aligned_to_contig" "${OUTDB}_nucl_orf_to_contig" --trim-to-one-column --filter-regex "^.*$" ${THREADS_PAR} \
            || fail "filterdb failed"
    fi

    if notExists "${OUTDB}_member_to_set.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" filterdb "${OUTDB}_nucl_orf_to_contig" "${OUTDB}_member_to_set" --mapping-file "${OUTDB}_nucl_contig_to_set.tsv" ${THREADS_PAR} \
            || fail "filterdb failed"
    fi

    if notExists "${OUTDB}_set_to_member.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" swapdb "${OUTDB}_member_to_set" "${OUTDB}_set_to_member" ${SWAPDB_PAR} \
            || fail "swapdb failed"
    fi

    if notExists "${OUTDB}_set_size.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" result2stats "${OUTDB}_nucl" "${OUTDB}_nucl" "${OUTDB}_set_to_member" "${OUTDB}_set_size" ${RESULT2STATS_PAR} \
            || fail "result2stats failed"
    fi

    if [ -n "${REVERSE_FRAGMENTS}" ]; then
        # shellcheck disable=SC2086
        echo  "Reversing database..."
        "${MMSEQS}" reverseseq "${OUTDB}" "${OUTDB}_reverse" "${THREADS_PAR}"\
            || fail "reverseseq died"
        mv -f "${OUTDB}_reverse" "${OUTDB}"
        mv -f "${OUTDB}_reverse.index" "${OUTDB}.index"
        mv -f "${OUTDB}_reverse.dbtype" "${OUTDB}.dbtype"
        #OUTDB="${OUTDB}_reverse"
        echo "Will base search on ${OUTDB}"
    fi
else
    fail "protein mode not implemented"
fi

if [ -n "${REMOVE_TMP}" ]; then
    echo "Remove temporary files"
    rmdir "${TMP_PATH}/search"
    if [ -n "${NUCL}" ]; then
        rm -f "${TMP_PATH}/nucl_set.tsv"
    fi
    rm -f "${TMP_PATH}/createsetdb.sh"
fi

