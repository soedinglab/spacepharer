#!/bin/sh -e

fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;

hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand wget
#hasCommand awk
#hasCommand zcat
hasCommand touch
#hasCommand tar

OUT_PATH="$1"

INPUTS=""
if notExists "${OUT_PATH}/phage_download.complete"; then
    echo "Download phage genomes"
    while read NAME URL; 
    do
        wget -nv -O "${OUT_PATH}/${NAME}" "${URL}" ;        
        INPUTS="${OUT_PATH}/${NAME} ${INPUTS}";
    done < /home/rosy/metaeuk/data/phagegenome.txt
    touch "${OUT_PATH}/phage_download.complete" \
        || fail "download failed"
    
else
    while read NAME URL;
    do
        INPUTS="${OUT_PATH}/${NAME} ${INPUTS}";
    done < /home/rosy/metaeuk/data/phagegenome.txt
fi

if notExists "${OUTDB}.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsetdb ${INPUTS} "${OUTDB}" "${TMP_PATH}" \
        || fail "createsetdb failed"
fi

if [ -n "$CREATE_REVERSE_SETDB" ] && notExists "${OUTDB}_reverse.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsetdb ${INPUTS} "${OUTDB}_reverse" "${TMP_PATH}" "--reverse-fragments" "1" \
        || fail "create reverse setdb failed"
fi

# if [ "$DOWNLOAD_MAPPING" -eq "1" ]; then
#     # Download the latest UniProt ID mapping to extract taxon identifiers
#     if notExists "${TMP_PATH}/mapping_download.complete"; then
#         echo "Download idmapping.dat.gz"
#         URL="ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
#         wget -nv -O - "$URL" | zcat | awk '$2 == "NCBI_TaxID" {print $1"\t"$3 }' > "${TMP_PATH}/taxidmapping"
#         touch "${TMP_PATH}/mapping_download.complete"
#     fi
#     MAPPINGFILE="${TMP_PATH}/taxidmapping"
# fi
# # create mapping
# if notExists "${TMP_PATH}/targetDB_mapping.complete"; then
#     awk 'NR == FNR { f[$1] = $2; next } $2 in f { print $1"\t"f[$2] }' \
#         "$MAPPINGFILE" "${TAXDBNAME}.lookup" > "${TMP_PATH}/targetDB_mapping"
#     touch "${TMP_PATH}/targetDB_mapping.complete"
# fi

# # finalize database
# cp -f "${TMP_PATH}/targetDB_mapping" "${TAXDBNAME}_mapping"
# cp -f "${NCBITAXINFO}/names.dmp"     "${TAXDBNAME}_names.dmp"
# cp -f "${NCBITAXINFO}/nodes.dmp"     "${TAXDBNAME}_nodes.dmp"
# cp -f "${NCBITAXINFO}/merged.dmp"    "${TAXDBNAME}_merged.dmp"
# cp -f "${NCBITAXINFO}/delnodes.dmp"  "${TAXDBNAME}_delnodes.dmp"
# echo "Database created"

# if [ -n "$REMOVE_TMP" ]; then
#    echo "Remove temporary files"
#    rm -f "${TMP_PATH}/names.dmp" "${TMP_PATH}/nodes.dmp" "${TMP_PATH}/merged.dmp" "${TMP_PATH}/delnodes.dmp"
#    rm -f "${TMP_PATH}/taxidmapping"
#    if [ "$DOWNLOAD_DATA" -eq "1" ]; then
#       rm -f "${TMP_PATH}/ncbi_download.complete" "${TMP_PATH}/targetDB_mapping.complete"
#    fi
#    rm -f "${TMP_PATH}/targetDB_mapping.complete"
#    rm -f "${TMP_PATH}/targetDB_mapping"
#    rm -f createtaxdb.sh
# fi
