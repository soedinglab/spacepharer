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

GENOME_FTP="$1"
OUT_PATH="$2"
OUTDB="$3"
TMP_PATH="$4"

#change path
INPUTS=""
if notExists "${OUT_PATH}/phage_download.complete"; then
    echo "Download phage genomes"
    while read NAME URL; 
    do
        wget -nv -O "${OUT_PATH}/${NAME}" "${URL}" ;        
        INPUTS="${OUT_PATH}/${NAME} ${INPUTS}";
    done <  "${GENOME_FTP}"
    touch "${OUT_PATH}/phage_download.complete" \
        || fail "download failed"
    
else
    while read NAME URL;
    do
        INPUTS="${OUT_PATH}/${NAME} ${INPUTS}";
    done < "${GENOME_FTP}"
fi


if notExists "${OUTDB}.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsetdb ${INPUTS} "${OUTDB}" "${TMP_PATH}" \
        || fail "createsetdb failed"
fi

if [ -n "$CREATE_REVERSE_SETDB" ] && notExists "${OUTDB}_rev.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsetdb ${INPUTS} "${OUTDB}_rev" "${TMP_PATH}" "--reverse-fragments" "1" \
        || fail "create reverse setdb failed"
fi
