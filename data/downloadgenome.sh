#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

ARR=""
push_back() {
    # shellcheck disable=SC1003
    CURR="$(printf '%s' "$1" | awk '{ gsub(/'\''/, "'\''\\'\'''\''"); print; }')"
    if [ -z "$ARR" ]; then
        ARR=''\'$CURR\'''
    else
        ARR=$ARR' '\'$CURR\'''
    fi
}

[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;

hasCommand () {
    command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

hasCommand wget
hasCommand touch

GENOME_FTP="$1"
OUT_PATH="$2"
OUTDB="$3"
TMP_PATH="$4"

#change path
if notExists "${OUT_PATH}/phage_download.complete"; then
    echo "Download phage genomes"
    while read -r NAME URL; do
        wget -nv -O "${OUT_PATH}/${NAME}" "${URL}"
        push_back "${OUT_PATH}/${NAME}"
    done <  "${GENOME_FTP}"
    touch "${OUT_PATH}/phage_download.complete" \
        || fail "download failed"
    
else
    while read -r NAME URL; do
        push_back "${OUT_PATH}/${NAME}"
    done < "${GENOME_FTP}"
fi

eval "set -- $ARR"
if notExists "${OUTDB}.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsetdb "${@}" "${OUTDB}" "${TMP_PATH}" ${VERBOSITY_PAR} \
        || fail "createsetdb failed"
fi

if [ -n "$CREATE_REVERSE_SETDB" ] && notExists "${OUTDB}_rev.index"; then
    # shellcheck disable=SC2086
    "${MMSEQS}" createsetdb "${@}" "${OUTDB}_rev" "${TMP_PATH}" --reverse-fragments 1 ${VERBOSITY_PAR} \
        || fail "create reverse setdb failed"
fi
