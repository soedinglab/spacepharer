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

abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

hasCommand wget

MMSEQS="$(abspath "${MMSEQS}")"
OUTDB="$(abspath "$2")"
TMP_PATH="$(abspath "$3")"
if [ -n "${GENOME_FTP}" ]; then
    GENOME_FTP="$(abspath "${GENOME_FTP}")"
fi

cd "${TMP_PATH}"
if [ -n "${GENOME_FTP}" ]; then
    if notExists downloaded.tsv; then
        : > downloaded.tsv.tmp
        : > downloaded.tax.tmp
        TAB="$(printf '\t')"
        while IFS="${TAB}" read -r URL TAXON; do
            NAME="$(basename "${URL}")"
            if wget -N -np -nv "${URL}"; then
                push_back "${NAME}"
                echo "${NAME}" >> downloaded.tsv.tmp
            fi
            if [ -n "${TAXON}" ]; then
                printf "%s\t%s\n" "${NAME}" "${TAXON}" >> downloaded.tax.tmp
            else
                INVALID_TAXON=1
            fi
        done < "${GENOME_FTP}"
        mv -f downloaded.tsv.tmp downloaded.tsv
        if [ -z "${INVALID_TAXON}" ]; then
            mv -f downloaded.tax.tmp downloaded.tax
        fi
    else
        while read -r NAME; do
            push_back "${NAME}"
        done < downloaded.tsv
    fi
    eval "set -- $ARR"
    if notExists "${OUTDB}.index"; then
        if [ -f "downloaded.tax" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" createsetdb "${@}" "${OUTDB}" "${TMP_PATH}" --reverse-fragments 0 --tax-mapping-file "downloaded.tax" ${CREATESETDB_PAR} \
                || fail "createsetdb failed"
        else
            # shellcheck disable=SC2086
            "${MMSEQS}" createsetdb "${@}" "${OUTDB}" "${TMP_PATH}" --reverse-fragments 0 ${CREATESETDB_PAR} \
                || fail "createsetdb failed"
        fi
    fi

    if [ -n "${CREATE_REVERSE_SETDB}" ] && notExists "${OUTDB}_rev.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createsetdb "${@}" "${OUTDB}_rev" "${TMP_PATH}" --reverse-fragments 1 ${CREATESETDB_PAR} \
            || fail "create reverse setdb failed"
    fi
else
    case "${DBNAME}" in
      "GenBank_phage_2018_09")
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/genbank_phages_2018_09.tar";
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/genbank_phages_2018_09.tsv";
          IN_TAR="genbank_phages_2018_09.tar"
          IN_TAX="genbank_phages_2018_09.tsv"
      ;;
      "GenBank_eukvir_2018_09")
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/genbank_eukvir_2018_09.tar";
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/genbank_eukvir_2018_09.tsv";
          IN_TAR="genbank_eukvir_2018_09.tar"
          IN_TAX="genbank_eukvir_2018_09.tsv"
      ;;
      "spacers_shmakov_et_al_2017")
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/spacers_shmakov_et_al_2017.tar.gz";
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/spacers_shmakov_et_al_2017.tsv";
          IN_TAR="spacers_shmakov_et_al_2017.tar.gz"
          IN_TAX="spacers_shmakov_et_al_2017.tsv"
          unset CREATE_REVERSE_SETDB
      ;;
      "spacers_dion_et_al_2021")
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/spacers_dion_et_al_2021.tar.gz";
          wget -N -np -nv "http://wwwuser.gwdg.de/~compbiol/spacepharer/2018_09/spacers_dion_et_al_2021.tsv";
          IN_TAR="spacers_dion_et_al_2021.tar.gz"
          IN_TAX="spacers_dion_et_al_2021.tsv"
          unset CREATE_REVERSE_SETDB
      ;;
    esac

    if notExists "${TMP_PATH}/tardb.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" tar2db "${IN_TAR}" "${TMP_PATH}/tardb" ${THREADS_PAR} \
            || fail "tar2db failed"
    fi

    if notExists "${TMP_PATH}/seqdb.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createdb "${TMP_PATH}/tardb" "${TMP_PATH}/seqdb" ${VERBOSITY_PAR} \
            || fail "createdb failed"
    fi

    if notExists "${OUTDB}.index"; then
        if [ -n "${IN_TAX}" ]; then
            # shellcheck disable=SC2086
            "${MMSEQS}" createsetdb "${TMP_PATH}/seqdb" "${OUTDB}" "${TMP_PATH}" --reverse-fragments 0 --tax-mapping-file "${IN_TAX}" ${CREATESETDB_PAR} \
                || fail "createsetdb failed"
        else
            # shellcheck disable=SC2086
            "${MMSEQS}" createsetdb "${TMP_PATH}/seqdb" "${OUTDB}" "${TMP_PATH}" --reverse-fragments 0 ${CREATESETDB_PAR} \
                || fail "createsetdb failed"
        fi
    fi

    if [ -n "${CREATE_REVERSE_SETDB}" ] && notExists "${OUTDB}_rev.index"; then
        # shellcheck disable=SC2086
        "${MMSEQS}" createsetdb "${TMP_PATH}/seqdb" "${OUTDB}_rev" "${TMP_PATH}" --reverse-fragments 1 ${CREATESETDB_PAR} \
            || fail "create reverse setdb failed"
    fi
fi
