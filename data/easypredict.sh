#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

if notExists "${TMP_PATH}/qdb.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" parsespacer "$@" "${TMP_PATH}/qdb" ${THREADS_COMP_PAR} \
        || fail "parsespacer failed"
fi

if notExists "${TMP_PATH}/qsetdb.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createsetdb "${TMP_PATH}/qdb" "${TMP_PATH}/qsetdb" "${TMP_PATH}/createsetdb" --extractorf-spacer 1 --reverse-fragments 1 ${CREATESETDB_PAR} \
        || fail "createsetdb failed"
fi

# shellcheck disable=SC2086
"$MMSEQS" predictmatch "${TMP_PATH}/qsetdb" "${TARGET}" "${TARGET}_rev" "${RESULTS}" "${TMP_PATH}/predict" ${PREDICTMATCH_PAR} \
    || fail "predictmatch failed"

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/qdb" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$MMSEQS" rmdb "${TMP_PATH}/qsetdb" ${VERBOSITY}
    rm -f "${TMP_PATH}/easypredict.sh"
fi
