include(MMseqsResourceCompiler)

set(COMPILED_RESOURCES
        easypredict.sh
        predictmatch.sh
        createsetdb.sh
        downloadgenome.sh
        genbank_phages_2018_09.tsv
        genbank_eukviruses_2018_09.tsv
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(local-generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
