include_directories(commons)
add_subdirectory(commons)
add_subdirectory(util)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(spacepharer
        ${commons_source_files}
        ${util_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        spacepharer.cpp
        )
mmseqs_setup_derived_target(spacepharer)
target_link_libraries(spacepharer version)

add_dependencies(spacepharer local-generated)

install(TARGETS spacepharer DESTINATION bin)
