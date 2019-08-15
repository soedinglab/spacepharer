include_directories(commons)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(whisper
        ${commons_source_files}
        ${workflow_source_files}
        LocalCommandDeclarations.h
        whisper.cpp
        )
mmseqs_setup_derived_target(whisper)
target_link_libraries(whisper version)

add_dependencies(whisper local-generated)

install(TARGETS whisper DESTINATION bin)
