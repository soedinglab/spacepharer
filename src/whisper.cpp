#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const int NO_CITATION = 0;
const char* binary_name = "whisper";
const char* tool_name = "whisper";
const char* tool_introduction = "CRISPR spacer phage-host finder";
const char* main_author = "Ruoshi Zhang, ruoshi.zhang@mpibpc.mpg.de";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<struct Command> commands = {
        //Main tools (workflows for non-experts)
        // {"spacerparser",             spacerparser,            &localPar.spacerparserworkflow,    COMMAND_MAIN,
        //         "Parse spacer files and create sequence database",
        //         "Parse output files from PILER-CR, CRISPRDetect and CRT",
        //         "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>", "Eli Levy Karin <eli.levy.karin@gmail.com>",
        //         "<i:spacerFile1[]> ... <i:spacerFileN[]> <o:setDB> <tmpDir>",
        //         NO_CITATION, {{"spacerFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile ,
        //                            {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
        //                            {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"createsetdb",             createsetdb,            &localPar.createsetdbworkflow,    COMMAND_MAIN,
                "Create sequence database from sets",
                "Create sequence database and associated metadata for predictmatch",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>", "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <o:setDB> <tmpDir>",
                NO_CITATION, {{"fast[a|q]File[.gz|bz]",  DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA | DbType::VARIADIC,  &DbValidator::flatfile },
                                   {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"predictmatch",             predictmatch,            &localPar.predictmatchworkflow,    COMMAND_MAIN,
                "Predict phage-host matches",
                "Create query sequence database and search with a query set of sequences against target set",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>", "Eli Levy Karin <eli.levy.karin@gmail.com>",
                "<i:queryfastaFile1[.gz]> ... <i:queryfastaFileN[.gz]>|<i:queryDB> <i:targetDB> <i:controltargetDB> <o:resultDB> <tmpDir>",
                NO_CITATION, {{"queryfast[a|q]File[.gz|bz]",  DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA | DbType::VARIADIC,  &DbValidator::flatfile }
                                         {"queryDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"targetDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"controltargetDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                         {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb  },
                                         {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}}

};
