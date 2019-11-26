#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const int NO_CITATION = 0;
const char* binary_name = "spacepharer";
const char* tool_name = "spacepharer";
const char* tool_introduction = "CRISPR spacer phage-host finder";
const char* main_author = "Ruoshi Zhang, ruoshi.zhang@mpibpc.mpg.de";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<Command> commands = {
        {"parsespacer",             parsespacer,            &localPar.parsespacer,    COMMAND_MAIN,
                "Parse spacer files and create sequence database",
                "Parse output files from PILER-CR, CRISPRDetect and CRT (CRISPRFinder still under construction)",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:spacerFile1[.txt]> ... <i:spacerFileN[.txt]> <o:spacerDB>",
                NO_CITATION, {{"spacerFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
                                    {"spacerDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}}},
        {"downloadgenome",             downloadgenome,            &localPar.downloadgenomeworkflow,    COMMAND_MAIN,
                "Download GenBank phage genomes",
                "The program will download a default list of phage genomes from NCBI GenBank ftp in provided directory, and create set db and optionally reversed set db",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i: genomeFtpFile> <i:outDir> <o:setDB> <tmpDir>",
                NO_CITATION, {{"genomeFtpFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA| DbType::VARIADIC, &DbValidator::flatfile},
                                {"outDir", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::directory},
                                {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

       {"createsetdb",             createsetdb,            &localPar.createsetdbworkflow,    COMMAND_MAIN,
                "Create sequence database from sets",
                "Create sequence database and associated metadata for predictmatch",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]>|<i:DB> <o:setDB> <tmpDir>",
                NO_CITATION, {{"fast[a|q]File[.gz|bz]",  DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA | DbType::VARIADIC,  &DbValidator::flatfile},
                                   {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                   {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},

        {"predictmatch",             predictmatch,            &localPar.predictmatchworkflow,    COMMAND_MAIN,
                "Predict phage-host matches",
                "Create query sequence database and search with a query set of sequences against target set",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:queryDB> <i:targetDB> <i:controltargetDB> <o:cEvalDB> <tmpDir>",
                NO_CITATION, {
                        //{"queryfast[a|q]File[.gz|bz]",  DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA | DbType::VARIADIC,  &DbValidator::flatfile },
                                         {"queryDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                         {"targetDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                         {"controltargetDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                         {"cEvalDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                         {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
                                         
       {"filtermatchbyfdr",             filtermatchbyfdr,            &localPar.filtermatchbyfdr,    COMMAND_MAIN,
                "Report matches based on FDR",
                "Report matches based on FDR(false discovey rate)",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:targetDB> <i:cEvalDB> <i:controltargetDB> <i:controlcEvalDB> <o:resultDB>",
                NO_CITATION, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"cEvalDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                   {"controltargetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"controlcEvalDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb}}},
       {"findpam",             findpam,            &localPar.findpam,    COMMAND_MAIN,
                "Report PAM upstream or downstream of protospacer",
                "Report if finds PAM(protospacer adjacent motif) upstream or downstream of protospacer and append to alignment result",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:targetDB> <i:alnResultDB>  <o:resultDB>",
                NO_CITATION, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                //convertalis produces genericDB
                                   {"alnResultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb}}},
       {"truncatebesthits",             truncatebesthits,            &localPar.truncatebesthits,    COMMAND_MAIN,
                "Truncate list of besthits",
                "Truncate list of besthits based on query set size",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:queryDB> <i:besthitsDB> <o:resultDB>",
                NO_CITATION, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                   {"besthitsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb}}},
       {"summarizeresults",             summarizeresults,            &localPar.summarizeresults,    COMMAND_MAIN,
                "Summarize results",
                "Summarize results on predicted matches (E-value) and hits (alignments and PAM)",
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:matchDB> <i:alnDB> <o:resultDB>",
                NO_CITATION, {{"matchDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                   {"alnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}, //&DbValidator::resultDb},
                                   {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb}}}


};

