#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "spacepharer";
const char* tool_name = "SpacePHARER";
const char* tool_introduction =
"SpacePHARER (CRISPR Spacer Phage-Host Pair Finder) is a sensitive and fast tool for de novo prediction of phage-host relationships via identifying phage genomes that match CRISPR spacers in genomic or metagenomic data.\n\nPlease cite: R. Zhang et al. SpacePHARER: Sensitive identification of phages from CRISPR spacers in prokaryotic hosts. biorxiv, doi:10.1101/2020.05.15.090266 (2020).";
const char* main_author = "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<Command> commands = {
        {"easy-predict",             easypredict,           &localPar.easypredictmatchworkflow, COMMAND_EASY,
                "Predict phage-host matches from common spacer files (PILER-CR, CRISPRDetect and CRT)",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:spacerFile1[.txt]> ... <i:spacerFileN[.txt]> <i:targetDB> <o:output[.tsv]> <tmpDir>",
                CITATION_SPACEPHARER, {{"spacerFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
                                       {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"result", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                       {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"parsespacer",             parsespacer,            &localPar.threadsandcompression,    COMMAND_MAIN,
                "Parse spacer files (PILER-CR, CRISPRDetect and CRT) and create sequence database",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:spacerFile1[.txt]> ... <i:spacerFileN[.txt]> <o:spacerDB>",
                CITATION_SPACEPHARER, {{"spacerFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
                                       {"spacerDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}}},
        {"downloaddb",              downloaddb,             &localPar.downloaddb,               COMMAND_MAIN,
                "Download spacers or phage genomes and create sequence database",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<name/downloadFile> <o:sequenceDB> <tmpDir>",
                CITATION_SPACEPHARER, {{"name", 0, DbType::ZERO_OR_ALL, &DbValidator::empty },
                                       {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                       {"tmpDir",     DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
       {"createsetdb",              createsetdb,            &localPar.createsetdbworkflow,      COMMAND_MAIN,
                "Create sequence database from FASTA input",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]>|<i:DB> <o:setDB> <tmpDir>",
                CITATION_SPACEPHARER, {{"fast[a|q]File[.gz|bz]", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC,  &DbValidator::flatfile},
                                       {"setDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                       {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"predictmatch",            predictmatch,           &localPar.predictmatchworkflow,     COMMAND_MAIN,
                "Predict phage-host matches",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:queryDB> <i:targetDB> <i:controltargetDB> <o:output[.tsv]> <tmpDir>",
                CITATION_SPACEPHARER, {{"queryDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"controlTargetDB",  DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"cScoreDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                                       {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
                                         
       {"filtermatchbyfdr",         filtermatchbyfdr,       &localPar.filtermatchbyfdr,         COMMAND_SPECIAL | COMMAND_EXPERT,
                "Report matches based on false discovey rate",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:cScoreDB> <i:controlcScoreDB> <o:resultDB>",
                CITATION_SPACEPHARER, {{"cScoreDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"controlcScoreDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
       {"findpam",                  findpam,                &localPar.findpam,    COMMAND_SPECIAL | COMMAND_EXPERT,
                "Report PAM upstream or downstream of protospacer",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:targetDB> <i:alnResultDB>  <o:resultDB>",
                CITATION_SPACEPHARER, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"alnResultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                       {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
       {"truncatebesthits",         truncatebesthits,       &localPar.threadsandcompression,    COMMAND_SPECIAL | COMMAND_EXPERT,
                "Truncate list of best hits based on query set size",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:queryDB> <i:besthitsDB> <o:resultDB>",
                CITATION_SPACEPHARER, {{"queryDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"besthitsDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb}}},
       {"summarizeresults",         summarizeresults,       &localPar.summarizeresults,         COMMAND_SPECIAL | COMMAND_EXPERT,
                "Summarize results on predicted matches (E-value) and hits (alignments and PAM)",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:targetDB> <i:matchDB> <i:alnDB> <o:output[.tsv]>",
                CITATION_SPACEPHARER, {{"targetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb},
                                       {"matchDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb},
                                       {"alnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, NULL}, //&DbValidator::resultDb},
                                       {"output", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
       {"combineprotnuclaln",         combineprotnuclaln,       &localPar.combineprotnuclaln,         COMMAND_SPECIAL | COMMAND_EXPERT,
                "Recompute bit score and E-value (protein and nucleotide search)",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:protAlnDB> <i:nuclAlnDB> <o:resultAlnDB>",
                CITATION_SPACEPHARER, {{"protAlnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"nuclAlnDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"resultAlnDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb}}},
       {"empiricalpval",         empiricalpval,       &localPar.empiricalpval,         COMMAND_SPECIAL | COMMAND_EXPERT,
                "Compute an empirical P-value for each query-target set pair given null model dataset",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:cScoreDB> <i:controlcScoreDB> <o:resultDB>",
                CITATION_SPACEPHARER, {{"cScoreDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"controlcScoreDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                       {"resultDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::genericDb}}},
        {"reverseseqbycodon",          reverseseqbycodon,            &localPar.reverseseqbycodon,          COMMAND_SPECIAL | COMMAND_EXPERT,
                "Reverse nucl sequences by codon to generate inverted ORFs",
                NULL,
                "Eli Levy Karin <eli.levy.karin@gmail.com> ",
                "<i:sequenceDB> <o:revSequenceDB>",
                CITATION_MMSEQS2, {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                        {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb }}},
       {"combinescore",         combinescore,       &localPar.combinescore,         COMMAND_SPECIAL | COMMAND_EXPERT,
                "Compute a combined score for each query-target set pair",
                NULL,
                "Ruoshi Zhang <ruoshi.zhang@mpibpc.mpg.de>",
                "<i:querySetDB> <i:targetSetDB> <i:resultDB> <o:pvalDB>",
                CITATION_SPACEPHARER, {{"querySetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"targetSetDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                                           {"resultDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"pvalDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                                                           {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"restrictranks",       restrictranks,      &localPar.restrictranks,        COMMAND_SPECIAL | COMMAND_EXPERT,
                "Restrict taxonomic result to ranks based on sequence identity",
                NULL,
                "Milot Mirdita <milot@mirdita.de>",
                "<i:taxSeqDB> <i:taxResult> <i:matchResult> <o:taxResult>",
                CITATION_SPACEPHARER, {{"taxSeqDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::taxSequenceDb},
                                              {"taxResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::taxResult},
                                              {"matchResult", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::resultDb},
                                              {"taxResult", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::taxResult}}},
};

