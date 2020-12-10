#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "LocalParameters.h"
#include "predictmatch.sh.h"

void setpredictmatchDefaults(Parameters *p) {
    //multihitsearch par
    p->sensitivity = 5.7;
    p->orfMinLength = 9;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV; 
    p->evalThr = 200;
    p->maxSequences = 1500;
    p->addBacktrace = 1;
    p->kmerSize = 6;
    p->spacedKmer = true;
    p->spacedKmerPattern = "11011101";
    //TODO: change path for VTML40
    p->scoringMatrixFile = MultiParam<char*>("VTML40.out", "nucleotide.out");
    p->gapExtend = MultiParam<int>(2, 2);
    p->gapOpen = MultiParam<int>(16, 5);

    //besthitpar
    p->simpleBestHit = true;
    //combinepval par
    p->aggregationMode = 3;
}

int predictmatch(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    setpredictmatchDefaults(&par);

    par.PARAM_MAX_REJECTED.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_OVERLAP.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.extractorfs.size(); i++){
        par.extractorfs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.translatenucs.size(); i++){
        par.translatenucs[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    for (size_t i = 0; i < par.result2profile.size(); i++){
        par.result2profile[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.parseParameters(argc, argv, command, true, 0, 0);

    // check if temp dir exists and if not, try to create it:
    std::string tmpDir = par.db5;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.predictmatchworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("REVERSE_FRAGMENTS", par.reverseFragments == 1 ? "TRUE" : NULL);
    cmd.addVariable("REPORT_PAM", par.reportPam == 1 ? "TRUE" : NULL);
    cmd.addVariable("PERFORM_NUCLALN", par.performNuclAln == 1 ? "TRUE" : NULL);
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("SWAPDB_PAR", par.createParameterString(par.swapdb).c_str());
    par.stat = "linecount";
    cmd.addVariable("RESULT2STATS_PAR", par.createParameterString(par.result2stats).c_str());
    cmd.addVariable("SEARCH_PAR", par.createParameterString(par.searchworkflow).c_str());
    cmd.addVariable("PROTALN2NUCL_PAR", par.createParameterString(par.proteinaln2nucl).c_str());
    cmd.addVariable("BESTHITBYSET_PAR", par.createParameterString(par.besthitbyset).c_str());
    cmd.addVariable("COMBINEPVALPERSET_PAR", par.createParameterString(par.combinepvalbyset).c_str());
    cmd.addVariable("FILTERMATCHBYFDR_PAR", par.createParameterString(par.filtermatchbyfdr).c_str());
    cmd.addVariable("FINDPAM_PAR", par.createParameterString(par.findpam).c_str());
    cmd.addVariable("SUMMARIZERESULTS_PAR", par.createParameterString(par.summarizeresults).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());

    FileUtil::writeFile(tmpDir + "/predictmatch.sh", predictmatch_sh, predictmatch_sh_len);
    std::string program(tmpDir + "/predictmatch.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
