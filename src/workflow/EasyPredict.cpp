#include "FileUtil.h"
#include "CommandCaller.h"
#include "Util.h"
#include "Debug.h"
#include "LocalParameters.h"
#include "easypredict.sh.h"
#include <cassert>

extern void setpredictmatchDefaults(Parameters *p);

int easypredict(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.PARAM_DB_OUTPUT.addCategory(MMseqsParameter::COMMAND_EXPERT);
    for (size_t i = 0; i < par.searchworkflow.size(); i++) {
        par.searchworkflow[i]->addCategory(MMseqsParameter::COMMAND_EXPERT);
    }
    par.PARAM_COMPRESSED.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_THREADS.removeCategory(MMseqsParameter::COMMAND_EXPERT);
    par.PARAM_V.removeCategory(MMseqsParameter::COMMAND_EXPERT);

    setpredictmatchDefaults(&par);
    par.dbOut = false;
    par.PARAM_DB_OUTPUT.wasSet = true;
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string tmpDir = par.filenames.back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, *command.params));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("RESULTS", par.filenames.back().c_str());
    par.filenames.pop_back();
    std::string target = par.filenames.back().c_str();
    cmd.addVariable("TARGET", target.c_str());
    par.filenames.pop_back();

    cmd.addVariable("PREDICTMATCH_PAR", par.createParameterString(par.predictmatchworkflow, true).c_str());
    par.extractorfsSpacer = 1;
    cmd.addVariable("CREATESETDB_PAR", par.createParameterString(par.createsetdbworkflow).c_str());
    cmd.addVariable("THREADS_COMP_PAR", par.createParameterString(par.threadsandcompression).c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY", par.createParameterString(par.onlyverbosity).c_str());

    std::string program = tmpDir + "/easypredict.sh";
    FileUtil::writeFile(program, easypredict_sh, easypredict_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    return EXIT_FAILURE;
}

