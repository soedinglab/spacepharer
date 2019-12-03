#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "LocalParameters.h"
#include "downloadgenome.sh.h"


int downloadgenome(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0 , 0);

    // check if temp dir exists and if not, try to create it:
    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.onlyverbosity));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    CommandCaller cmd;
    cmd.addVariable("CREATE_REVERSE_SETDB", par.reverseSetDb == 1 ? "TRUE" : NULL);
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program(tmpDir + "/downloadgenome.sh");
    FileUtil::writeFile(program.c_str(), downloadgenome_sh, downloadgenome_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
