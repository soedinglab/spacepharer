#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include <cassert>

#include "downloadgenome.sh.h"

struct GenomesDownload {
    const char *name;
    const char *description;
    const char *citation;
};

std::vector<GenomesDownload> genomesDownloads = {
    {
       "GenBank_phage_2018_09",
       "GenBank phage genomes from September 2018",
       "NCBI Resource Coordinators: Database resources of the National Center for Biotechnology Information. Nucleic Acids Res 46(D1), D8-D13 (2018)",

    },
    {
       "GenBank_eukvir_2018_09",
       "GenBank eukaryotic viral genomes from September 2018",
       "NCBI Resource Coordinators: Database resources of the National Center for Biotechnology Information. Nucleic Acids Res 46(D1), D8-D13 (2018)",
    }
};

extern void appendPadded(std::string& dst, const std::string& value, size_t n, int direction, char padding);

std::string listGenomeDatabases(const Command &command, bool detailed) {
    size_t nameWidth = 4;
    for (size_t i = 0; i < genomesDownloads.size(); ++i) {
        nameWidth = std::max(nameWidth, strlen(genomesDownloads[i].name));
    }

    std::string description;
    description.reserve(1024);
    if (detailed) {
        description += " By ";
        description += command.author;
        description += "\n";
    }

    description += "\n  ";
    appendPadded(description, "Name", nameWidth, 0, ' ');
    description.append(1, '\n');

    for (size_t i = 0; i < genomesDownloads.size(); ++i) {
        description.append("- ");
        appendPadded(description, genomesDownloads[i].name, nameWidth, 0, ' ');
        description.append(1, '\n');
        if (strlen(genomesDownloads[i].description) > 0) {
            description.append(2, ' ');
            description.append(genomesDownloads[i].description);
            description.append(1, '\n');
        }
        if (strlen(genomesDownloads[i].citation) > 0) {
            description.append("  Cite: ");
            description.append(genomesDownloads[i].citation);
            description.append(1, '\n');
        }
    }
    description.append(1, '\n');

    return description;
}

int downloadgenome(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    std::string description = listGenomeDatabases(command, par.help);
    if (par.filenames.size() == 0 || par.help) {
        par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
        EXIT(EXIT_SUCCESS);
    }

    par.printParameters(command.cmd, argc, argv, par.downloadgenome);
    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.downloadgenome));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);
    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    ssize_t downloadIdx = -1;
    for (size_t i = 0; i < genomesDownloads.size(); ++i) {
        if (par.db1 == std::string(genomesDownloads[i].name)) {
            downloadIdx = i;
            break;
        }
    }

    CommandCaller cmd;
    if (downloadIdx == -1) {
        if (FileUtil::fileExists(par.db1.c_str()) == false) {
            par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
            Debug(Debug::ERROR) << "Selected database " << par.db1 << " was not found\n";
            EXIT(EXIT_FAILURE);
        }
        cmd.addVariable("GENOME_FTP", par.db1.c_str());
        cmd.addVariable("DBNAME", NULL);
    } else {
        cmd.addVariable("GENOME_FTP", NULL);
        cmd.addVariable("DBNAME", genomesDownloads[downloadIdx].name);
    }
    cmd.addVariable("CREATE_REVERSE_SETDB", par.reverseSetDb == 1 ? "TRUE" : NULL);
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program(tmpDir + "/downloadgenome.sh");
    FileUtil::writeFile(program, downloadgenome_sh, downloadgenome_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    EXIT(EXIT_FAILURE);
}
