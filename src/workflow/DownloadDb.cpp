#include "Util.h"
#include "LocalParameters.h"
#include "Debug.h"
#include "FileUtil.h"
#include "CommandCaller.h"

#include <cassert>

#include "downloaddb.sh.h"

struct SpacePharerDownload {
    const char *name;
    const char *description;
    const char *citation;
    enum DbType {
        SPACER,
        GENOME
    } type;

    static const char* formatType(DbType type) {
        switch (type) {
            case SPACER:
                return "Spacer";
            case GENOME:
                return "Genome";
            default:
                return "-";
        }
    }
};

std::vector<SpacePharerDownload> genomesDownloads = {
    {
        "GenBank_phage_2018_09",
        "GenBank phage genomes from September 2018",
        "NCBI Resource Coordinators: Database resources of the National Center for Biotechnology Information. Nucleic Acids Res 46(D1), D8-D13 (2018)",
            SpacePharerDownload::GENOME,
    },
    {
        "GenBank_eukvir_2018_09",
        "GenBank eukaryotic viral genomes from September 2018",
        "NCBI Resource Coordinators: Database resources of the National Center for Biotechnology Information. Nucleic Acids Res 46(D1), D8-D13 (2018)",
            SpacePharerDownload::GENOME,
    },
    {
        "spacers_shmakov_et_al_2017",
        "Spacers extracted from Shmakov et al",
        "The CRISPR Spacer Space Is Dominated by Sequences from Species-Specific Mobilomes. mBio 8(5), e01397-17 (2017)",
            SpacePharerDownload::SPACER,
    }
};

extern void appendPadded(std::string& dst, const std::string& value, size_t n, int direction, char padding);

std::string listGenomeDatabases(const Command &command, bool detailed) {
    size_t nameWidth = 4;
    size_t typeWidth = 4;
    for (size_t i = 0; i < genomesDownloads.size(); ++i) {
        nameWidth = std::max(nameWidth, strlen(genomesDownloads[i].name));
        nameWidth = std::max(nameWidth, strlen(SpacePharerDownload::formatType(genomesDownloads[i].type)));
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
    description.append(1, '\t');
    appendPadded(description, "Type", typeWidth, 0, ' ');
    description.append(1, '\n');

    for (size_t i = 0; i < genomesDownloads.size(); ++i) {
        description.append("- ");
        appendPadded(description, genomesDownloads[i].name, nameWidth, 0, ' ');
        description.append(1, '\t');
        appendPadded(description, SpacePharerDownload::formatType(genomesDownloads[i].type), typeWidth, 0, ' ');
        description.append(1, '\n');
        if (detailed == false) {
            continue;
        }
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


extern const char* binary_name;

int downloaddb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    std::string description = listGenomeDatabases(command, par.help);
    if (par.filenames.size() == 0 || par.help) {
        if (par.help == false) {
            description.append("Show a detailed list of databases by calling '");
            description.append(binary_name);
            description.append(1, ' ');
            description.append(command.cmd);
            description.append(" -h'\n\n");
        }
        par.printUsageMessage(command, par.help ? MMseqsParameter::COMMAND_EXPERT : 0, description.c_str());
        EXIT(EXIT_SUCCESS);
    }

    par.printParameters(command.cmd, argc, argv, *command.params);
    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.downloaddb));
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
        if (genomesDownloads[downloadIdx].type == SpacePharerDownload::SPACER) {
            par.extractorfsSpacer = 1;
            if (par.PARAM_REVERSE_SETDB.wasSet == false) {
                par.reverseSetDb = 0;
            }
        }
    }
    std::vector<MMseqsParameter*> createsetdbWithoutRevTax;
    for (size_t i = 0; i < par.createsetdbworkflow.size(); i++) {
        if (par.createsetdbworkflow[i]->uniqid != par.PARAM_REVERSE_FRAGMENTS.uniqid && par.createsetdbworkflow[i]->uniqid != par.PARAM_TAX_MAPPING_FILE.uniqid) {
            createsetdbWithoutRevTax.push_back(par.createsetdbworkflow[i]);
        }
    }

    cmd.addVariable("CREATESETDB_PAR", par.createParameterString(createsetdbWithoutRevTax).c_str());
    cmd.addVariable("CREATE_REVERSE_SETDB", par.reverseSetDb == 1 ? "TRUE" : NULL);
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());

    std::string program(tmpDir + "/downloadb.sh");
    FileUtil::writeFile(program, downloaddb_sh, downloaddb_sh_len);
    cmd.execProgram(program.c_str(), par.filenames);

    // Should never get here
    assert(false);
    EXIT(EXIT_FAILURE);
}
