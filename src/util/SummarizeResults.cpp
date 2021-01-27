#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "NcbiTaxonomy.h"

#ifdef OPENMP
#include <omp.h>
#endif

int summarizeresults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy *t = NULL;
    std::string nodesDb = par.db1 + "_nodes.dmp";
    if (FileUtil::fileExists(nodesDb.c_str())) {
        t = NcbiTaxonomy::openTaxonomy(par.db1);
    }
    std::vector<std::string> ranks = NcbiTaxonomy::parseRanks(par.lcaRanks);


    DBReader<unsigned int> matchReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    matchReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> alnReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::NOSORT);

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), par.threads, shouldCompress, dbType);
    dbw.open();

    Debug::Progress progress(matchReader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        const char *entry[255];

        std::string cScore;
        cScore.reserve(255);

        std::string buffer;
        buffer.reserve(1024 * 1024);

        std::string tmpBuffer;
        tmpBuffer.reserve(1024 * 1024);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < matchReader.getSize(); ++id) {
            progress.updateProgress();

            unsigned int matchKey = matchReader.getDbKey(id);
            size_t alnId = alnReader.getId(matchKey);
            if (alnId == UINT_MAX) {
                Debug(Debug::WARNING) << "Missing alignment result for key " << matchKey << " \n";
                continue;
            }
            char *alnData = alnReader.getData(alnId, thread_idx);
            char *matchData = matchReader.getData(id, thread_idx);
            while (*matchData != '\0') {
                size_t columns = Util::getWordsOfLine(matchData, entry, 255);
                matchData = Util::skipLine(matchData);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid alignment result record\n";
                    EXIT(EXIT_FAILURE);
                }
                unsigned int targetSetKey = Util::fast_atoi<unsigned int>(entry[0]);

                TaxID taxId = 0;
                size_t lineCount = 0;
                char* alnCurrent = alnData;
                while (*alnCurrent != '\0') {
                    // read only key
                    unsigned int dbKey = Util::fast_atoi<unsigned int>(alnCurrent);
                    if (dbKey != targetSetKey) {
                        alnCurrent = Util::skipLine(alnCurrent);
                        continue;
                    }

                    // call before getWordsOfLine or entry is overwritten
                    if (lineCount == 0) {
                        cScore.assign(entry[1], entry[2] - entry[1]);
                    }

                    columns = Util::getWordsOfLine(alnCurrent, entry, 255);
                    if (columns < 10) {
                        Debug(Debug::ERROR) << "Invalid alignment result record\n";
                        EXIT(EXIT_FAILURE);
                    }
                    alnCurrent = Util::skipLine(alnCurrent);

                    if (lineCount == 0) {
                        buffer.append("#");
                        buffer.append(entry[2], entry[3] - entry[2] - 1);
                        buffer.append("\t");
                        buffer.append(entry[3], entry[4] - entry[3] - 1);
                        buffer.append("\t");
                        buffer.append(cScore);
                        //buffer.append("\t");
                        if (t != NULL) {
                            taxId = Util::fast_atoi<TaxID>(entry[12]);
                        }
                    }
                    lineCount++;

                    if (par.formatType != LocalParameters::FORMAT_TYPE_LONG && par.formatType != LocalParameters::FORMAT_TYPE_ALN) {
                        continue;
                    }
                    tmpBuffer.append(">");
                    tmpBuffer.append(entry[1], entry[2] - entry[1] - 1); //spacername
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[3], entry[4] - entry[3] - 1); //genomename
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[4], entry[5] - entry[4] - 1); //best-hit pval
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[5], entry[6] - entry[5] - 1); //qstart
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[6], entry[7] - entry[6] - 1); //qend
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[8], entry[9] - entry[8] - 1); //tstart
                    tmpBuffer.append("\t");
                    tmpBuffer.append(entry[9], entry[10] - entry[9] - 1); //tend
                    if (columns == 15) {
                        tmpBuffer.append("\t");
                        tmpBuffer.append(entry[13], entry[14] - entry[13] - 1); //PAM
                        tmpBuffer.append("\t");
                        tmpBuffer.append(entry[14], entry[15] - entry[14]); //PAM reverse complement
                    }
                    if (par.formatType == LocalParameters::FORMAT_TYPE_ALN) {
                        tmpBuffer.append("\n");
                        tmpBuffer.append(entry[10], entry[11] - entry[10] - 1); //qaln
                        tmpBuffer.append("\n");
                        tmpBuffer.append(entry[11], entry[12] - entry[11] - 1); //taln
                    }
                    tmpBuffer.append("\n");
                }
                if (lineCount > 0) {
                    buffer.append(SSTR(lineCount));
                    if (t != NULL) {
                        TaxonNode const * node = t->taxonNode(taxId, false);
                        if (node != NULL) {
                            buffer.append(1, '\t');
                            buffer.append(SSTR(node->taxId));
                            buffer.append(1, '\t');
                            buffer.append(t->getString(node->rankIdx));
                            buffer.append(1, '\t');
                            buffer.append(t->getString(node->nameIdx));
                            if (!ranks.empty()) {
                                buffer.append(1, '\t');
                                buffer.append(Util::implode(t->AtRanks(node, ranks), ';'));
                            }
                            if (par.showTaxLineage == 1) {
                                buffer.append(1, '\t');
                                buffer.append(t->taxLineage(node, true));
                            }
                            if (par.showTaxLineage == 2) {
                                buffer.append(1, '\t');
                                buffer.append(t->taxLineage(node, false));
                            }
                        }
                    }
                    buffer.append(1, '\n');
                    buffer.append(tmpBuffer);
                }
                tmpBuffer.clear();
            }
            dbw.writeData(buffer.c_str(), buffer.length(), matchKey, thread_idx, par.dbOut);
            buffer.clear();
        }
    }
    dbw.close(par.dbOut == false);
    if (par.dbOut == false) {
        FileUtil::remove(par.db4Index.c_str());
    }
    matchReader.close();
    alnReader.close();

    return EXIT_SUCCESS;
}

