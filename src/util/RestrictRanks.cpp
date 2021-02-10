#include "LocalParameters.h"
#include "Debug.h"
#include "NcbiTaxonomy.h"
#include "DBReader.h"
#include "DBWriter.h"

#ifdef OPENMP
#include <omp.h>
#endif

const char *maxRank(double seqId, const double thresholds[8]) {
    if (seqId > thresholds[0]) {
        return "species";
    } else if (seqId > thresholds[1]) {
        return "genus";
    } else if (seqId > thresholds[2]) {
        return "family";
    } else if (seqId > thresholds[3]) {
        return "order";
    } else if (seqId > thresholds[4]) {
        return "class";
    } else if (seqId > thresholds[5]) {
        return "phylum";
    } else if (seqId > thresholds[6]) {
        return "kingdom";
    } else if (seqId > thresholds[7]) {
        return "superkingdom";
    } else {
        return NULL;
    }
}

int restrictranks(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy *t = NcbiTaxonomy::openTaxonomy(par.db1.c_str());
    std::vector<std::string> ranks = NcbiTaxonomy::parseRanks(par.lcaRanks);
    // will be used when no hits
    std::string noTaxResult = "0\tno rank\tunclassified";
    if (!ranks.empty()) {
        noTaxResult += '\t';
    }
    if (par.showTaxLineage > 0) {
        noTaxResult += '\t';
    }
    noTaxResult += '\n';

    std::vector<std::string> splitThresholds = Util::split(par.rankMinSeqIds, ",");
    if (splitThresholds.size() != 8) {
        Debug(Debug::ERROR) << "Expected 8 sequence identity thresholds, got " <<  splitThresholds.size() << " instead\n";
        EXIT(EXIT_FAILURE);
    }
    double thresholds[8];
    for (size_t i = 0; i < 8; ++i) {
        thresholds[i] = strtod(splitThresholds[i].c_str(), NULL);
    }

    DBReader<unsigned int> taxReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    taxReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> matchReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    matchReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    Debug::Progress progress;
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        const char *entry[255];
        std::string buffer;
        buffer.reserve(2048);
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < taxReader.getSize(); ++i) {
            progress.updateProgress();

            unsigned int key = taxReader.getDbKey(i);

            size_t matchId = matchReader.getId(key);
            char *data = matchReader.getData(matchId, thread_idx);
            double seqId = 0;
            size_t count = 0;
            while (*data != '\0') {
                Util::getWordsOfLine(data, entry, 255);
                data = Util::skipLine(data);
                // 0% seqid hits are an artifact by filtermatchbyfdr when the control database is too small
                // ignore them here, if we ignore them in filtermatchbyfdr we cant annotate many sequences as unassigned
                // as they disappear completely due to how swapdb works
                double currSeqId = strtod(entry[2], NULL);
                if (currSeqId > 0.0f) {
                    seqId += currSeqId;
                    count++;
                }
            }
            seqId = seqId / count;

            data = taxReader.getData(i, thread_idx);
            writer.writeStart(thread_idx);
            while (*data != '\0') {
                Util::getWordsOfLine(data, entry, 255);
                data = Util::skipLine(data);

                TaxID taxon = Util::fast_atoi<unsigned int>(entry[0]);

                const char *bestRank = maxRank(seqId, thresholds);
                if (bestRank == NULL) {
                    writer.writeAdd(noTaxResult.c_str(), noTaxResult.length(), thread_idx);
                    continue;
                }

                const TaxonNode* node = t->taxonNode(taxon, false);
                if (node == NULL) {
                    writer.writeAdd(noTaxResult.c_str(), noTaxResult.length(), thread_idx);
                    continue;
                }

                const char* rank = t->getString(node->rankIdx);
                int rankLevel = t->findRankIndex(rank);
                int bestLevel = t->findRankIndex(bestRank);
                if (rankLevel >= bestLevel) {
                    buffer.append(SSTR(taxon));
                    buffer.append(1, '\t');
                    buffer.append(rank);
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
                    buffer.append(1, '\n');
                    writer.writeAdd(buffer.c_str(), buffer.length(), thread_idx);
                    buffer.clear();
                    continue;
                }

                while (node->parentTaxId != node->taxId) {
                    node = t->taxonNode(node->parentTaxId, false);
                    std::string currRank = t->getString(node->rankIdx);
                    if (currRank != "no rank" && t->findRankIndex(currRank) >= bestLevel) {
                        break;
                    }
                }
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
                buffer.append(1, '\n');
                writer.writeAdd(buffer.c_str(), buffer.length(), thread_idx);
                buffer.clear();
            }
            writer.writeEnd(key, thread_idx);
        }
    }
    writer.close();
    taxReader.close();
    matchReader.close();
    delete t;
    return EXIT_SUCCESS;
}
