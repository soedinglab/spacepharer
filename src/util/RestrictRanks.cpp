#include "Debug.h"
#include "NcbiTaxonomy.h"
#include "DBReader.h"
#include "DBWriter.h"

#ifdef OPENMP
#include <omp.h>
#endif

const char *maxRank(double seqId) {
    if (seqId > 0.87) {
        return "species";
    } else if (seqId > 0.83) {
        return "genus";
    } else if (seqId > 0.78) {
        return "family";
    } else if (seqId > 0.70) {
        return "order";
    } else if (seqId > 0.60) {
        return "class";
    } else if (seqId > 0.50) {
        return "phylum";
    } else if (seqId > 0.40) {
        return "kingdom";
    } else if (seqId > 0.30) {
        return "superkingdom";
    } else {
        return NULL;
    }
}

int restrictranks(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy *t = NcbiTaxonomy::openTaxonomy(par.db1.c_str());

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
        char buffer[2048];
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
                seqId += strtod(entry[2], NULL);
                count++;
            }
            seqId = seqId / count;

            data = taxReader.getData(i, thread_idx);
            writer.writeStart(thread_idx);
            while (*data != '\0') {
                Util::getWordsOfLine(data, entry, 255);
                data = Util::skipLine(data);

                TaxID taxon = Util::fast_atoi<unsigned int>(entry[0]);

                const char *bestRank = maxRank(seqId);
                if (bestRank == NULL) {
                    int len = snprintf(buffer, sizeof(buffer), "0\tno rank\tunclassified\t%.3f\n", seqId);
                    writer.writeAdd(buffer, len, thread_idx);
                    continue;
                }

                const TaxonNode* node = t->taxonNode(taxon, false);
                if (node == NULL) {
                    int len = snprintf(buffer, sizeof(buffer), "0\tno rank\tunclassified\t%.3f\n", seqId);
                    writer.writeAdd(buffer, len, thread_idx);
                    continue;
                }

                const char* rank = t->getString(node->rankIdx);
                int rankLevel = t->findRankIndex(rank);
                int bestLevel = t->findRankIndex(bestRank);
                if (rankLevel >= bestLevel) {
                    int len = snprintf(buffer, sizeof(buffer), "%i\t%s\t%s\t%.3f\n", taxon, rank, t->getString(node->nameIdx), seqId);
                    writer.writeAdd(buffer, len, thread_idx);
                    continue;
                }

                while (node->parentTaxId != node->taxId) {
                    node = t->taxonNode(node->parentTaxId, false);
                    std::string currRank = t->getString(node->rankIdx);
                    if (currRank != "no rank" && t->findRankIndex(currRank) >= bestLevel) {
                        break;
                    }
                }
                int len = snprintf(buffer, sizeof(buffer), "%i\t%s\t%s\t%.3f\n", node->taxId, t->getString(node->rankIdx), t->getString(node->nameIdx), seqId);
                writer.writeAdd(buffer, len, thread_idx);
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
