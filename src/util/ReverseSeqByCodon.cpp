#include "LocalParameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FileUtil.h"
#include <string>

#ifdef OPENMP
#include <omp.h>
#endif

int reverseseqbycodon(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> seqReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    seqReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter revSeqWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, seqReader.getDbtype());
    revSeqWriter.open();
    Debug::Progress progress(seqReader.getSize());

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::string revStr;
        revStr.reserve(32000);

#pragma omp for schedule(dynamic, 100)
        for (size_t id = 0; id < seqReader.getSize(); id++) {
            progress.updateProgress();
            unsigned int seqKey = seqReader.getDbKey(id);
            char *seq = seqReader.getData(id, thread_idx);
            size_t lenSeq = seqReader.getSeqLen(id);
            std::string codon;

            for (size_t i = 0; i < lenSeq; ++i) {
                size_t revInd = lenSeq - i - 1;
                codon.push_back(seq[revInd]);
                if (codon.size() == 3){
                    std::reverse(codon.begin(), codon.end());
                    revStr.append(codon);
                    codon.clear();
                }
            }

            revStr.push_back('\n');
            
            revSeqWriter.writeData(revStr.c_str(), revStr.size(), seqKey, thread_idx, true);
            revStr.clear();
        }
    }
    revSeqWriter.close(true);
    seqReader.close();
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SEQUENCE_ANCILLARY);

    return EXIT_SUCCESS;
}
