#include "Debug.h"
#include "LocalParameters.h"
#include "Aggregation.h"
#include "itoa.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif


class PvalueAggregator : public Aggregation {
public:
    PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                     const std::string &outputDbName, unsigned int threads, unsigned int compressed) :
            Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed) {

        std::string sizeDBName = queryDbName + "_set_size";
        std::string sizeDBIndex = queryDbName + "_set_size.index";
        querySizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        querySizeReader->open(DBReader<unsigned int>::NOSORT);

        sizeDBName = targetDbName + "_set_size";
        sizeDBIndex = targetDbName + "_set_size.index";
        targetSizeReader = new DBReader<unsigned int>(sizeDBName.c_str(), sizeDBIndex.c_str(), threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);
    }

    ~PvalueAggregator() {
        targetSizeReader->close();
        delete targetSizeReader;

        querySizeReader->close();
        delete querySizeReader;
    }

    void prepareInput(unsigned int querySetKey, unsigned int thread_idx) {
        unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
    }

    //Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for it
    std::string aggregateEntry(std::vector<std::vector<std::string> > &dataToAggregate, unsigned int querySetKey,
                               unsigned int targetSetKey, unsigned int thread_idx) {
        

        std::string buffer;
        char keyBuffer[255];
        char *tmpBuff = Itoa::u32toa_sse2(targetSetKey, keyBuffer);
        buffer.append(keyBuffer, tmpBuff - keyBuffer - 1);
        buffer.append("\t");

        //the P-values of the (modified) truncated product method 
        //new theory: taking the best hit regardless of threshold and (from second hit on)sum of how much it surpassed threshold
        unsigned int orfCount = Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
        double logPvalThreshold = log(1.0/ (orfCount + 1));
        double minLogPval = 0;
        double sumLogPval = 0; 
        size_t k = 0;
        for (size_t i = 0; i < dataToAggregate.size(); ++i) {
            double logPvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
            double seqId = strtod(dataToAggregate[i][2].c_str(), NULL);
            int qStart = Util::fast_atoi<int>(dataToAggregate[i][4].c_str());
            int qEnd = Util::fast_atoi<int>(dataToAggregate[i][5].c_str());
            int qLen = Util::fast_atoi<int>(dataToAggregate[i][6].c_str());
            double qCov = 1.0* (qEnd - qStart + 1) / qLen;
            if(seqId == 1.0 && qCov == 1.0){
                logPvalThreshold = (logPvalue < logPvalThreshold) ? logPvalThreshold : logPvalue;
            };

        }
        for (size_t i = 0; i < dataToAggregate.size(); ++i) {
            double logPvalue = std::strtod(dataToAggregate[i][1].c_str(), NULL);
            if (logPvalue < minLogPval) {
                if (logPvalue == 0) {
                    //to avoid -0.0
                    minLogPval = logPvalue;
                }
                else {minLogPval = -logPvalue;}
            }
            if (logPvalue < logPvalThreshold) {
                //sum up the part exceeding logThreshold, add a minus to make score positive
                sumLogPval -= logPvalue - logPvalThreshold;
                k++;
            }
        }
        if(k == 0){
            //if no hit passed thr, take the -log of best hit pval as score
            buffer.append(SSTR(minLogPval));
            return buffer;
        }
        else {
            //if one or more hits passed thr
            buffer.append(SSTR(sumLogPval - logPvalThreshold));
            return buffer;
        }
    }

private:

    DBReader<unsigned int> *querySizeReader;
    DBReader<unsigned int> *targetSizeReader;
};

int combinescore(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);


    PvalueAggregator aggregation(par.db1, par.db2, par.db3, par.db4, (unsigned int) par.threads, par.compressed);
    return aggregation.run();
}
