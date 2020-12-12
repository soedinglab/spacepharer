#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>
#endif

int empiricalpval(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // std::string sizeDBName = par.db1 + "_set_size";
    // std::string sizeDBIndex = par.db1 + "_set_size.index";
    // DBReader<unsigned int> targetSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    // targetSizeReader.open(DBReader<unsigned int>::NOSORT);
    // const size_t numTargetSets = targetSizeReader.getSize();  
    // targetSizeReader.close();

    // sizeDBName = par.db3 + "_set_size";
    // sizeDBIndex = par.db3 + "_set_size.index";
    // DBReader<unsigned int> controlSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    // controlSizeReader.open(DBReader<unsigned int>::NOSORT);
    // const size_t numControlSets = controlSizeReader.getSize();
    // controlSizeReader.close();

    DBReader<unsigned int> posEvalDb(par.db1.c_str(),par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	posEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug::Progress progress(posEvalDb.getSize());

    //sort cEval of neg evalDB
    DBReader<unsigned int> negEvalDb(par.db2.c_str(),par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    negEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::vector<double> negToSort;

#pragma omp parallel
        {
            int thread_idx = 0;
#ifdef OPENMP
            thread_idx = omp_get_thread_num();
#endif

            std::vector<double> threadNegToSort;
            const char *entry[255];
#pragma omp for schedule(dynamic, 10)
            //TODO: check if control cEval list is empty
            for (size_t id = 0; id < negEvalDb.getSize(); id++) {
                progress.updateProgress();
                char *data = negEvalDb.getData(id, thread_idx);
                while (*data != '\0'){
                    size_t columns = Util::getWordsOfLine(data, entry, 255);
                    if (columns < 2) {
                        Debug(Debug::ERROR) << "Invalid result record.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    double eval = strtod(entry[1], NULL);
                    threadNegToSort.emplace_back(eval);
                    data = Util::skipLine(data);
                }
            }
#pragma omp critical
            negToSort.insert(negToSort.end(), threadNegToSort.begin(), threadNegToSort.end());
        }
    SORT_PARALLEL(negToSort.begin(), negToSort.end());

    negEvalDb.close();

    //compute empirical P-val for each pos entry based on neg evalDB
    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(),par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();
    progress.reset(posEvalDb.getSize());
    #pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(1000000);
        const char *entry[255];
        // const double evalCoef = 1.0 ;// * numTargetSets / numControlSets;
        const size_t nNull = negToSort.size();
        double empPval;
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < posEvalDb.getSize(); id++) {
            progress.updateProgress();
            char *data = posEvalDb.getData(id, thread_idx);
            while (*data != '\0') {
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid result record.\n";
                    EXIT(EXIT_FAILURE);
                }
                double eval = strtod(entry[1], NULL);
                //count how many negEval are lower than eval
                int left = 0;
                int right = nNull - 1;
                int neg_counter = nNull;
                while (left <= right) {
                    int mid = left + (right - left) / 2;
                    // Check if middle element is less than or equal to key 
                    if (negToSort[mid] >= eval) { 
                        // At least (mid + 1) elements are there whose values are less than or equal to key 
                        neg_counter = mid;
                        right = mid - 1;
                    } 
                    // If key is smaller, ignore right half 
                    else
                        left = mid + 1;
                    }
                    //empPval = 1.0 * ( neg_counter) / ( nNull);
                    empPval = ( (nNull - neg_counter) + 0.5 ) / ( nNull + 1 );
                const char *current = data;
                data = Util::skipLine(data);
                size_t length = data - current;
                std::string line(current, length - 1);
                buffer.append(line);
                buffer.append("\t");
                buffer.append(SSTR(empPval));
                buffer.append("\n");
            }    
            writer.writeData(buffer.c_str(), buffer.length(), posEvalDb.getDbKey(id), thread_idx);
            buffer.clear();
        }
    }
    writer.close();
    posEvalDb.close();

    return EXIT_SUCCESS;
}

