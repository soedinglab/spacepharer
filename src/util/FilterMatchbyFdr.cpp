#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include "omptl/omptl_algorithm"

#ifdef OPENMP
#include <omp.h>
#endif

int filtermatchbyfdr(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string sizeDBName = par.db1 + "_set_size";
    std::string sizeDBIndex = par.db1 + "_set_size.index";
    DBReader<unsigned int> targetSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    targetSizeReader.open(DBReader<unsigned int>::NOSORT);
    const size_t numTargetSets = targetSizeReader.getSize();  
    targetSizeReader.close();

    sizeDBName = par.db3 + "_set_size";
    sizeDBIndex = par.db3 + "_set_size.index";
    DBReader<unsigned int> controlSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    controlSizeReader.open(DBReader<unsigned int>::NOSORT);
    const size_t numControlSets = controlSizeReader.getSize();
    controlSizeReader.close();

    DBReader<unsigned int> posEvalDb(par.db2.c_str(),par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	posEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug::Progress progress(posEvalDb.getSize());

    //get the cEvals from pos and neg evalDB
    std::vector<double> posToSort;
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::vector<double> threadPosToSort;
        const char *entry[255];
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < posEvalDb.getSize(); id++){
            progress.updateProgress();
            char *data = posEvalDb.getData(id, thread_idx);
            while (*data != '\0'){
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid result record.\n";
                    EXIT(EXIT_FAILURE);
                }
                double eval = strtod(entry[1], NULL);
                threadPosToSort.emplace_back(eval);
                data = Util::skipLine(data);
            }
        }

#pragma omp critical
        posToSort.insert(posToSort.end(), threadPosToSort.begin(), threadPosToSort.end());
    }
    omptl::sort(posToSort.begin(), posToSort.end());

    DBReader<unsigned int> negEvalDb(par.db4.c_str(),par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    negEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    double threshold;
    if (negEvalDb.getSize() > 0) {
        progress.reset(negEvalDb.getSize());
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
            //check if control cEval list is empty
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
        omptl::sort(negToSort.begin(), negToSort.end());    

        size_t neg_counter = 0;
        size_t pos_counter = 0;
        double currentFDR = 0;
        double expectedFPcoef = 1.0 * (numTargetSets * posEvalDb.getSize() + numControlSets * negEvalDb.getSize()) / (numControlSets * negEvalDb.getSize());
        while (currentFDR <= par.fdrCutoff) {
            pos_counter++;
            while (posToSort[pos_counter] >= negToSort[neg_counter]) {
                neg_counter++;
                //check if neg_counter has run to the end
                if(neg_counter == negEvalDb.getSize()){
                    Debug(Debug::WARNING) << "FDR of " << par.fdrCutoff << " cannot be reached, combined e-value threshold is " << numTargetSets << "\n";
                    threshold = numTargetSets;
                    break;
                }
            }
            currentFDR = neg_counter * expectedFPcoef / (pos_counter + neg_counter);
        }
        threshold = posToSort[pos_counter];
        //check if pos_counter has run to the end
        if (pos_counter == posEvalDb.getSize()) {
            Debug(Debug::WARNING) << "Combined e-value threshold with FDR of " << par.fdrCutoff << "cannot be determined.\n";
        } else {
            if (neg_counter < negEvalDb.getSize()){
                Debug(Debug::INFO) << "Combined e-value threshold is " << threshold << " with FDR of " << par.fdrCutoff << ".\n"; 
            }
            Debug(Debug::INFO) << pos_counter << " matches passed combined e-value threshold.\n";   
        }
    } else {
        Debug(Debug::WARNING) << "Combined e-value list of control set is empty, combined e-value threshold is " << numTargetSets << "\n";  
        threshold = numTargetSets;
    }
    negEvalDb.close();
    posToSort.clear();

    DBWriter writer(par.db5.c_str(), par.db5Index.c_str(),par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
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
                const char* start = data;
                data = Util::skipLine(data);
                if(eval < threshold) {
                    buffer.append(start, data - start);
                }
            }    
            writer.writeData(buffer.c_str(), buffer.length(), posEvalDb.getDbKey(id), thread_idx);
            buffer.clear();
        }
    }
    writer.close();
    posEvalDb.close();

    return EXIT_SUCCESS;
}
