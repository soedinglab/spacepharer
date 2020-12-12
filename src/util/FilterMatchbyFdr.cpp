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

int filtermatchbyfdr(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> posScoreDb(par.db1.c_str(),par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	posScoreDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug::Progress progress(posScoreDb.getSize());

    //get the cScores from pos and neg scoreDB
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
        for (size_t id = 0; id < posScoreDb.getSize(); id++){
            progress.updateProgress();
            char *data = posScoreDb.getData(id, thread_idx);
            while (*data != '\0'){
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid result record.\n";
                    EXIT(EXIT_FAILURE);
                }
                double score = strtod(entry[1], NULL);
                threadPosToSort.emplace_back(score);
                data = Util::skipLine(data);
            }
        }

#pragma omp critical
        posToSort.insert(posToSort.end(), threadPosToSort.begin(), threadPosToSort.end());
    }
    SORT_PARALLEL(posToSort.rbegin(), posToSort.rend());

    DBReader<unsigned int> negScoreDb(par.db2.c_str(),par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    negScoreDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    double threshold;
    if (negScoreDb.getSize() > 0) {
        progress.reset(negScoreDb.getSize());
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
            //check if control cScore list is empty
            for (size_t id = 0; id < negScoreDb.getSize(); id++) {
                progress.updateProgress();
                char *data = negScoreDb.getData(id, thread_idx);
                while (*data != '\0'){
                    size_t columns = Util::getWordsOfLine(data, entry, 255);
                    if (columns < 2) {
                        Debug(Debug::ERROR) << "Invalid result record.\n";
                        EXIT(EXIT_FAILURE);
                    }
                    double score = strtod(entry[1], NULL);
                    threadNegToSort.emplace_back(score);
                    data = Util::skipLine(data);
                }
            }
#pragma omp critical
            negToSort.insert(negToSort.end(), threadNegToSort.begin(), threadNegToSort.end());
        }
        SORT_PARALLEL(negToSort.rbegin(), negToSort.rend());

        //cumsum for TP+FP(pos_counter) and FP(neg_counter) and x = FP/nNeg, y = (TP + FP)/nPos
        std::vector<double> uniqueScoreList;
        std::vector<double> x;
        std::vector<double> y;
        size_t cnt = 0;
        double currentScore = DBL_MAX;
        size_t neg_counter = 0;
        size_t pos_counter = 0;
        //set the first point [0,0], set uniqueScoreList first element to 0.0
        //go through posToSort and find unique values, while counting # of predictions in pos and neg lists above this value
        for (size_t i = 0; i < posToSort.size(); i++) {
            if (posToSort[pos_counter] < currentScore) {
                currentScore = posToSort[pos_counter];
                while (currentScore < (negToSort[neg_counter]) && neg_counter < negToSort.size()) {
                    neg_counter++;
                }
                uniqueScoreList.push_back(currentScore);
                y.push_back(1.0 * pos_counter / posToSort.size());
                x.push_back(1.0 * (neg_counter + 0.5)/ (negToSort.size() + 1));
                cnt++ ;
            }
            pos_counter++;
        }
        x.push_back(1.0);
        y.push_back(1.0);

        size_t i = 0;
        double slopeMax = 0;
        std::vector<double> slopeList;
        std::vector<size_t> idxList;
        double slope;
        while (i < (x.size() - 1)) {
            size_t jMax = i + 1;
            for(size_t j = i + 1; j <  x.size() ; j++){
                slope = 1.0 * (y[j] - y[i]) / (x[j] - x[i]);
                if (slope >= slopeMax){
                    jMax = j;
                    slopeMax = slope;
                } 
            } 
            i = jMax;
            slopeList.push_back(slopeMax);
            idxList.push_back(jMax);
            slopeMax = 0;
        }

        i = 0;
        //select the average of last and second last slope as pi0
        double pi0 = (slopeList.end()[-2] + slopeList.end()[-1])/2;
        double deltaX = 0;
        double currentFDR = 0;
        while (currentFDR <= par.fdrCutoff) {
            currentFDR = x[idxList[i]] * pi0/ y[idxList[i]];
            i++;
        }

        if (i < 2) {
            threshold = posToSort[posToSort.size()];
            Debug(Debug::WARNING) << "Combined score list too short. Using threshold " << threshold << "\n";
        } else {
            size_t j = idxList[i-2];
            double TPFP = y[j];
            double FP = x[j]* pi0;
            currentFDR = 0;
            while (currentFDR <= par.fdrCutoff){
                j++;
                deltaX = x[j] - x[j-1];
                TPFP += deltaX * slopeList[i-1];
                FP += deltaX * pi0;
                currentFDR = FP / TPFP;
            }

            threshold = uniqueScoreList[j];
            Debug(Debug::INFO) << "Combined score threshold is " << threshold << " with FDR of " << par.fdrCutoff << ".\n";
            Debug(Debug::INFO) <<  y[j]* posToSort.size() << " matches passed combined score threshold.\n";
        }
    } else {
        Debug(Debug::WARNING) << "Combined score list of control set is empty\n";
        threshold = posToSort[posToSort.size()];
    }
    negScoreDb.close();
    posToSort.clear();

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();
    progress.reset(posScoreDb.getSize());
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
        for (size_t id = 0; id < posScoreDb.getSize(); id++) {
            progress.updateProgress();
            char *data = posScoreDb.getData(id, thread_idx);
            while (*data != '\0') {
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid result record.\n";
                    EXIT(EXIT_FAILURE);
                }
                double score = strtod(entry[1], NULL);
                const char* start = data;
                data = Util::skipLine(data);
                if(score >= threshold) {
                    buffer.append(start, data - start);
                }
            }    
            writer.writeData(buffer.c_str(), buffer.length(), posScoreDb.getDbKey(id), thread_idx);
            buffer.clear();
        }
    }
    writer.close();
    posScoreDb.close();

    return EXIT_SUCCESS;
}

