#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <climits>

#ifdef OPENMP
#include <omp.h>
#endif

int truncatebesthits(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string sizeDBName = std::string(par.db1) + "_set_size";
    std::string sizeDBIndex = std::string(par.db1) + "_set_size.index";
    DBReader<unsigned int> sizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    sizeReader.open(DBReader<unsigned int>::NOSORT);

    std::string setDBName = std::string(par.db1) + "_member_to_set";
    std::string setDBIndex = std::string(par.db1) + "_member_to_set.index";
    DBReader<unsigned int> setReader(setDBName.c_str(), setDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    setReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, resultReader.getDbtype());
    dbw.open();

    Debug::Progress progress(resultReader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            progress.updateProgress();
            char *data = resultReader.getData(i, thread_idx);
            unsigned int resultKey = resultReader.getDbKey(i);
            size_t setKey = Util::fast_atoi<size_t>(setReader.getDataByDBKey(resultKey, thread_idx));
            unsigned int setSize = Util::fast_atoi<unsigned int>(sizeReader.getDataByDBKey(setKey, thread_idx));
            while (*data != '\0') {
                char *current = data;
                data = Util::skipLine(data);
                size_t length = data - current;
                std::string line(current, length - 1);
                if (line.empty() == true) {
                    continue;
                }

                std::vector<std::string> columns = Util::split(line, "\t");
                //compute besthit P-val, column [3] now at positon [1]
                double logPval = strtod(columns[1].c_str(), NULL);
                double logPvalThr = log(1.0/(setSize + 1));
                if (logPval >= logPvalThr){
                    continue;
                }
                buffer.append(columns[0]);
                buffer.append("\t");
                buffer.append(columns[3]);
                buffer.append("\t");
                buffer.append(columns[2]);
                buffer.append("\t");
                buffer.append(SSTR(exp(logPval)));
                buffer.append("\t");
                buffer.append(columns[4]);
                buffer.append("\t");
                buffer.append(columns[5]);
                buffer.append("\t");
                buffer.append(columns[6]);
                buffer.append("\t");
                buffer.append(columns[7]);
                buffer.append("\t");
                buffer.append(columns[8]);
                buffer.append("\t");
                buffer.append(columns[9]);
                buffer.append("\t");
                buffer.append(columns[10]);
                buffer.append("\n");
            }
            dbw.writeData(buffer.c_str(), buffer.length(), resultKey, thread_idx);
            buffer.clear();
        }
    }
    dbw.close();
    resultReader.close();
    setReader.close();
    sizeReader.close();

    return EXIT_SUCCESS;
}

