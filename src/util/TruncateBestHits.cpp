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
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string sizeDBName = std::string(par.db1) + "_set_size";
    std::string sizeDBIndex = std::string(par.db1) + "_set_size.index";
    DBReader<unsigned int> sizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    sizeReader.open(DBReader<unsigned int>::NOSORT);

    std::string setDBName = std::string(par.db1) + "_member_to_set";
    std::string setDBIndex = std::string(par.db1) + "_member_to_set.index";
    DBReader<unsigned int> setReader(setDBName.c_str(), setDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    setReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> resultReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    resultReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, resultReader.getDbtype());
    dbw.open();

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(10 * 1024);
        char dbKey[255];
#pragma omp for schedule(static)
        for (size_t i = 0; i < resultReader.getSize(); ++i) {
            char *data = resultReader.getData(i, thread_idx);
            // char *current = data;
            // data = Util::skipLine(data);
            // size_t length = data - current;
            // std::string line(current, length - 1);
            // Util::parseKey(data, dbKey);
            // unsigned int key = Util::fast_atoi<unsigned int>(dbKey);
            // size_t id = resultReader.getId(key);
            // go through the results in the cluster and add them to one entry
            while (*data != '\0'){
                char *current = data;
                data = Util::skipLine(data);
                size_t length = data - current;
                std::string line(current, length - 1);
                if (line.empty() == true) {
                    continue;
                }
                std::vector<std::string> columns = Util::split(line, "\t");
                double logPval = strtod(columns[1].c_str(), NULL);
                size_t id = resultReader.getId(i);
                size_t key = resultReader.getDbKey(id);
                size_t setKey = Util::fast_atoi<size_t>(setReader.getDataByDBKey(key, thread_idx));
                unsigned int setSize = Util::fast_atoi<unsigned int>(sizeReader.getDataByDBKey(setKey, thread_idx));
                double logPvalThr = log(1.0/(setSize + 1));
                if(logPval < logPvalThr) {
                    buffer.append(line);
                    if (buffer.back() != '\n'){
                        buffer.append("\n");
                    }
                }
            }
        dbw.writeData(buffer.c_str(), buffer.length(), resultReader.getDbKey(i), thread_idx);
        buffer.clear();
        }
    }
    dbw.close();
    resultReader.close();
    setReader.close();

    return EXIT_SUCCESS;
}

