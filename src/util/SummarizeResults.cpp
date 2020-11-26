#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int summarizeresults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> matchReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    matchReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> alnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::NOSORT);

    const bool shouldCompress = par.dbOut == true && par.compressed == true;
    const int dbType = par.dbOut == true ? Parameters::DBTYPE_GENERIC_DB : Parameters::DBTYPE_OMIT_FILE;
    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, shouldCompress, dbType);
    dbw.open();

    Debug::Progress progress(matchReader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string buffer;
        buffer.reserve(1024 * 1024);
        std::string tmpBuffer;
        tmpBuffer.reserve(1024 * 1024);
        std::string querySetName;
        std::string targetSetName;

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < matchReader.getSize(); ++i) {
            progress.updateProgress();
            char *data = matchReader.getData(i, thread_idx);
            while (*data != '\0'){
                char *current = data;
                data = Util::skipLine(data);
                size_t length = data - current;
                std::string line(current, length - 1);
                if (line.empty() == true) {
                    continue;
                }
                std::vector<std::string> columns = Util::split(line, "\t");
                size_t targetSetId = Util::fast_atoi<size_t>(columns[0].c_str());
                std::string cScore = columns[1].c_str();
                size_t key = matchReader.getDbKey(i);
                char *alnData = alnReader.getDataByDBKey(key, thread_idx);

                size_t lineCount = 0;
                while (*alnData != '\0') {
                    char *alnCurrent = alnData;
                    alnData = Util::skipLine(alnData);
                    size_t length = alnData - alnCurrent;
                    std::string line(alnCurrent, length - 1);
                    if (line.empty() == true) {
                        continue;
                    }
                    std::vector<std::string> alnColumns = Util::split(line, "\t");
                    querySetName = alnColumns[2];
                    if (targetSetId == Util::fast_atoi<size_t>(alnColumns[0].c_str())){
                        targetSetName = alnColumns[3];
                        if (par.formatType == LocalParameters::FORMAT_TYPE_LONG || par.formatType == LocalParameters::FORMAT_TYPE_ALN){
                            tmpBuffer.append(">");
                            tmpBuffer.append(alnColumns[1].c_str()); //spacername
                            tmpBuffer.append("\t");
                            tmpBuffer.append(alnColumns[3].c_str()); //genomename
                            tmpBuffer.append("\t");
                            tmpBuffer.append(alnColumns[4].c_str()); //best-hit pval
                            tmpBuffer.append("\t");
                            tmpBuffer.append(alnColumns[5].c_str()); //qstart
                            tmpBuffer.append("\t");
                            tmpBuffer.append(alnColumns[6].c_str()); //qend
                            tmpBuffer.append("\t");
                            tmpBuffer.append(alnColumns[8].c_str()); //tstart
                            tmpBuffer.append("\t");
                            tmpBuffer.append(alnColumns[9].c_str()); //tend
                            if(alnColumns.size() == 13){
                                tmpBuffer.append("\t");
                                tmpBuffer.append(alnColumns[12].c_str()); //PAM                       
                            }
                            if (par.formatType == LocalParameters::FORMAT_TYPE_ALN){
                                tmpBuffer.append("\n");
                                tmpBuffer.append(alnColumns[10].c_str()); //qaln
                                tmpBuffer.append("\n");
                                tmpBuffer.append(alnColumns[11].c_str()); //taln
                            }
                            tmpBuffer.append("\n");
                        }
                        lineCount++;
                    }
                }
                if(lineCount > 0){
                    buffer.append("#");
                    buffer.append(querySetName);
                    buffer.append("\t");
                    buffer.append(targetSetName);
                    buffer.append("\t");
                    buffer.append(cScore.c_str());
                    buffer.append("\t");
                    buffer.append(SSTR(lineCount));
                    buffer.append("\n");
                    if (par.formatType == LocalParameters::FORMAT_TYPE_LONG || par.formatType == LocalParameters::FORMAT_TYPE_ALN){
                        buffer.append(tmpBuffer);
                    }
                }

                // if (buffer.back() != '\n'){
                //     buffer.append("\n");
                // }
                tmpBuffer.clear();
                targetSetName.clear();
                querySetName.clear();
            }
            dbw.writeData(buffer.c_str(), buffer.length(), matchReader.getDbKey(i), thread_idx, par.dbOut);
            buffer.clear();
        }
    }
    dbw.close(true);
    if (par.dbOut == false) {
        FileUtil::remove(par.db3Index.c_str());
    }
    matchReader.close();
    alnReader.close();

    return EXIT_SUCCESS;
}

