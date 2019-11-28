#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <climits>


int summarizeresults(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    DBReader<unsigned int> matchReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    matchReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> alnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, alnReader.getDbtype());
    dbw.open();

    std::string buffer;
    buffer.reserve(1024 * 1024);
    for (size_t i = 0; i < matchReader.getSize(); ++i) {
        char *data = matchReader.getData(i, 0);
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
            std::string cEval = columns[1].c_str();
            size_t key = matchReader.getDbKey(i);
            char *alnData = alnReader.getDataByDBKey(key ,0);
            std::string tmpBuffer;
            std::string querySetName;
            std::string targetSetName;
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
                querySetName= alnColumns[2].c_str();
                if (targetSetId == Util::fast_atoi<size_t>(alnColumns[0].c_str())){
                    targetSetName = alnColumns[3].c_str();
                    if (par.formatType == LocalParameters::FORMAT_TYPE_LONG || par.formatType == LocalParameters::FORMAT_TYPE_ALN){
                        tmpBuffer.append(">");
                        tmpBuffer.append(alnColumns[1].c_str()); //spacername
                        tmpBuffer.append("\t");
                        tmpBuffer.append(alnColumns[3].c_str()); //genomename
                        tmpBuffer.append("\t");
                        tmpBuffer.append(alnColumns[4].c_str()); //eval
                        tmpBuffer.append("\t");
                        tmpBuffer.append(alnColumns[5].c_str()); //qstart
                        tmpBuffer.append("\t");
                        tmpBuffer.append(alnColumns[6].c_str()); //qend
                        //tmpBuffer.append("\t");
                        //tmpBuffer.append(alnColumns[10].c_str()); //qaln
                        //tmpBuffer.append(">");
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
                //buffer.append(SSTR(i));
                //buffer.append("\t");
                buffer.append(querySetName.c_str());
                //buffer.append("\t");
                //buffer.append(SSTR(targetSetId));
                buffer.append("\t");
                buffer.append(targetSetName.c_str());
                buffer.append("\t");
                buffer.append(cEval.c_str());
                buffer.append("\t");
                buffer.append(SSTR(lineCount));
                buffer.append("\n");
                if (par.formatType == LocalParameters::FORMAT_TYPE_LONG || par.formatType == LocalParameters::FORMAT_TYPE_ALN){
                    buffer.append(tmpBuffer);
                }
            }
            //size_t setKey = Util::fast_atoi<size_t>(alnReader.getDataByDBKey(key, 0));
            //unsigned int setSize = Util::fast_atoi<unsigned int>(alnReader.getDataByDBKey(setKey, 0));
            //double logPvalThr = log(1.0/(setSize + 1));
            // if(logPval < logPvalThr) {
            //     buffer.append(line);
            if (buffer.back() != '\n'){
                buffer.append("\n");
            }
            // }
            tmpBuffer.clear();
        }
        dbw.writeData(buffer.c_str(), buffer.length(), matchReader.getDbKey(i), 0);
        buffer.clear();
        }
    dbw.close(true);
    matchReader.close();
    alnReader.close();

    return EXIT_SUCCESS;
}

