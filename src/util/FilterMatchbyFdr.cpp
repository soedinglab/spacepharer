#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#include <algorithm>
#include <climits>

int filtermatchbyfdr(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);
    std::string sizeDBName = std::string(par.db1) + "_set_size";
    std::string sizeDBIndex = std::string(par.db1) + "_set_size.index";
    DBReader<unsigned int> targetSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    targetSizeReader.open(DBReader<unsigned int>::NOSORT);

    sizeDBName = std::string(par.db3) + "_set_size";
    sizeDBIndex = std::string(par.db3) + "_set_size.index";
    DBReader<unsigned int> controlSizeReader(sizeDBName.c_str(), sizeDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    controlSizeReader.open(DBReader<unsigned int>::NOSORT);

    const size_t numTargetSets = targetSizeReader.getSize();  
    const size_t numControlSets = controlSizeReader.getSize();

    targetSizeReader.close();
    controlSizeReader.close();

    DBReader<unsigned int> posEvalDb(par.db2.c_str(),(std::string(par.db2).append(".index")).c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
	posEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> negEvalDb(par.db4.c_str(),(std::string(par.db4).append(".index")).c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    negEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    std::string outDb = par.db5;
    DBWriter writer(outDb.c_str(), (outDb + ".index").c_str(), 1, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();


//get the cEvals from pos and neg evalDB
    std::vector<double> posToSort;
    for (size_t id = 0; id < posEvalDb.getSize(); id++){
        char *data = posEvalDb.getData(id, 0);
        while (*data != '\0'){
            char *current = data;
            data = Util::skipLine(data);
            size_t length = data - current;
            std::string line(current, length - 1);
            if (line.empty() == true) {
                continue;
            }
            std::vector<std::string> columns = Util::split(line, "\t");
            double eval = strtod(columns[1].c_str(), NULL);
            posToSort.push_back(eval);
            }    
    }
    std::sort(posToSort.begin(), posToSort.end());

    std::vector<double> negToSort;
    for (size_t id = 0; id < negEvalDb.getSize(); id++){
        char *data = negEvalDb.getData(id, 0);
        while (*data != '\0'){
            char *current = data;
            data = Util::skipLine(data);
            size_t length = data - current;
            std::string line(current, length - 1);
            if (line.empty() == true) {
                continue;
            }
            std::vector<std::string> columns = Util::split(line, "\t");
            double eval = strtod(columns[1].c_str(), NULL);
            negToSort.push_back(eval);
            }        
    }
    std::sort(negToSort.begin(), negToSort.end());    


    size_t neg_counter = 0;
    size_t pos_counter = 0;
    double currentFDR = 0;
    double expectedFPcoef = (numTargetSets * posEvalDb.getSize() + numControlSets * negEvalDb.getSize()) / (numControlSets * negEvalDb.getSize());
    //TODO: factorize FDR, default 0.05
    while (currentFDR <= par.fdrCutoff) {
        pos_counter++;
        while (posToSort[pos_counter] >= negToSort[neg_counter]){
            neg_counter++;
        }
        currentFDR = neg_counter * expectedFPcoef / (pos_counter + neg_counter);
    }
    double threshold = posToSort[pos_counter];
    std::cout << "cEval threshold cutoff is " << threshold << " with FDR of " << par.fdrCutoff << ".\n";
    negToSort.clear();
    posToSort.clear();

    std::string buffer = "";
	buffer.reserve(1000000);
    for (size_t id = 0; id < posEvalDb.getSize(); id++){
        char *data = posEvalDb.getData(id, 0);
        while (*data != '\0'){
            char *current = data;
            data = Util::skipLine(data);
            size_t length = data - current;
            std::string line(current, length - 1);
            if (line.empty() == true) {
                continue;
            }
            int nomatch = 0;
            std::vector<std::string> columns = Util::split(line, "\t");
            double eval = strtod(columns[1].c_str(), NULL);
            nomatch = !(eval < threshold);
            if(!nomatch) {
                buffer.append(line);
                if (buffer.back() != '\n'){
                    buffer.append("\n");
                }
            }

        }    
            writer.writeData(buffer.c_str(), buffer.length(), posEvalDb.getDbKey(id));
            buffer.clear();
    }


    writer.close(true);
    posEvalDb.close();
    negEvalDb.close();


    return EXIT_SUCCESS;
}
