#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int combineprotnuclaln(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> protAlnReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    protAlnReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBReader<unsigned int> nuclAlnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    nuclAlnReader.open(DBReader<unsigned int>::NOSORT);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, protAlnReader.getDbtype());
    dbw.open();
    
    Debug::Progress progress(protAlnReader.getSize());

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string buffer;
        buffer.reserve(1024 * 1024);

#pragma omp for schedule(dynamic, 10)

        for (size_t i = 0; i < protAlnReader.getSize(); ++i) {
            
            progress.updateProgress();
            char *protData = protAlnReader.getData(i, thread_idx);
            while (*protData != '\0'){
                char *protCurrent = protData;
                protData = Util::skipLine(protData);
                size_t length = protData - protCurrent;
                std::string line(protCurrent, length - 1);
                if (line.empty() == true) {
                    continue;
                }
                std::vector<std::string> columns = Util::split(line, "\t");
                size_t targetId = Util::fast_atoi<size_t>(columns[0].c_str());
                double protEval = strtod(columns[3].c_str(), NULL);
                size_t key = protAlnReader.getDbKey(i);
                char *nuclData = nuclAlnReader.getDataByDBKey(key, thread_idx);

                while (*nuclData != '\0') {
                    char *nuclCurrent = nuclData;
                    nuclData = Util::skipLine(nuclData);
                    size_t length = nuclData - nuclCurrent;
                    std::string nuclLine(nuclCurrent, length - 1);
                    if (nuclLine.empty() == true) {
                        continue;
                    }
                    std::vector<std::string> nuclColumns = Util::split(nuclLine, "\t");
                    double nuclEval = strtod(nuclColumns[3].c_str(), NULL);
                    if (targetId == Util::fast_atoi<size_t>(nuclColumns[0].c_str())){   
                        columns[3] = SSTR((protEval < nuclEval) ? protEval : nuclEval);
                    }
                }
                buffer.append(columns[0]);
                buffer.append("\t");
                buffer.append(columns[1]);
                buffer.append("\t");
                buffer.append(columns[2]);
                buffer.append("\t");
                buffer.append(columns[3]);
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
            dbw.writeData(buffer.c_str(), buffer.length(), protAlnReader.getDbKey(i), thread_idx);
            buffer.clear();

        }

    }
    dbw.close();
    protAlnReader.close();
    nuclAlnReader.close();

    return EXIT_SUCCESS;
}

// #include "LocalParameters.h"
// #include "FileUtil.h"
// #include "Matcher.h"
// #include "NucleotideMatrix.h"
// #include "SubstitutionMatrix.h"
// #include "DBReader.h"
// #include "DBWriter.h"
// #include "Debug.h"
// #include "Util.h"

// #ifdef OPENMP
// #include <omp.h>
// #endif

// int combineprotnuclaln(int argc, const char **argv, const Command& command) {
//     LocalParameters& par = LocalParameters::getLocalInstance();
//     par.parseParameters(argc, argv, command, true, 0, 0);

//     DBReader<unsigned int> queryReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
//     queryReader.open(DBReader<unsigned int>::NOSORT);

//     DBReader<unsigned int> targetReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
//     targetReader.open(DBReader<unsigned int>::NOSORT);

//     std::string nuclDBName = std::string(par.db2) + "_nucl_orf";
//     std::string nuclDBIndex = std::string(par.db2) + "_nucl_orf.index";
//     DBReader<unsigned int> targetNuclReader(nuclDBName.c_str(), nuclDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
//     targetNuclReader.open(DBReader<unsigned int>::NOSORT);

//     DBReader<unsigned int> protAlnReader(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
//     protAlnReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

//     DBReader<unsigned int> nuclAlnReader(par.db4.c_str(), par.db4Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
//     nuclAlnReader.open(DBReader<unsigned int>::NOSORT);

//     DBWriter dbw(par.db5.c_str(), par.db5Index.c_str(), par.threads, par.compressed, protAlnReader.getDbtype());
//     dbw.open();

//     Debug::Progress progress(protAlnReader.getSize());
//     BaseMatrix *protSubMat;
//     BaseMatrix *nuclSubMat;
//     protSubMat = new SubstitutionMatrix(par.scoringMatrixFile.aminoacids, 2.0, 0.0);
//     EvalueComputation protEvaluer(targetReader.getAminoAcidDBSize(), protSubMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
//     nuclSubMat = new NucleotideMatrix(par.scoringMatrixFile.nucleotides, 1.0, 0.0);
//     EvalueComputation nuclEvaluer(targetNuclReader.getAminoAcidDBSize(), nuclSubMat, par.gapOpen.nucleotides, par.gapExtend.nucleotides);
// #pragma omp parallel
//     {
//         int thread_idx = 0;
// #ifdef OPENMP
//         thread_idx = omp_get_thread_num();
// #endif

//         std::string buffer;
//         buffer.reserve(1024 * 1024);

// #pragma omp for schedule(dynamic, 10)

//         for (size_t i = 0; i < protAlnReader.getSize(); ++i) {
//             progress.updateProgress();
//             char *protData = protAlnReader.getData(i, thread_idx);
//             while (*protData != '\0'){
//                 char *protCurrent = protData;
//                 protData = Util::skipLine(protData);
//                 size_t length = protData - protCurrent;
//                 std::string line(protCurrent, length - 1);
//                 if (line.empty() == true) {
//                     continue;
//                 }
//                 std::vector<std::string> columns = Util::split(line, "\t");
//                 size_t targetId = Util::fast_atoi<size_t>(columns[0].c_str());
//                 int protBitScore = Util::fast_atoi<int>(columns[1].c_str());
//                 int protRawScore = static_cast<int>(protEvaluer.computeRawScoreFromBitScore(protBitScore) + 0.5);
//                 size_t key = protAlnReader.getDbKey(i);
//                 char *nuclData = nuclAlnReader.getDataByDBKey(key, thread_idx);

//                 while (*nuclData != '\0') {
//                     char *nuclCurrent = nuclData;
//                     nuclData = Util::skipLine(nuclData);
//                     size_t length = nuclData - nuclCurrent;
//                     std::string nuclLine(nuclCurrent, length - 1);
//                     if (nuclLine.empty() == true) {
//                         continue;
//                     }
//                     std::vector<std::string> nuclColumns = Util::split(nuclLine, "\t");
//                     int nuclBitScore = Util::fast_atoi<int>(nuclColumns[1].c_str());
//                     int nuclRawScore = static_cast<int>(nuclEvaluer.computeRawScoreFromBitScore(nuclBitScore) + 0.5);
//                     if (targetId == Util::fast_atoi<size_t>(nuclColumns[0].c_str())){   
//                         //S_aa = 13.5503 + 0.9588 * S_nt
//                         int recomputedRawScore = (1.0 * protRawScore + 1.0 * nuclRawScore) / 2;
//                         int recomputedBitScore = static_cast<int>(protEvaluer.computeBitScore(recomputedRawScore) + 0.5);
//                         size_t querySeqId = queryReader.getId(key);
//                         int queryLen = queryReader.getSeqLen(querySeqId);
//                         double recomputedEvalue = protEvaluer.computeEvalue(recomputedRawScore, queryLen);;
//                         columns[1] = SSTR(recomputedBitScore);
//                         columns[3] = SSTR(recomputedEvalue);
//                         // buffer.append(columns[0]);
//                         // buffer.append("\t");
//                         // buffer.append(SSTR(protRawScore));
//                         // buffer.append("\t");
//                         // buffer.append(SSTR(nuclRawScore));
//                         // buffer.append("\n");
//                     }
//                 }
//                 buffer.append(columns[0]);
//                 buffer.append("\t");
//                 buffer.append(columns[1]);
//                 buffer.append("\t");
//                 buffer.append(columns[2]);
//                 buffer.append("\t");
//                 buffer.append(columns[3]);
//                 buffer.append("\t");
//                 buffer.append(columns[4]);
//                 buffer.append("\t");
//                 buffer.append(columns[5]);
//                 buffer.append("\t");
//                 buffer.append(columns[6]);
//                 buffer.append("\t");
//                 buffer.append(columns[7]);
//                 buffer.append("\t");
//                 buffer.append(columns[8]);
//                 buffer.append("\t");
//                 buffer.append(columns[9]);
//                 buffer.append("\t");
//                 buffer.append(columns[10]);
//                 buffer.append("\n");

//                 // if (buffer.back() != '\n'){
//                 //     buffer.append("\n");
//                 // }
//             }
//             dbw.writeData(buffer.c_str(), buffer.length(), protAlnReader.getDbKey(i), thread_idx);
//             buffer.clear();

//         }
//     }
//     delete protSubMat;
//     delete nuclSubMat;
//     dbw.close();
    
//     queryReader.close();
//     targetReader.close();
//     targetNuclReader.close();
//     protAlnReader.close();
//     nuclAlnReader.close();

//     return EXIT_SUCCESS;
// }
