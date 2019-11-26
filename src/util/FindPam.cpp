#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Orf.h"
#include "PatternCompiler.h"

#include <algorithm>
#include <climits>
//#include <regex>

//degenrate nucleotide code
// A 	Adenine 	A
// C 	Cytosine 	C
// G 	Guanine 	G
// T 	Thymine (DNA) 	T
// W 	Weak 	A/T
// S 	Strong 	C/G
// M 	Amino 	A/C
// K 	Keto 	G/T
// R 	Purine 	A/G
// Y 	Pyrimidine 	C/T
// B 	Not A 	C/G/T
// D 	Not C 	A/G/T
// H 	Not G 	A/C/T
// V 	Not T 	A/C/G
// N 	Any 	A/C/G/T

//adapted from Leenay & Beisel, 2017. Table 1, represented in guide-centric orientation (PAM is from the strand that matches crRNA)
//classification    PAM location    Consensus sequence
//1-I-A     5'  YCN
//1-I-B     5'  CCW
//1-I-C     5'  YYC
//1-I-E     5'  AWG
//1-I-F     5'  CC
//1-III-B   3'  MMA (RNA PAM)
//2-II-A    3'  NGG
//2-II-A    3'  NNAGAA
//2-II-A    3'  NNGRRT
//2-II-C    3'  NNNNGWWT
//2-V-A     5'  TTN
//2-V-B     5'  TTN
//2-VI      3'  D (protospacer flanking sequence)

std::pair<std::string, std::string> searchpamlist (std::string threePrimeStrand, std::string fivePrimeStrand){
    char fivePrime[11];
    strcpy(fivePrime, fivePrimeStrand.c_str()); 
    char threePrime[11];
    strcpy(threePrime, threePrimeStrand.c_str()); 
    regmatch_t pmatch[1];
    std::pair<std::string, std::string> result;

    PatternCompiler YCN("[TC]C[ACGT]");
    PatternCompiler CCW ("CC[AT]");
    PatternCompiler YYC ("[TC][TC]C");
    PatternCompiler AWG ("A[AT]G");
    PatternCompiler CC ("CC");
    PatternCompiler TTN ("TT[ACGT]");
    if(YCN.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == 10) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(CCW.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == 10) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(YYC.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == 10) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(AWG.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == 10) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(CC.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == 10) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(TTN.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == 10) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else {
        result.first.append("-");
    }

    PatternCompiler MMA ("[AC][AC]A");
    PatternCompiler NGG ("[ACGT]GG");
    PatternCompiler NNAGAA ("[ACGT][ACGT]AGAA");
    PatternCompiler NNGRRT ("[ACGT][ACGT]G[AG][AG]T");
    PatternCompiler NNNNGWWT ("[ACGT][ACGT][ACGT][ACGT]G[AT][AT]T");
    PatternCompiler D ("[AGT]");
    // if(MMA.isMatch(threePrime, 1, pmatch)) {
    //     if(pmatch[0].rm_so < matchStartThreePrime){
    //         matchStartThreePrime = pmatch[0].rm_so;
    //         matchEndThreePrime = pmatch[0].rm_eo - 1;
    //         matchStringThreePrime = threePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so);
    //     }
    // } 
    if(NGG.isMatch(threePrime, 1, pmatch) && pmatch[0].rm_so == 0) {
        result.second.append(threePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    }  else if(NNAGAA.isMatch(threePrime, 1, pmatch) && pmatch[0].rm_so == 0) {
        result.second.append(threePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(NNGRRT.isMatch(threePrime, 1, pmatch) && pmatch[0].rm_so == 0) {
        result.second.append(threePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(NNNNGWWT.isMatch(threePrime, 1, pmatch) && pmatch[0].rm_so == 0) {
        result.second.append(threePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else{
        result.second.append("-");
    }
    return result;
}

int findpam(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string targetDBName = std::string(par.db1);
    std::string targetDBIndex = std::string(par.db1) + ".index";
    DBReader<unsigned int> targetReader(targetDBName.c_str(), targetDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    targetReader.open(DBReader<unsigned int>::NOSORT);

    std::string setName = std::string(par.db1) + "_set_to_contig";
    std::string setIndex = std::string(par.db1) + "_set_to_contig.index";
    DBReader<unsigned int> setReader(setName.c_str(), setIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    setReader.open(DBReader<unsigned int>::NOSORT);

    std::string alnDBName = std::string(par.db2);
    std::string alnDBIndex = std::string(par.db2) + ".index";
    DBReader<unsigned int> alnReader(alnDBName.c_str(), alnDBIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::NOSORT);

    std::string outDb = par.db3;
    DBWriter writer(outDb.c_str(), (outDb + ".index").c_str(), 1, par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    writer.open();

    std::string buffer = "";
	buffer.reserve(1000000);
    for (size_t id = 0; id < alnReader.getSize(); id++){
        char *data = alnReader.getData(id, 0);
        while (*data != '\0'){
            char *current = data;
            data = Util::skipLine(data);
            size_t length = data - current;
            std::string line(current, length - 1);
            if (line.empty() == true) {
                continue;
            }
            std::vector<std::string> columns = Util::split(line, "\t");
            size_t tsetid = Util::fast_atoi<size_t>(columns[0].c_str());
            size_t contigid = Util::fast_atoi<size_t>(setReader.getData(tsetid, 0));
            char *data = targetReader.getData(contigid, 0);
            size_t qstart = Util::fast_atoi<size_t>(columns[5].c_str()) - 1;
            size_t qend = Util::fast_atoi<size_t>(columns[6].c_str()) - 1;
            size_t qlen = Util::fast_atoi<size_t>(columns[7].c_str());
            size_t tstart = Util::fast_atoi<size_t>(columns[8].c_str()) - 1;
            size_t tend = Util::fast_atoi<size_t>(columns[9].c_str()) - 1;
            bool isqReverse = (qstart > qend) ? true : false;
            bool istReverse = (tstart > tend) ? true : false;
            std::string fivePrimeStrand = "";
            std::string threePrimeStrand = "";
            size_t fivePrimeEndPos;
            size_t threePrimeEndPos;
            if (isqReverse == false && istReverse == false) {
                fivePrimeEndPos = tstart - qstart;
                threePrimeEndPos = tend + (qlen - qend);
                for (size_t i = fivePrimeEndPos - 10; i < fivePrimeEndPos; ++i){
                    char seqChar = data[i];
                    fivePrimeStrand.append(1, seqChar);
                }
                for (size_t i = threePrimeEndPos; i < threePrimeEndPos + 10; ++i){
                    char seqChar = data[i];
                    threePrimeStrand.append(1, seqChar);
                }
            } else if (isqReverse == false && istReverse == true){   
                threePrimeEndPos = tend - (qlen - qend);
                fivePrimeEndPos = tstart + qstart;
                for (size_t i = fivePrimeEndPos + 10; i > fivePrimeEndPos; i--){
                    char seqChar = Orf::complement(data[i]);
                    fivePrimeStrand.append(1, seqChar);
                }
                for (size_t i = threePrimeEndPos; i > threePrimeEndPos -10; i--){
                    char seqChar = Orf::complement(data[i]);
                    threePrimeStrand.append(1, seqChar);
                }
            } else if (isqReverse == true && istReverse == false){
                threePrimeEndPos = tend + qend + 1;
                fivePrimeEndPos = tstart - (qlen - qstart) + 1;
                for (size_t i = fivePrimeEndPos -10; i < fivePrimeEndPos; ++i){
                    char seqChar = data[i];
                    fivePrimeStrand.append(1, seqChar);
                }
                for (size_t i = threePrimeEndPos; i < threePrimeEndPos + 10; ++i){
                    char seqChar = data[i];
                    threePrimeStrand.append(1, seqChar);
                }
            } else if (isqReverse == true && istReverse == true){
                threePrimeEndPos = tend - qend - 1;
                fivePrimeEndPos = tstart + (qlen - qstart) - 1;
                for (size_t i = fivePrimeEndPos + 10; i > fivePrimeEndPos; i--){
                    char seqChar = Orf::complement(data[i]);
                    fivePrimeStrand.append(1, seqChar);
                }
                for (size_t i = threePrimeEndPos; i > threePrimeEndPos -10; i--){
                    char seqChar = Orf::complement(data[i]);
                    threePrimeStrand.append(1, seqChar);
                }
            }
            std::string na = "-";
            std::pair<std::string, std::string> pamResult = searchpamlist(threePrimeStrand, fivePrimeStrand);
            buffer.append(line);
            buffer.append("\t");
            if(pamResult.first.c_str() == na){
                buffer.append("-");
            } else {
                buffer.append(SSTR(pamResult.first));
                // buffer.append("|");
                // buffer.append(SSTR(pamResult.first.second.first));
                // buffer.append("|");
                // buffer.append(SSTR(pamResult.first.second.second));
            }

            buffer.append("|");

            if(pamResult.second.c_str() == na){
                buffer.append("-");
            } else {
                buffer.append(SSTR(pamResult.second));
                // buffer.append("|");
                // buffer.append(SSTR(pamResult.second.second.first));
                // buffer.append("|");
                // buffer.append(SSTR(pamResult.second.second.second));
            }                    
            if (buffer.back() != '\n'){
                buffer.append("\n");
            }
        }    
            writer.writeData(buffer.c_str(), buffer.length(), alnReader.getDbKey(id));
            buffer.clear();
    }


    writer.close();
    alnReader.close();
    targetReader.close();

    return EXIT_SUCCESS;
}
