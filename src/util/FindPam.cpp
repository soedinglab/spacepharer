#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "Orf.h"
#include "PatternCompiler.h"

#ifdef OPENMP
#include <omp.h>
#endif

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
//classification    PAM location    Consensus sequence  Rep. Genome
//1-I-A     5'  YCN     Sulfolobus solfataricus P2, Pyrococcus. furiosus
//1-I-B     5'  CCW     Clostridium difficile
//1-I-C     5'  YYC     Bacillus halodurans
//1-I-E     5'  AWG     Escherichia coli K12, salmonella enterica
//1-I-F     5'  CC      Pseudomonas aeruginosa
//1-III-B   3'  MMA (RNA PAM)
//2-II-A    3'  NGG     Streptococcus pyogenes
//2-II-A    3'  NNAGAA  Streptococcus thermophilus(CRISPR1)
//2-II-A    3'  NNGRRT  Staphylococcus aureus
//2-II-C    3'  NNNNGWWT    Neisseria meningitidis
//2-V-A     5'  TTN     Francisella novicida
//2-V-B     5'  TTN     Alicyclobacillus acidoterrestris
//2-VI      3'  D (protospacer flanking sequence)

std::pair<std::string, std::string> searchpamlist(const std::string &threePrimeStrand, const std::string &fivePrimeStrand, int flankingSeqLen){
    char fivePrime[11];
    strcpy(fivePrime, fivePrimeStrand.c_str()); 
    char threePrime[11];
    strcpy(threePrime, threePrimeStrand.c_str()); 
    regmatch_t pmatch[1];
    std::pair<std::string, std::string> result;

    static PatternCompiler YCN("[TC]C[ACGT]");
    static PatternCompiler CCW ("CC[GAT]");
    static PatternCompiler YYC ("[TC][TC]C");
    static PatternCompiler AWG ("A[AT]G");
    static PatternCompiler CC ("CC");
    static PatternCompiler TTN ("TT[ACGT]");
    if(YCN.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == flankingSeqLen) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(CCW.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == flankingSeqLen) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(YYC.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == flankingSeqLen) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(CC.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == flankingSeqLen) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(AWG.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == flankingSeqLen) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else if(TTN.isMatch(fivePrime, 1, pmatch) && pmatch[0].rm_eo == flankingSeqLen) {
        result.first.append(fivePrimeStrand.substr(pmatch[0].rm_so, pmatch[0].rm_eo - pmatch[0].rm_so));
    } else {
        result.first.append("-");
    }

    // static PatternCompiler MMA ("[AC][AC]A");
    static PatternCompiler NGG ("[ACGT]GG");
    static PatternCompiler NNAGAA ("[ACGT][ACGT]AGAA");
    static PatternCompiler NNGRRT ("[ACGT][ACGT]G[AG][AG]T");
    static PatternCompiler NNNNGWWT ("[ACGT][ACGT][ACGT][ACGT]G[AT][AT]T");
    // static PatternCompiler D ("[AGT]");

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

std::string reverseComplement(std::string seq){
    unsigned int lenSeq = seq.length();
    std::string revStr = "";
    for (size_t i = 0; i < lenSeq; ++i) {
        size_t revInd = lenSeq - i - 1;
        char seqChar = Orf::complement(seq[revInd]);
        revStr.append(1, seqChar);
    }
    return revStr;
}

int findpam(int argc, const char **argv, const Command& command) {
    LocalParameters& par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> targetReader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    targetReader.open(DBReader<unsigned int>::NOSORT);

    std::string setName = std::string(par.db1) + "_set_to_contig";
    std::string setIndex = std::string(par.db1) + "_set_to_contig.index";
    DBReader<unsigned int> setReader(setName.c_str(), setIndex.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    setReader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> alnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alnReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    writer.open();

    Debug::Progress progress(alnReader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::string buffer;
        buffer.reserve(1000000);
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < alnReader.getSize(); id++) {
            progress.updateProgress();
            char *data = alnReader.getData(id, thread_idx);
            while (*data != '\0') {
                char *current = data;
                data = Util::skipLine(data);
                size_t length = data - current;
                std::string line(current, length - 1);
                if (line.empty() == true) {
                    continue;
                }
                std::vector<std::string> columns = Util::split(line, "\t");
                size_t tsetid = Util::fast_atoi<size_t>(columns[0].c_str());
                size_t contigid = Util::fast_atoi<size_t>(setReader.getDataByDBKey(tsetid, thread_idx));
                char *data = targetReader.getDataByDBKey(contigid, thread_idx);
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
                    for (size_t i = fivePrimeEndPos - par.flankingSeqLen; i < fivePrimeEndPos; ++i){
                        char seqChar = data[i];
                        fivePrimeStrand.append(1, seqChar);
                    }
                    for (size_t i = threePrimeEndPos; i < threePrimeEndPos + par.flankingSeqLen; ++i){
                        char seqChar = data[i];
                        threePrimeStrand.append(1, seqChar);
                    }
                } else if (isqReverse == false && istReverse == true){   
                    threePrimeEndPos = tend - (qlen - qend);
                    fivePrimeEndPos = tstart + qstart;
                    for (size_t i = fivePrimeEndPos + par.flankingSeqLen; i > fivePrimeEndPos; i--){
                        char seqChar = Orf::complement(data[i]);
                        fivePrimeStrand.append(1, seqChar);
                    }
                    for (size_t i = threePrimeEndPos; i > threePrimeEndPos - par.flankingSeqLen; i--){
                        char seqChar = Orf::complement(data[i]);
                        threePrimeStrand.append(1, seqChar);
                    }
                } else if (isqReverse == true && istReverse == false){
                    fivePrimeEndPos = tend + qend;
                    threePrimeEndPos = tstart - (qlen - qstart);
                    for (size_t i = fivePrimeEndPos + par.flankingSeqLen; i > fivePrimeEndPos; i--){
                        char seqChar = Orf::complement(data[i]);
                        fivePrimeStrand.append(1, seqChar);
                    }
                    for (size_t i = threePrimeEndPos; i > threePrimeEndPos - par.flankingSeqLen; i--){
                        char seqChar = Orf::complement(data[i]);
                        threePrimeStrand.append(1, seqChar);
                    }
                } else if (isqReverse == true && istReverse == true){
                    fivePrimeEndPos = tend - qend;
                    threePrimeEndPos = tstart + (qlen - qstart);
                    for (size_t i = fivePrimeEndPos - par.flankingSeqLen; i < fivePrimeEndPos; ++i){
                        char seqChar = data[i];
                        fivePrimeStrand.append(1, seqChar);
                    }
                    for (size_t i = threePrimeEndPos; i < threePrimeEndPos + par.flankingSeqLen; ++i){
                        char seqChar = data[i];
                        threePrimeStrand.append(1, seqChar);
                    }
                }
                std::string na = "-";
                std::pair<std::string, std::string> pamResultFwdStrand = searchpamlist(threePrimeStrand, fivePrimeStrand, par.flankingSeqLen);
                std::string threePrimeReverseStrand = reverseComplement(fivePrimeStrand);
                std::string fivePrimeReverseStrand = reverseComplement(threePrimeStrand);
                std::pair<std::string, std::string> pamResultRevStrand = searchpamlist(threePrimeReverseStrand, fivePrimeReverseStrand, par.flankingSeqLen);
                buffer.append(line);
                buffer.append("\t");
                if(pamResultFwdStrand.first.c_str() == na){
                    buffer.append("-");
                } else {
                    buffer.append(SSTR(pamResultFwdStrand.first));
                }

                buffer.append("|");

                if(pamResultFwdStrand.second.c_str() == na){
                    buffer.append("-");
                } else {
                    buffer.append(SSTR(pamResultFwdStrand.second));
                }                        
                                
                buffer.append("\t");
                if(pamResultRevStrand.first.c_str() == na){
                    buffer.append("-");
                } else {
                    buffer.append(SSTR(pamResultRevStrand.first));
                }

                buffer.append("|");

                if(pamResultRevStrand.second.c_str() == na){
                    buffer.append("-");
                } else {
                    buffer.append(SSTR(pamResultRevStrand.second));
                }                    
                if (buffer.back() != '\n'){
                    buffer.append("\n");
                }


            }    
            writer.writeData(buffer.c_str(), buffer.length(), alnReader.getDbKey(id), thread_idx);
            buffer.clear();
        }
    }
    writer.close();
    alnReader.close();
    targetReader.close();

    return EXIT_SUCCESS;
}
