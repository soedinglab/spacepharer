#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

const int CITATION_SPACEPHARER = CITATION_END;

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> filtermatchbyfdr; 
    std::vector<MMseqsParameter*> summarizeresults;    
    std::vector<MMseqsParameter*> predictmatchworkflow;
    std::vector<MMseqsParameter*> createsetdbworkflow;
    std::vector<MMseqsParameter*> downloadgenome;
    std::vector<MMseqsParameter*> combinescore;
    std::vector<MMseqsParameter*> combineprotnuclaln;
    std::vector<MMseqsParameter*> empiricalpval;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    PARAMETER(PARAM_REVERSE_SETDB)
    PARAMETER(PARAM_EXTRACTORF_SPACER)
    PARAMETER(PARAM_FDR_CUTOFF)
    PARAMETER(PARAM_FORMAT_TYPE)
    PARAMETER(PARAM_REPORT_PAM)
    PARAMETER(PARAM_PERFORM_NUCLALN)

    static const int FORMAT_TYPE_SHORT = 0;
    static const int FORMAT_TYPE_LONG = 1;
    static const int FORMAT_TYPE_ALN = 2;

    int reverseFragments;
    int extractorfsSpacer;
    int reverseSetDb;
    float fdrCutoff;
    int formatType;
    int reportPam;
    int performNuclAln;
    
private:
    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "Reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_REVERSE_SETDB(PARAM_REVERSE_SETDB_ID,"--reverse-setdb", "Create reversed setdb", "Create additional setDB with reversed fragments to compute under null [0,1]", typeid(int), (void *) &reverseSetDb, "^[0-1]{1}$"),
        PARAM_EXTRACTORF_SPACER(PARAM_EXTRACTORF_SPACER_ID,"--extractorf-spacer", "Extract orfs from spacers", "change parameter settings for extractorfs when createsetdb for spacer db [0,1]", typeid(int), (void *) &extractorfsSpacer, "^[0-1]{1}$"),
        PARAM_FDR_CUTOFF(PARAM_FDR_CUTOFF_ID,"--fdr", "FDR cutoff", "FDR cutoff for filtering matches[0.0, 1.0]", typeid(float), (void *) &fdrCutoff, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_FORMAT_TYPE(PARAM_FORMAT_TYPE_ID,"--fmt", "Output format", "0: short (only matches); 1: long (matches and hits); 2: long with nucleotide alignment", typeid(int), (void *) &formatType, "^[0-2]{1}$"),
        PARAM_REPORT_PAM(PARAM_REPORT_PAM_ID,"--report-pam", "Report PAM", "Report protospacer adjacent motifs up and downstream of hits [0,1]", typeid(int), (void *) &reportPam, "^[0-1]{1}$"),
        PARAM_PERFORM_NUCLALN(PARAM_PERFORM_NUCLALN_ID,"--perform-nucl-aln", "Perform nucl-nucl alignment", "Perform a nucl-nucl alignment in addition to the hits reported in 6-frame translated search, and assign the smallest E-value between the two alignments [0,1]", typeid(int), (void *) &performNuclAln, "^[0-1]{1}$")
    {
        downloadgenome.push_back(&PARAM_HELP);
        downloadgenome.push_back(&PARAM_HELP_LONG);
        downloadgenome.push_back(&PARAM_REVERSE_SETDB);
        downloadgenome.push_back(&PARAM_V);

        filtermatchbyfdr.push_back(&PARAM_FDR_CUTOFF);
        filtermatchbyfdr.push_back(&PARAM_COMPRESSED);
        filtermatchbyfdr.push_back(&PARAM_THREADS);
        filtermatchbyfdr.push_back(&PARAM_V);

        summarizeresults.push_back(&PARAM_FORMAT_TYPE);
        summarizeresults.push_back(&PARAM_COMPRESSED);
        summarizeresults.push_back(&PARAM_THREADS);
        summarizeresults.push_back(&PARAM_V);
        
        combinescore.push_back(&PARAM_COMPRESSED);
        combinescore.push_back(&PARAM_THREADS);
        combinescore.push_back(&PARAM_V);
                
        empiricalpval.push_back(&PARAM_COMPRESSED);
        empiricalpval.push_back(&PARAM_THREADS);
        empiricalpval.push_back(&PARAM_V);

        combineprotnuclaln.push_back(&PARAM_COMPRESSED);
        combineprotnuclaln.push_back(&PARAM_THREADS);
        combineprotnuclaln.push_back(&PARAM_V);
        
        createsetdbworkflow = combineList(createdb, translatenucs);
        createsetdbworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        createsetdbworkflow.push_back(&PARAM_EXTRACTORF_SPACER);
    
        predictmatchworkflow = combineList(searchworkflow, besthitbyset);
        predictmatchworkflow = combineList(predictmatchworkflow, combinepvalbyset);
        predictmatchworkflow = combineList(predictmatchworkflow, filtermatchbyfdr);
        predictmatchworkflow = combineList(predictmatchworkflow, summarizeresults);
        //predictmatchworkflow.push_back(&PARAM_FDR_CUTOFF);
        //predictmatchworkflow.push_back(&PARAM_FORMAT_TYPE);
        predictmatchworkflow.push_back(&PARAM_REPORT_PAM);
        predictmatchworkflow.push_back(&PARAM_PERFORM_NUCLALN);

        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;
        extractorfsSpacer = 0;
        reverseSetDb = 1;
        fdrCutoff = 0.05;
        formatType = 1;
        reportPam = 1;
        performNuclAln = 1;
        
        citations.emplace(CITATION_SPACEPHARER, "Zhang R, Mirdita M, Levy Karin E, Norroy C, Galiez C, and Soding J: SpacePHARER: Sensitive identification of phages from CRISPR spacers in prokaryotic hosts. (2020)");
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif
