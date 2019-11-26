#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

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

    std::vector<MMseqsParameter*> parsespacer;
    std::vector<MMseqsParameter*> filtermatchbyfdr;
    std::vector<MMseqsParameter*> findpam;
    std::vector<MMseqsParameter*> truncatebesthits;    
    std::vector<MMseqsParameter*> summarizeresults;    
    std::vector<MMseqsParameter*> predictmatchworkflow;
    std::vector<MMseqsParameter*> createsetdbworkflow;
    std::vector<MMseqsParameter*> downloadgenomeworkflow;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    PARAMETER(PARAM_REVERSE_SETDB)
    PARAMETER(PARAM_EXTRACTORF_SPACER)
    PARAMETER(PARAM_FDR_CUTOFF)
    PARAMETER(PARAM_FORMAT_TYPE)
    PARAMETER(PARAM_REPORT_PAM)

    static const int FORMAT_TYPE_SHORT = 0;
    static const int FORMAT_TYPE_LONG = 1;
    static const int FORMAT_TYPE_ALN = 2;

    int reverseFragments;
    int extractorfsSpacer;
    int reverseSetDb;
    float fdrCutoff;
    int formatType;
    int reportPam;
    
private:
    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_REVERSE_SETDB(PARAM_REVERSE_SETDB_ID,"--reverse-setDB", "Create setDB with reversed fragments", "Create additional setDB with reversed fragments to compute under null [0,1]", typeid(int), (void *) &reverseSetDb, "^[0-1]{1}$"),
        PARAM_EXTRACTORF_SPACER(PARAM_EXTRACTORF_SPACER_ID,"--extractorf-spacer", "Extract orfs from spacers", "change parameter settings for extractorfs when createsetdb for spacer db [0,1]", typeid(int), (void *) &extractorfsSpacer, "^[0-1]{1}$"),
        PARAM_FDR_CUTOFF(PARAM_FDR_CUTOFF_ID,"--fdr", "FDR cutoff", "FDR cutoff for filtering matches[0.0, 1.0]", typeid(float), (void *) &fdrCutoff, "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PARAM_FORMAT_TYPE(PARAM_FORMAT_TYPE_ID,"--fmt", "output format", "0: short (only matches); 1: long (matches and hits); 2: long with nucleotide alignment", typeid(int), (void *) &formatType, "^[0-2]{1}$"),
        PARAM_REPORT_PAM(PARAM_REPORT_PAM_ID,"--report-pam", "report PAM", "report protospacer adjacent motifs up and downstream of hits [0,1]", typeid(int), (void *) &reportPam, "^[0-1]{1}$")
       
    {


        createsetdbworkflow = combineList(createdb, extractorfs);
        createsetdbworkflow = combineList(createsetdbworkflow, translatenucs);
        createsetdbworkflow = combineList(createsetdbworkflow, result2stats);
        createsetdbworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        createsetdbworkflow.push_back(&PARAM_EXTRACTORF_SPACER);
    
        predictmatchworkflow = combineList(predictmatchworkflow, searchworkflow);
        predictmatchworkflow = combineList(predictmatchworkflow, besthitbyset);
        predictmatchworkflow = combineList(predictmatchworkflow, combinepvalbyset);
        predictmatchworkflow = combineList(predictmatchworkflow, filtermatchbyfdr);
        predictmatchworkflow = combineList(predictmatchworkflow, summarizeresults);
        predictmatchworkflow.push_back(&PARAM_FDR_CUTOFF);
        predictmatchworkflow.push_back(&PARAM_FORMAT_TYPE);
        predictmatchworkflow.push_back(&PARAM_REPORT_PAM);
        
        downloadgenomeworkflow = combineList(createdb, extractorfs);
        downloadgenomeworkflow = combineList(downloadgenomeworkflow, translatenucs);
        downloadgenomeworkflow = combineList(downloadgenomeworkflow, result2stats);
        downloadgenomeworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        downloadgenomeworkflow.push_back(&PARAM_REVERSE_SETDB);

        parsespacer.push_back(&PARAM_COMPRESSED);
        parsespacer.push_back(&PARAM_THREADS);
        parsespacer.push_back(&PARAM_V);

        filtermatchbyfdr.push_back(&PARAM_FDR_CUTOFF);
        filtermatchbyfdr.push_back(&PARAM_COMPRESSED);
        filtermatchbyfdr.push_back(&PARAM_THREADS);
        filtermatchbyfdr.push_back(&PARAM_V);

        truncatebesthits.push_back(&PARAM_COMPRESSED);
        truncatebesthits.push_back(&PARAM_THREADS);
        truncatebesthits.push_back(&PARAM_V);

        summarizeresults.push_back(&PARAM_FORMAT_TYPE);
        summarizeresults.push_back(&PARAM_COMPRESSED);
        summarizeresults.push_back(&PARAM_THREADS);
        summarizeresults.push_back(&PARAM_V);

        findpam.push_back(&PARAM_COMPRESSED);
        findpam.push_back(&PARAM_THREADS);
        findpam.push_back(&PARAM_V);
        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;
        extractorfsSpacer = 0;
        reverseSetDb = 1;
        fdrCutoff = 0.05;
        formatType = 1;
        reportPam = 1;
        
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
