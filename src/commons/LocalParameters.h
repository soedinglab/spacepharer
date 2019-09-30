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
    std::vector<MMseqsParameter*> truncatebesthits;    
    std::vector<MMseqsParameter*> predictmatchworkflow;
    std::vector<MMseqsParameter*> createsetdbworkflow;
    std::vector<MMseqsParameter*> downloadgbphageworkflow;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    PARAMETER(PARAM_REVERSE_SETDB)
    PARAMETER(PARAM_EXTRACTORF_SPACER)
    PARAMETER(PARAM_FDR_CUTOFF)

    int reverseFragments;
    int extractorfsSpacer;
    int reverseSetDb;
    float fdrCutoff;
    
private:
    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
        PARAM_REVERSE_SETDB(PARAM_REVERSE_SETDB_ID,"--reverse-setDB", "Create setDB with reversed fragments", "Create additional setDB with reversed fragments to compute under null [0,1]", typeid(int), (void *) &reverseSetDb, "^[0-1]{1}$"),
        PARAM_EXTRACTORF_SPACER(PARAM_EXTRACTORF_SPACER_ID,"--extractorf-spacer", "Extract orfs from spacers", "change parameter settings for extractorfs when createsetdb for spacer db [0,1]", typeid(int), (void *) &extractorfsSpacer, "^[0-1]{1}$"),
        PARAM_FDR_CUTOFF(PARAM_EXTRACTORF_SPACER_ID,"--fdr", "FDR cutoff", "FDR cutoff for filtering matches[0.0, 1.0]", typeid(float), (void *) &fdrCutoff, "^0(\\.[0-9]+)?|1(\\.0+)?$")

    {

        createsetdbworkflow = combineList(createdb, extractorfs);
        createsetdbworkflow = combineList(createsetdbworkflow, translatenucs);
        createsetdbworkflow = combineList(createsetdbworkflow, result2stats);
        createsetdbworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        createsetdbworkflow.push_back(&PARAM_EXTRACTORF_SPACER);
    
        predictmatchworkflow = combineList(predictmatchworkflow, searchworkflow);
        predictmatchworkflow = combineList(predictmatchworkflow, besthitbyset);
        predictmatchworkflow = combineList(predictmatchworkflow, combinepvalbyset);
        
        downloadgbphageworkflow = combineList(createdb, extractorfs);
        downloadgbphageworkflow = combineList(downloadgbphageworkflow, translatenucs);
        downloadgbphageworkflow = combineList(downloadgbphageworkflow, result2stats);
        downloadgbphageworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);
        downloadgbphageworkflow.push_back(&PARAM_REVERSE_SETDB);

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

        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;
        extractorfsSpacer = 0;
        reverseSetDb = 1;
        fdrCutoff = 0.05;
        
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
