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

    std::vector<MMseqsParameter*> predictmatchworkflow;
    std::vector<MMseqsParameter*> createsetdbworkflow;

    PARAMETER(PARAM_REVERSE_FRAGMENTS)
    int reverseFragments;

private:
    LocalParameters() : 
        Parameters(),
        PARAM_REVERSE_FRAGMENTS(PARAM_REVERSE_FRAGMENTS_ID,"--reverse-fragments", "Reverse AA Fragments", "reverse AA fragments to compute under null [0,1]", typeid(int), (void *) &reverseFragments, "^[0-1]{1}$"),
    
    {

        createsetdbworkflow = combineList(createdb, extractorfs);
        createsetdbworkflow = combineList(createsetdbworkflow, translatenucs);
        createsetdbworkflow = combineList(createsetdbworkflow, result2stats);
        createsetdbworkflow.push_back(&PARAM_REVERSE_FRAGMENTS);

        predictmatchworkflow = combineList(predictmatchworkflow, searchworkflow);
        predictmatchworkflow = combineList(predictmatchworkflow, besthitbyset);
        predictmatchworkflow = combineList(predictmatchworkflow, combinepvalperset);

        // default value 0 means no reverse of AA fragments
        reverseFragments = 0;

        
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
