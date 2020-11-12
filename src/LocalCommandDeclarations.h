#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int easypredict(int argc, const char **argv, const Command& command);
extern int predictmatch(int argc, const char **argv, const Command& command);
extern int createsetdb(int argc, const char **argv, const Command& command);
extern int parsespacer(int argc, const char **argv, const Command& command);
extern int filtermatchbyfdr(int argc, const char **argv, const Command& command);
extern int findpam(int argc, const char **argv, const Command& command);
extern int truncatebesthits(int argc, const char **argv, const Command& command);
extern int summarizeresults(int argc, const char **argv, const Command& command);
extern int downloadgenome(int argc, const char **argv, const Command& command);
extern int combinescore(int argc, const char **argv, const Command& command);
extern int combineprotnuclaln(int argc, const char **argv, const Command& command);
extern int empiricalpval(int argc, const char **argv, const Command& command);

#endif