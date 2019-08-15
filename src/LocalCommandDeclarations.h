#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int predictmatch(int argc, const char **argv, const Command& command);
extern int createsetdb(int argc, const char **argv, const Command& command);
extern int spacerparser(int argc, const char **argv, const Command& command);

#endif
