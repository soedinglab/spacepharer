#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "MemoryMapped.h"
#include "KSeqWrapper.h"

#include <climits>

enum Format {
    PILERCR,
    CRT,
    CRISPRFINDER,
    MINCED,
    CRISPRDETECT,
    FASTA
};

const char *formatNames[] = {
    "PILERCR", "CRT", "CRISPRFINDER", "MINCED", "CRISPRDETECT", "FASTA"
};

Format detectFileType(char* data) {
    //check first 5 char for format type
    std::string line(data, 5);
    if (line == "piler") {
        return PILERCR;
    } else if (line == "ORGAN") {
        return CRT;
    } else if (line == "#####") {
        return CRISPRFINDER;
    } else if (line == "Array") {
        return CRISPRDETECT;
    } else if (line == "Seque") {
        return MINCED;
    } else {
        return FASTA;
    }
}

std::string getCurrentLine(char* data, size_t pos = 0) {
    std::string line;
    while (data[pos] != '\n' && data[pos] != '\0') {
        line.push_back(data[pos]);
        pos++;
    }
    return line;
}

bool isNucl(const std::string& sequence) {
    for (size_t i = 0; i < sequence.size(); i++) {
        if (sequence[i] != 'A' && sequence[i] != 'T' && sequence[i] != 'C' && sequence[i] != 'G') {
            return false;
        }
    }
    return true;
}

int parsespacer(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_VARIADIC, 0);

    std::string outDb = par.filenames.back();
    par.filenames.pop_back();

    DBWriter writer(outDb.c_str(), (outDb + ".index").c_str(), 1, par.compressed, Parameters::DBTYPE_NUCLEOTIDES);
    writer.open();

    DBWriter headers((outDb + "_h").c_str(), (outDb + "_h.index").c_str(), 1, par.compressed, Parameters::DBTYPE_GENERIC_DB);
    headers.open();

    FILE *lookup = FileUtil::openAndDelete((outDb + ".lookup").c_str(), "w");
    FILE *source = FileUtil::openAndDelete((outDb + ".source").c_str(), "w");

    Debug::Progress progress;
    size_t key = 0;
    for (size_t i = 0; i < par.filenames.size(); ++i) {
        std::string &file = par.filenames[i];
        std::string filename = FileUtil::baseName(file);
        fprintf(source, "%zu\t%s\n", i, filename.c_str());

        std::string sequence = "";
        sequence.reserve(par.maxSeqLen);

        std::string arrayHeader;
        std::string header;
        header.reserve(10000);
        std::string accession;

        unsigned int spacerStartPos;
        unsigned int spacerEndPos;
        size_t minSpacerLen = 20;

        Format type;

        size_t entry = 0;
        size_t arrayEntry = 0;
        size_t headerEntry = 0;
        size_t arrayNum = 0;
        size_t spacerNum = 0;
        bool isArrayReverse = false;

        MemoryMapped input(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        if (!input.isValid()){
            Debug(Debug::ERROR) << "Can map open file " << file << "\n";
            EXIT(EXIT_FAILURE);
        }
        char *data = (char *) input.getData();
        while (*data != '\0') {
            // pointer to the beginning of any line

            // if (entry == 0) {
            //     type = detectFileType(data);
            //     if (type != CRT||CRISPRDetect) {
            //         continue;
            //     }
            // }
            if (entry == 0) {
                type = detectFileType(data);
                Debug(Debug::INFO) << "Detected input file " << file << " is of type " << formatNames[type] << "\n";
                if (type == CRT) {
                    arrayHeader = getCurrentLine(data, 11);
                    accession = Util::parseFastaHeader(arrayHeader.c_str());
                }
                if (type == MINCED) {
                    arrayHeader = getCurrentLine(data);
                    std::vector<std::string> values = Util::split(arrayHeader, "'");
                    accession = Util::parseFastaHeader(values[1].c_str());
                }
            }

            switch (type) {
                case PILERCR:
                    if (*data == 'A') {
                        arrayEntry = entry;
                        arrayNum++;
                    }
                    if (*data == '>' && (entry == arrayEntry + 1)) {
                        arrayHeader = getCurrentLine(data, 1);
                        accession = Util::parseFastaHeader(arrayHeader.c_str());
                        headerEntry = entry;
                    }
                    if (*data == ' ' && (entry == headerEntry + 4)) {
                        while (*data != '=') {
                            std::string line = getCurrentLine(data);
                            std::vector<std::string> values = Util::split(line, " ");
                            // magic?
                            sequence = Util::removeWhiteSpace(sequence = values.back().c_str());
                            if (sequence.size() >= minSpacerLen && isNucl(sequence)) {
                                spacerNum++;
                                //position + repeat length
                                spacerStartPos = (unsigned int) strtoul(values[0].c_str(), NULL, 10) + (unsigned int) strtoul(values[1].c_str(), NULL, 10);
                                spacerEndPos = spacerStartPos + sequence.size();
                                //"accession_Array_#_spacer_#_spacerStartPos_spacerEndPos_spacerLen"
                                header.append(accession.c_str());
                                header.append("_Array_");
                                header.append(SSTR(arrayNum));
                                header.append("_spacer_");
                                header.append(SSTR(spacerNum));
                                header.append("_");
                                header.append(SSTR(spacerStartPos));
                                header.append("_");
                                header.append(SSTR(spacerEndPos));
                                header.append("_");
                                header.append(SSTR(sequence.size()));
                                fprintf(lookup, "%zu\t%s\t%zu\n", key, header.c_str(), i);
                                sequence += "\n";
                                header += "\n";
                                writer.writeData(sequence.c_str(), sequence.size(), key, 0);
                                headers.writeData(header.c_str(), header.size(), key, 0);
                                key++;
                                progress.updateProgress();
                                // TODO shouldnt the clear be here??
                                sequence.clear();
                                header.clear();
                            }
                            data = Util::skipLine(data);
                            entry++;

                        }
                    }
                    data = Util::skipLine(data);
                    entry++;
                    break;
                case MINCED:
                case CRT:
                    if (*data == 'C') {
                        arrayEntry = entry;
                        arrayNum++;
                    }
                    if (entry > 3 && entry == arrayEntry + 3) {
                        while (*data != '-') {
                            std::string line = getCurrentLine(data);
                            std::vector<std::string> values = Util::split(line, "\t");
                            //avoid last row empty spacer column
                            if (values.size() == 3) {
                                sequence = Util::removeWhiteSpace(values[2].c_str());
                                if (sequence.size() >= minSpacerLen && isNucl(sequence)) {
                                    spacerNum++;
                                    spacerStartPos = (unsigned int) strtoul(values[0].c_str(), NULL, 10) + values[1].size();
                                    spacerEndPos = spacerStartPos + sequence.size();
                                    //"accession_Array_#_spacer_#_spacerStartPos_spacerEndPos_spacerLen"
                                    header.append(accession.c_str());
                                    header.append("_Array_");
                                    header.append(SSTR(arrayNum));
                                    header.append("_spacer_");
                                    header.append(SSTR(spacerNum));
                                    header.append("_");
                                    header.append(SSTR(spacerStartPos));
                                    header.append("_");
                                    header.append(SSTR(spacerEndPos));
                                    header.append("_");
                                    header.append(SSTR(sequence.size()));
                                    fprintf(lookup, "%zu\t%s\t%zu\n", key, header.c_str(), i);
                                    sequence += std::string("\n");
                                    header += std::string("\n");
                                    writer.writeData(sequence.c_str(), sequence.size(), key, 0);
                                    headers.writeData(header.c_str(), header.size(), key, 0);
                                    key++;
                                    progress.updateProgress();
                                    sequence.clear();
                                    header.clear();
                                }
                            }
                            data = Util::skipLine(data);
                            entry++;
                        }
                    }
                    data = Util::skipLine(data);
                    entry++;
                    break;

                case CRISPRDETECT:
                    if (*data == 'A') {
                        arrayEntry = entry;
                        arrayNum++;
                    }
                    if (*data == '>' && (entry == arrayEntry + 1)) {
                        arrayHeader = getCurrentLine(data, 1);
                        std::vector<std::string> headerValues = Util::split(arrayHeader, "\t");
                        isArrayReverse = ((headerValues[1].find("Reverse") != std::string::npos) ? true : false);
                        accession = Util::parseFastaHeader(arrayHeader.c_str());
                        headerEntry = entry;
                    }
                    if (*data == ' ' && (entry == headerEntry + 4)) {
                        while (*data != '=') {
                            std::string line = getCurrentLine(data);
                            std::vector<std::string> values = Util::split(line, "\t");
                            sequence = Util::removeWhiteSpace(values[5].c_str());
                            if (sequence.size() >= minSpacerLen && isNucl(sequence)) {
                                spacerNum++;
                                if (isArrayReverse) {
                                    spacerStartPos = (unsigned int) strtoul(values[0].c_str(), NULL, 10) - (unsigned int) strtoul(values[1].c_str(), NULL, 10);
                                    spacerEndPos = spacerStartPos - sequence.size();
                                } else {
                                    spacerStartPos = (unsigned int) strtoul(values[0].c_str(), NULL, 10) + (unsigned int) strtoul(values[1].c_str(), NULL, 10);
                                    spacerEndPos = spacerStartPos + sequence.size();
                                }
                                //"accession_Array_#_spacer_#_spacerStartPos_spacerEndPos_spacerLen"
                                header.append(accession.c_str());
                                header.append("_Array_");
                                header.append(SSTR(arrayNum));
                                header.append("_spacer_");
                                header.append(SSTR(spacerNum));
                                header.append("_");
                                header.append(SSTR(spacerStartPos));
                                header.append("_");
                                header.append(SSTR(spacerEndPos));
                                header.append("_");
                                header.append(SSTR(sequence.size()));
                                fprintf(lookup, "%zu\t%s\t%zu\n", key, header.c_str(), i);
                                sequence += std::string("\n");
                                header += std::string("\n");
                                writer.writeData(sequence.c_str(), sequence.size(), key, 0);
                                headers.writeData(header.c_str(), header.size(), key, 0);
                                key++;
                                progress.updateProgress();
                                sequence.clear();
                                header.clear();
                            }
                            data = Util::skipLine(data);
                            entry++;
                        }
                    }
                    data = Util::skipLine(data);
                    entry++;
                    break;

                case CRISPRFINDER:
                    Debug(Debug::ERROR) << "SpacePHARER currently does not support CRISPRfinder format.\n";
                    EXIT(EXIT_FAILURE);
                    break;

                default:
                    input.close();
                    KSeqWrapper* kseq = KSeqFactory(file.c_str());
                    while (kseq->ReadEntry()) {
                        progress.updateProgress();
                        const KSeqWrapper::KSeqEntry &e = kseq->entry;
                        if (e.name.l == 0) {
                            Debug(Debug::ERROR) << "Fasta entry " << key << " is invalid.\n";
                            EXIT(EXIT_FAILURE);
                        }

                        header.append(e.name.s, e.name.l);
                        if (e.comment.l > 0) {
                            header.append(" ", 1);
                            header.append(e.comment.s, e.comment.l);
                        }
                        header.push_back('\n');

                        std::string headerId = Util::parseFastaHeader(header.c_str());
                        if (headerId.empty()) {
                            // An identifier is necessary for these two cases, so we should just give up
                            Debug(Debug::WARNING) << "Cannot extract identifier from entry " << key << "\n";
                        }

                        headers.writeData(header.c_str(), header.length(), key, 0);
                        writer.writeStart(0);
                        writer.writeAdd(e.sequence.s, e.sequence.l, 0);
                        const char newline = '\n';
                        writer.writeAdd(&newline, 1, 0);
                        writer.writeEnd(key, 0, true);
                        fprintf(lookup, "%zu\t%s\t%zu\n", key, headerId.c_str(), i);
                        key++;
                        progress.updateProgress();
                        header.clear();
                    }
                    delete kseq;
                // process next file instead of continuing with the loop
                char null = '\0';
                data = &null;
                break;
            }
        }
        input.close();
    }
    fclose(lookup);
    fclose(source);
    headers.close(true);
    writer.close(true);

    return EXIT_SUCCESS;
}
