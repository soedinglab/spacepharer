#include "PrefilteringIndexReader.h"
#include "DBWriter.h"
#include "Prefiltering.h"
#include "ExtendedSubstitutionMatrix.h"
#include "FileUtil.h"
#include "IndexBuilder.h"
#include "Parameters.h"

const char*  PrefilteringIndexReader::CURRENT_VERSION = "15";
unsigned int PrefilteringIndexReader::VERSION = 0;
unsigned int PrefilteringIndexReader::META = 1;
unsigned int PrefilteringIndexReader::SCOREMATRIXNAME = 2;
unsigned int PrefilteringIndexReader::SCOREMATRIX2MER = 3;
unsigned int PrefilteringIndexReader::SCOREMATRIX3MER = 4;
unsigned int PrefilteringIndexReader::DBR1INDEX = 5;
unsigned int PrefilteringIndexReader::DBR1DATA  = 6;
unsigned int PrefilteringIndexReader::DBR2INDEX = 7;
unsigned int PrefilteringIndexReader::DBR2DATA  = 8;
unsigned int PrefilteringIndexReader::ENTRIES = 9;
unsigned int PrefilteringIndexReader::ENTRIESOFFSETS = 10;
unsigned int PrefilteringIndexReader::ENTRIESGRIDSIZE = 11;
unsigned int PrefilteringIndexReader::ENTRIESNUM = 12;
unsigned int PrefilteringIndexReader::SEQCOUNT = 13;
unsigned int PrefilteringIndexReader::SEQINDEXDATA = 14;
unsigned int PrefilteringIndexReader::SEQINDEXDATASIZE = 15;
unsigned int PrefilteringIndexReader::SEQINDEXSEQOFFSET = 16;
unsigned int PrefilteringIndexReader::HDR1INDEX = 18;
unsigned int PrefilteringIndexReader::HDR1DATA = 19;
unsigned int PrefilteringIndexReader::HDR2INDEX = 20;
unsigned int PrefilteringIndexReader::HDR2DATA = 21;
unsigned int PrefilteringIndexReader::GENERATOR = 22;
unsigned int PrefilteringIndexReader::SPACEDPATTERN = 23;

extern const char* version;

bool PrefilteringIndexReader::checkIfIndexFile(DBReader<unsigned int>* reader) {
    char * version = reader->getDataByDBKey(VERSION, 0);
    if(version == NULL){
        return false;
    }
    return (strncmp(version, CURRENT_VERSION, strlen(CURRENT_VERSION)) == 0 ) ? true : false;
}

std::string PrefilteringIndexReader::indexName(const std::string &outDB) {
    std::string result(outDB);
    result.append(".idx");
    return result;
}

void PrefilteringIndexReader::createIndexFile(const std::string &outDB,
                                              DBReader<unsigned int> *dbr1, DBReader<unsigned int> *dbr2,
                                              DBReader<unsigned int> *hdbr1, DBReader<unsigned int> *hdbr2,
                                              BaseMatrix *subMat, int maxSeqLen,
                                              bool hasSpacedKmer, const std::string &spacedKmerPattern,
                                              bool compBiasCorrection, int alphabetSize, int kmerSize,
                                              int maskMode, int maskLowerCase, int kmerThr, int splits) {
    DBWriter writer(outDB.c_str(), std::string(outDB).append(".index").c_str(), splits, Parameters::WRITER_ASCII_MODE, Parameters::DBTYPE_INDEX_DB);
    writer.open();

    Debug(Debug::INFO) << "Write VERSION (" << VERSION << ")\n";
    writer.writeData((char *) CURRENT_VERSION, strlen(CURRENT_VERSION) * sizeof(char), VERSION, 0);
    writer.alignToPageSize();

    Debug(Debug::INFO) << "Write META (" << META << ")\n";
    const int biasCorr = compBiasCorrection ? 1 : 0;
    const int mask = maskMode > 0;
    const int spacedKmer = (hasSpacedKmer) ? 1 : 0;
    const int headers1 = (hdbr1 != NULL) ? 1 : 0;
    const int headers2 = (hdbr2 != NULL) ? 1 : 0;
    const int seqType = dbr1->getDbtype();
    const int srcSeqType = (dbr2 !=NULL) ? dbr2->getDbtype() : seqType;
    int metadata[] = {maxSeqLen, kmerSize, biasCorr, alphabetSize, mask, spacedKmer, kmerThr, seqType, srcSeqType, headers1, headers2, splits};
    char *metadataptr = (char *) &metadata;
    writer.writeData(metadataptr, sizeof(metadata), META, 0);
    writer.alignToPageSize();

    if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE) == false &&
        Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_PROFILE_STATE_SEQ) == false) {
        int alphabetSize = subMat->alphabetSize;
        subMat->alphabetSize = subMat->alphabetSize-1;
        ScoreMatrix s3 = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
        ScoreMatrix s2 = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        subMat->alphabetSize = alphabetSize;

        char* serialized3mer = ScoreMatrix::serialize(s3);
        Debug(Debug::INFO) << "Write SCOREMATRIX3MER (" << SCOREMATRIX3MER << ")\n";
        writer.writeData(serialized3mer, ScoreMatrix::size(s3), SCOREMATRIX3MER, 0);
        writer.alignToPageSize();
        ExtendedSubstitutionMatrix::freeScoreMatrix(s3);
        free(serialized3mer);

        char* serialized2mer = ScoreMatrix::serialize(s2);
        Debug(Debug::INFO) << "Write SCOREMATRIX2MER (" << SCOREMATRIX2MER << ")\n";
        writer.writeData(serialized2mer, ScoreMatrix::size(s2), SCOREMATRIX2MER, 0);
        writer.alignToPageSize();
        ExtendedSubstitutionMatrix::freeScoreMatrix(s2);
        free(serialized2mer);
    }

    Debug(Debug::INFO) << "Write SCOREMATRIXNAME (" << SCOREMATRIXNAME << ")\n";
    char* subData = BaseMatrix::serialize(subMat);
    writer.writeData(subData, BaseMatrix::memorySize(subMat), SCOREMATRIXNAME, 0);
    writer.alignToPageSize();
    free(subData);

    if (spacedKmerPattern.empty() != false) {
        Debug(Debug::INFO) << "Write SPACEDPATTERN (" << SPACEDPATTERN << ")\n";
        writer.writeData(spacedKmerPattern.c_str(), spacedKmerPattern.length(), SPACEDPATTERN, 0);
        writer.alignToPageSize();
    }

    Debug(Debug::INFO) << "Write DBR1INDEX (" << DBR1INDEX << ")\n";
    char* data = DBReader<unsigned int>::serialize(*dbr1);
    size_t offsetIndex = writer.getOffset(0);
    writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*dbr1), DBR1INDEX, 0);
    writer.alignToPageSize();

    Debug(Debug::INFO) << "Write DBR1DATA (" << DBR1DATA << ")\n";
    size_t offsetData = writer.getOffset(0);
    writer.writeStart(0);
    for(size_t fileIdx = 0; fileIdx < dbr1->getDataFileCnt(); fileIdx++) {
        writer.writeAdd(dbr1->getDataForFile(fileIdx), dbr1->getDataSizeForFile(fileIdx), 0);
    }
    writer.writeEnd(DBR1DATA, 0);
    writer.alignToPageSize();
    free(data);

    if (dbr2 == NULL) {
        writer.writeIndexEntry(DBR2INDEX, offsetIndex, DBReader<unsigned int>::indexMemorySize(*dbr1)+1, 0);
        writer.writeIndexEntry(DBR2DATA,  offsetData,  dbr1->getTotalDataSize()+1, 0);
    } else {
        Debug(Debug::INFO) << "Write DBR2INDEX (" << DBR2INDEX << ")\n";
        data = DBReader<unsigned int>::serialize(*dbr2);
        writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*dbr2), DBR2INDEX, 0);
        writer.alignToPageSize();
        Debug(Debug::INFO) << "Write DBR2DATA (" << DBR2DATA << ")\n";
        writer.writeStart(0);
        for(size_t fileIdx = 0; fileIdx < dbr2->getDataFileCnt(); fileIdx++) {
            writer.writeAdd(dbr2->getDataForFile(fileIdx), dbr2->getDataSizeForFile(fileIdx), 0);
        }
        writer.writeEnd(DBR2DATA, 0);
        writer.alignToPageSize();
        free(data);
    }

    if (hdbr1 != NULL) {
        Debug(Debug::INFO) << "Write HDR1INDEX (" << HDR1INDEX << ")\n";
        data = DBReader<unsigned int>::serialize(*hdbr1);
        size_t offsetIndex = writer.getOffset(0);
        writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*hdbr1), HDR1INDEX, 0);
        writer.alignToPageSize();

        Debug(Debug::INFO) << "Write HDR1DATA (" << HDR1DATA << ")\n";
        size_t offsetData = writer.getOffset(0);
        writer.writeStart(0);
        for(size_t fileIdx = 0; fileIdx < hdbr1->getDataFileCnt(); fileIdx++) {
            writer.writeAdd(hdbr1->getDataForFile(fileIdx), hdbr1->getDataSizeForFile(fileIdx), 0);
        }
        writer.writeEnd(HDR1DATA, 0);
        writer.alignToPageSize();
        free(data);
        if (hdbr2 == NULL) {
            writer.writeIndexEntry(HDR2INDEX, offsetIndex, DBReader<unsigned int>::indexMemorySize(*hdbr1)+1, 0);
            writer.writeIndexEntry(HDR2DATA,  offsetData, hdbr1->getTotalDataSize()+1, 0);
        }
    }
    if (hdbr2 != NULL) {
        Debug(Debug::INFO) << "Write HDR2INDEX (" << HDR2INDEX << ")\n";
        data = DBReader<unsigned int>::serialize(*hdbr2);
        writer.writeData(data, DBReader<unsigned int>::indexMemorySize(*hdbr2), HDR2INDEX, 0);
        writer.alignToPageSize();
        Debug(Debug::INFO) << "Write HDR2DATA (" << HDR2DATA << ")\n";
        writer.writeStart(0);
        for(size_t fileIdx = 0; fileIdx < hdbr2->getDataFileCnt(); fileIdx++) {
            writer.writeAdd(hdbr2->getDataForFile(fileIdx), hdbr2->getDataSizeForFile(fileIdx), 0);
        }
        writer.writeEnd(HDR2DATA, 0);
        writer.alignToPageSize();
        free(data);
    }
    Debug(Debug::INFO) << "Write GENERATOR (" << GENERATOR << ")\n";
    writer.writeData(version, strlen(version), GENERATOR, 0);
    writer.alignToPageSize();

    Sequence seq(maxSeqLen, seqType, subMat, kmerSize, hasSpacedKmer, compBiasCorrection, true, spacedKmerPattern);
    // remove x (not needed in index)
    const int adjustAlphabetSize =
            (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_NUCLEOTIDES) || Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_AMINO_ACIDS))
                ? alphabetSize -1: alphabetSize;

    for (int s = 0; s < splits; s++) {
        size_t dbFrom = 0;
        size_t dbSize = 0;
        Util::decomposeDomainByAminoAcid(dbr1->getDataSize(), dbr1->getSeqLens(), dbr1->getSize(), s, splits, &dbFrom, &dbSize);
        if (dbSize == 0) {
            continue;
        }

        IndexTable indexTable(adjustAlphabetSize, kmerSize, false);
        SequenceLookup *sequenceLookup = NULL;
        IndexBuilder::fillDatabase(&indexTable,
                                   (maskMode == 1 || maskLowerCase == 1) ? &sequenceLookup : NULL,
                                   (maskMode == 0 ) ? &sequenceLookup : NULL,
                                   *subMat, &seq, dbr1, dbFrom, dbSize, kmerThr, maskMode, maskLowerCase);
        indexTable.printStatistics(subMat->int2aa);

        if (sequenceLookup == NULL) {
            Debug(Debug::ERROR) << "Invalid mask mode. No sequence lookup created!\n";
            EXIT(EXIT_FAILURE);
        }

        // save the entries
        unsigned int keyOffset = 1000 * s;
        Debug(Debug::INFO) << "Write ENTRIES (" << (keyOffset + ENTRIES) << ")\n";
        char *entries = (char *) indexTable.getEntries();
        size_t entriesSize = indexTable.getTableEntriesNum() * indexTable.getSizeOfEntry();
        writer.writeData(entries, entriesSize, (keyOffset + ENTRIES), s);
        writer.alignToPageSize(s);

        // save the size
        Debug(Debug::INFO) << "Write ENTRIESOFFSETS (" << (keyOffset + ENTRIESOFFSETS) << ")\n";
        char *offsets = (char*)indexTable.getOffsets();
        size_t offsetsSize = (indexTable.getTableSize() + 1) * sizeof(size_t);
        writer.writeData(offsets, offsetsSize, (keyOffset + ENTRIESOFFSETS), s);
        writer.alignToPageSize(s);
        indexTable.deleteEntries();

        Debug(Debug::INFO) << "Write SEQINDEXDATASIZE (" << (keyOffset + SEQINDEXDATASIZE) << ")\n";
        int64_t seqindexDataSize = sequenceLookup->getDataSize();
        char *seqindexDataSizePtr = (char *) &seqindexDataSize;
        writer.writeData(seqindexDataSizePtr, 1 * sizeof(int64_t), (keyOffset + SEQINDEXDATASIZE), s);
        writer.alignToPageSize(s);

        size_t *sequenceOffsets = sequenceLookup->getOffsets();
        size_t sequenceCount = sequenceLookup->getSequenceCount();
        Debug(Debug::INFO) << "Write SEQINDEXSEQOFFSET (" << (keyOffset + SEQINDEXSEQOFFSET) << ")\n";
        writer.writeData((char *) sequenceOffsets, (sequenceCount + 1) * sizeof(size_t), (keyOffset + SEQINDEXSEQOFFSET), s);
        writer.alignToPageSize(s);

        Debug(Debug::INFO) << "Write SEQINDEXDATA (" << (keyOffset + SEQINDEXDATA) << ")\n";
        writer.writeData(sequenceLookup->getData(), (sequenceLookup->getDataSize() + 1) * sizeof(char), (keyOffset + SEQINDEXDATA), s);
        writer.alignToPageSize(s);
        delete sequenceLookup;

        // ENTRIESNUM
        Debug(Debug::INFO) << "Write ENTRIESNUM (" << (keyOffset + ENTRIESNUM) << ")\n";
        uint64_t entriesNum = indexTable.getTableEntriesNum();
        char *entriesNumPtr = (char *) &entriesNum;
        writer.writeData(entriesNumPtr, 1 * sizeof(uint64_t), (keyOffset + ENTRIESNUM), s);
        writer.alignToPageSize(s);

        // SEQCOUNT
        Debug(Debug::INFO) << "Write SEQCOUNT (" << (keyOffset + SEQCOUNT) << ")\n";
        size_t tablesize = indexTable.getSize();
        char *tablesizePtr = (char *) &tablesize;
        writer.writeData(tablesizePtr, 1 * sizeof(size_t), (keyOffset + SEQCOUNT), s);
        writer.alignToPageSize(s);
    }

    writer.close(false);
}

DBReader<unsigned int> *PrefilteringIndexReader::openNewHeaderReader(DBReader<unsigned int>*dbr, unsigned int dataIdx, unsigned int indexIdx, int threads,  bool touchIndex, bool touchData) {
    size_t indexId = dbr->getId(indexIdx);
    char *indexData = dbr->getData(indexId, 0);
    if (touchIndex) {
        dbr->touchData(indexId);
    }

    size_t dataId = dbr->getId(dataIdx);
    char *data = dbr->getData(dataId, 0);

    size_t currDataOffset = dbr->getOffset(dataId);
    size_t nextDataOffset = dbr->findNextOffsetid(dataId);
    size_t dataSize = nextDataOffset-currDataOffset;

    if (touchData) {
        dbr->touchData(dataId);
    }

    DBReader<unsigned int> *reader = DBReader<unsigned int>::unserialize(indexData, threads);
    reader->open(DBReader<unsigned int>::NOSORT);
    reader->setData(data, dataSize);
    reader->setMode(DBReader<unsigned int>::USE_DATA);
    return reader;
}

DBReader<unsigned int> *PrefilteringIndexReader::openNewReader(DBReader<unsigned int>*dbr, unsigned int dataIdx, unsigned int indexIdx, bool includeData, int threads, bool touchIndex, bool touchData) {
    size_t id = dbr->getId(indexIdx);
    char *data = dbr->getDataUncompressed(id);
    if (touchIndex) {
        dbr->touchData(id);
    }

    if (includeData) {
        id = dbr->getId(dataIdx);
        if (id == UINT_MAX) {
            return NULL;
        }
        if (touchData) {
            dbr->touchData(id);
        }

        DBReader<unsigned int> *reader = DBReader<unsigned int>::unserialize(data, threads);
        reader->open(DBReader<unsigned int>::NOSORT);
        size_t currDataOffset = dbr->getOffset(id);
        size_t nextDataOffset = dbr->findNextOffsetid(id);
        size_t dataSize = nextDataOffset-currDataOffset;
        reader->setData(dbr->getDataUncompressed(id), dataSize);
        reader->setMode(DBReader<unsigned int>::USE_DATA);
        return reader;
    }

    DBReader<unsigned int> *reader = DBReader<unsigned int>::unserialize(data, threads);
    reader->open(DBReader<unsigned int>::NOSORT);
    return reader;
}

SequenceLookup *PrefilteringIndexReader::getSequenceLookup(unsigned int split, DBReader<unsigned int> *dbr, int preloadMode) {
    PrefilteringIndexData data = getMetadata(dbr);
    if (split >= (unsigned int)data.splits) {
        Debug(Debug::ERROR) << "Invalid split " << split << " out of " << data.splits << " chosen.\n";
        EXIT(EXIT_FAILURE);
    }

    unsigned int splitOffset = split * 1000;

    size_t id = dbr->getId(splitOffset + SEQINDEXDATA);
    if (id == UINT_MAX) {
        return NULL;
    }

    char * seqData = dbr->getDataUncompressed(id);

    size_t seqOffsetsId = dbr->getId(splitOffset + SEQINDEXSEQOFFSET);
    char * seqOffsetsData = dbr->getDataUncompressed(seqOffsetsId);

    size_t seqDataSizeId = dbr->getId(splitOffset + SEQINDEXDATASIZE);
    int64_t seqDataSize = *((int64_t *)dbr->getDataUncompressed(seqDataSizeId));

    size_t sequenceCountId = dbr->getId(splitOffset + SEQCOUNT);
    size_t sequenceCount = *((size_t *)dbr->getDataUncompressed(sequenceCountId));

    if (preloadMode == Parameters::PRELOAD_MODE_FREAD) {
        SequenceLookup *sequenceLookup = new SequenceLookup(sequenceCount, seqDataSize);
        sequenceLookup->initLookupByExternalDataCopy(seqData, seqDataSize, (size_t *) seqOffsetsData);
        return sequenceLookup;
    }

    if (preloadMode == Parameters::PRELOAD_MODE_MMAP_TOUCH) {
        dbr->touchData(id);
        dbr->touchData(seqOffsetsId);
    }

    SequenceLookup *sequenceLookup = new SequenceLookup(sequenceCount);
    sequenceLookup->initLookupByExternalData(seqData, seqDataSize, (size_t *) seqOffsetsData);
    return sequenceLookup;
}

IndexTable *PrefilteringIndexReader::getIndexTable(unsigned int split, DBReader<unsigned int> *dbr, int preloadMode) {
    PrefilteringIndexData data = getMetadata(dbr);
    if (split >= (unsigned int)data.splits) {
        Debug(Debug::ERROR) << "Invalid split " << split << " out of " << data.splits << " chosen.\n";
        EXIT(EXIT_FAILURE);
    }

    unsigned int splitOffset = split * 1000;
    size_t entriesNumId = dbr->getId(splitOffset + ENTRIESNUM);
    int64_t entriesNum = *((int64_t *)dbr->getDataUncompressed(entriesNumId));
    size_t sequenceCountId = dbr->getId(splitOffset +SEQCOUNT);
    size_t sequenceCount = *((size_t *)dbr->getDataUncompressed(sequenceCountId));

    size_t entriesDataId = dbr->getId(splitOffset + ENTRIES);
    char *entriesData = dbr->getDataUncompressed(entriesDataId);

    size_t entriesOffsetsDataId = dbr->getId(splitOffset + ENTRIESOFFSETS);
    char *entriesOffsetsData = dbr->getDataUncompressed(entriesOffsetsDataId);

    int adjustAlphabetSize;
    if (Parameters::isEqualDbtype(data.seqType, Parameters::DBTYPE_NUCLEOTIDES) || Parameters::isEqualDbtype(data.seqType, Parameters::DBTYPE_AMINO_ACIDS)) {
        adjustAlphabetSize = data.alphabetSize - 1;
    } else {
        adjustAlphabetSize = data.alphabetSize;
    }

    if (preloadMode == Parameters::PRELOAD_MODE_FREAD) {
        IndexTable* table = new IndexTable(adjustAlphabetSize, data.kmerSize, true);
        table->initTableByExternalDataCopy(sequenceCount, entriesNum, (IndexEntryLocal*) entriesData, (size_t *)entriesOffsetsData);
        return table;
    }

    if (preloadMode == Parameters::PRELOAD_MODE_MMAP_TOUCH) {
        dbr->touchData(entriesNumId);
        dbr->touchData(sequenceCountId);
        dbr->touchData(entriesDataId);
        dbr->touchData(entriesOffsetsDataId);
    }

    IndexTable* table = new IndexTable(adjustAlphabetSize, data.kmerSize, true);
    table->initTableByExternalData(sequenceCount, entriesNum, (IndexEntryLocal*) entriesData, (size_t *)entriesOffsetsData);
    return table;
}

void PrefilteringIndexReader::printSummary(DBReader<unsigned int> *dbr) {
    Debug(Debug::INFO) << "Index version: " << dbr->getDataByDBKey(VERSION, 0) << "\n";

    size_t id;
    if ((id = dbr->getId(GENERATOR)) != UINT_MAX) {
        Debug(Debug::INFO)
                       << "Generated by:  " << dbr->getDataUncompressed(id) << "\n";
    }

    char * subMatData = dbr->getDataByDBKey(SCOREMATRIXNAME, 0);
    size_t pos = 0;
    while(subMatData[pos]!='\0') {
        if (subMatData[pos] == '.'
            && subMatData[pos + 1] == 'o'
            && subMatData[pos + 2] == 'u'
            && subMatData[pos + 3] == 't'
            && subMatData[pos + 4] == ':') {
            break;
        }
        pos++;
    }
    Debug(Debug::INFO) << "ScoreMatrix:  " << std::string(subMatData, pos+4) << "\n";
}

void PrefilteringIndexReader::printMeta(int *metadata_tmp) {
    Debug(Debug::INFO) << "MaxSeqLength: " << metadata_tmp[0] << "\n";
    Debug(Debug::INFO) << "KmerSize:     " << metadata_tmp[1] << "\n";
    Debug(Debug::INFO) << "CompBiasCorr: " << metadata_tmp[2] << "\n";
    Debug(Debug::INFO) << "AlphabetSize: " << metadata_tmp[3] << "\n";
    Debug(Debug::INFO) << "Masked:       " << metadata_tmp[4] << "\n";
    Debug(Debug::INFO) << "Spaced:       " << metadata_tmp[5] << "\n";
    Debug(Debug::INFO) << "KmerScore:    " << metadata_tmp[6] << "\n";
    Debug(Debug::INFO) << "SequenceType: " << Parameters::getDbTypeName(metadata_tmp[7]) << "\n";
    Debug(Debug::INFO) << "SourcSeqType: " << Parameters::getDbTypeName(metadata_tmp[8]) << "\n";
    Debug(Debug::INFO) << "Headers1:     " << metadata_tmp[9] << "\n";
    Debug(Debug::INFO) << "Headers2:     " << metadata_tmp[10] << "\n";
    // Keep compatible to index version 15
    Debug(Debug::INFO) << "Splits:       " << (metadata_tmp[11] == 0 ? 1 : metadata_tmp[11]) << "\n";
}

PrefilteringIndexData PrefilteringIndexReader::getMetadata(DBReader<unsigned int> *dbr) {
    int *meta = (int *)dbr->getDataByDBKey(META, 0);

    PrefilteringIndexData data;
    data.maxSeqLength = meta[0];
    data.kmerSize = meta[1];
    data.compBiasCorr = meta[2];
    data.alphabetSize = meta[3];
    data.mask = meta[4];
    data.spacedKmer = meta[5];
    data.kmerThr = meta[6];
    data.seqType = meta[7];
    data.srcSeqType = meta[8];
    data.headers1 = meta[9];
    data.headers2 = meta[10];
    // Keep compatible to index version 15, where meta[11] would have been zero due to the alignment padding
    data.splits = meta[11] == 0 ? 1 : meta[11];

    return data;
}

std::string PrefilteringIndexReader::getSubstitutionMatrixName(DBReader<unsigned int> *dbr) {
    unsigned int key = dbr->getDbKey(SCOREMATRIXNAME);
    if (key == UINT_MAX) {
        return "";
    }
    const char *data = dbr->getData(key, 0);
    size_t len = dbr->getSeqLens(key) - 1;
    std::string matrixName;
    bool found = false;
    for (size_t pos = 0; pos < std::max(len, (size_t)4) - 4 && found == false; pos++) {
        if (data[pos] == '.'
            && data[pos + 1] == 'o'
            && data[pos + 2] == 'u'
            && data[pos + 3] == 't'
            && data[pos + 4] == ':') {
            matrixName = std::string(data, pos + 4);
            found = true;
        }
    }
    if (found == false) {
        matrixName = std::string(data);
    }
    return matrixName;
}

std::string PrefilteringIndexReader::getSubstitutionMatrix(DBReader<unsigned int> *dbr) {
    return std::string(dbr->getDataByDBKey(SCOREMATRIXNAME, 0));
}

std::string PrefilteringIndexReader::getSpacedPattern(DBReader<unsigned int> *dbr) {
    size_t id = dbr->getId(SPACEDPATTERN);
    if (id == UINT_MAX) {
        return "";
    }
    return std::string(dbr->getDataUncompressed(id));
}

ScoreMatrix PrefilteringIndexReader::get2MerScoreMatrix(DBReader<unsigned int> *dbr, int preloadMode) {
    size_t id = dbr->getId(SCOREMATRIX2MER);
    if (id == UINT_MAX) {
        return ScoreMatrix();
    }

    // read alphabetsize to remove x (not needed in index)
    PrefilteringIndexData meta = getMetadata(dbr);

    char *data = dbr->getDataUncompressed(id);
    if (preloadMode == Parameters::PRELOAD_MODE_FREAD) {
        return ScoreMatrix::unserializeCopy(data, meta.alphabetSize-1, 2);
    }

    if (preloadMode == Parameters::PRELOAD_MODE_MMAP_TOUCH) {
        dbr->touchData(id);
    }
    return ScoreMatrix::unserialize(data, meta.alphabetSize-1, 2);
}

ScoreMatrix PrefilteringIndexReader::get3MerScoreMatrix(DBReader<unsigned int> *dbr, int preloadMode) {
    size_t id = dbr->getId(SCOREMATRIX3MER);
    if (id == UINT_MAX) {
        return ScoreMatrix();
    }

    // read alphabetsize to remove x (not needed in index)
    PrefilteringIndexData meta = getMetadata(dbr);

    char *data = dbr->getDataUncompressed(id);
    if (preloadMode == Parameters::PRELOAD_MODE_FREAD) {
        return ScoreMatrix::unserializeCopy(data, meta.alphabetSize-1, 3);
    }

    if (preloadMode == Parameters::PRELOAD_MODE_MMAP_TOUCH) {
        dbr->touchData(id);
    }
    return ScoreMatrix::unserialize(data, meta.alphabetSize-1, 3);
}

std::string PrefilteringIndexReader::searchForIndex(const std::string &pathToDB) {
    std::string outIndexName = pathToDB;
    outIndexName.append(".idx");
    if (FileUtil::fileExists(outIndexName.c_str()) == true) {
        return outIndexName;
    }
    return "";
}


