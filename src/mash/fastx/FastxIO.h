/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTXREADER
#define H_FASTXREADER

#include <vector>
#include <string>

#include "Globals.h"
#include "Common.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "FastxStream.h"
#include "FastxChunk.h"
//#include "Sequence.h"
#include "../Sketch.h"

namespace mash
{

namespace fa
{


class FastaReader
{
public:
    FastaReader(FastaFileReader& reader_, FastaDataPool& pool_)
	    :   recordsPool(pool_)
		,	fileReader(reader_)
		,	numParts(0)
	{};

	FastaChunk* readNextChunk();
	
	int64 Read(byte* memory_, uint64 size_)
	{
		int64 n = fileReader.Read(memory_, size_);
		return n;
	}

public:	
	SeqInfos seqInfos;

private:

	FastaDataPool&      recordsPool;
	FastaFileReader&	fileReader;
	uint32 numParts;
};

//int chunkFormat(FastaDataChunk* &chunk, std::vector<Sequence*>&, bool);


std::string getSequence(FastaDataChunk* &chunk, uint64 &pos);	//addbyxxm
std::string getLine(FastaDataChunk* &chunk, uint64 &pos);
int chunkFormat(FastaChunk & fachunk, std::vector<Sketch::Reference> & refs);
int chunkFormat(FastaChunk & fachunk, std::vector<Sketch::Reference> & refs, int kmerSize);
Sketch::Reference getNextSeq(FastaChunk & fachunk, bool & done, uint64 & pos);

} // namespace fa

namespace fq
{

//moved to FastxChunk.h
//typedef core::TDataQueue<FastqDataChunk> FastqDataQueue;
//typedef core::TDataPool<FastqDataChunk> FastqDataPool;


class FastqReader //: public IFastqIoOperator
{
public:
    FastqReader(FastqFileReader& reader_, FastqDataPool& pool_)
	    :   recordsPool(pool_)
		,	fileReader(reader_)
		,	numParts(0)
	{};

	void readChunk();
	FastqDataChunk* readNextChunk();
	//single pe file
	FastqDataChunk* readNextPairedChunk();
	
	int64 Read(byte* memory_, uint64 size_)
	{
		int64 n = fileReader.Read(memory_, size_);
		return n;
	}

private:

	FastqDataPool&      recordsPool;
	FastqFileReader&	fileReader;
	uint32 numParts;
};

int chunkFormat(FastqChunk* &chunk, std::vector< Sketch::Reference > &,bool);

//single pe file 
//int pairedChunkFormat(FastqDataChunk* &chunk, std::vector<ReadPair*>&,bool mHasQuality);
//Read* getOnePairedRead(FastqDataChunk* &chunk,int &pos_, bool mHasQuality);
//end single pe file
std::string getLine(FastqDataChunk* &chunk, int &pos);

} // namespace fq

} // namespace mash

#endif
