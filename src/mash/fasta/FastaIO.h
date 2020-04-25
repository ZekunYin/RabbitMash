/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQREADER
#define H_FASTQREADER

#include <vector>
#include <string>

#include "Globals.h"
#include "Common.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "FastaStream.h"
#include "FastaChunk.h"
#include "Sequence.h"
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

int chunkFormat(FastaDataChunk* &chunk, std::vector<Sequence*>&, bool);


string getSequence(FastaDataChunk* &chunk, uint64 &pos);	//addbyxxm
string getLine(FastaDataChunk* &chunk, uint64 &pos);
int chunkFormat(FastaChunk & fachunk, vector<Sketch::Reference> & refs);
int chunkFormat(FastaChunk & fachunk, vector<Sketch::Reference> & refs, int kmerSize);
Sketch::Reference getNextSeq(FastaChunk & fachunk, bool & done, uint64 & pos);

} // namespace fa

} // namespace mash

#endif
