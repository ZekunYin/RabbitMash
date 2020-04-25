/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include <cstdio>
#include <vector>
#include <map>
#include <cstdio>

#include "../Sketch.h"
#include "FastaIO.h"
#include "Buffer.h"
#include "FastaStream.h" 
#include "Sequence.h"

namespace mash
{

namespace fa
{

FastaChunk* FastaReader::readNextChunk(){
	FastaDataChunk* part = NULL;
	recordsPool.Acquire(part);
	FastaChunk *dataPart = new FastaChunk;
	dataPart->chunk = part;
	if(fileReader.ReadNextChunk(dataPart, this->seqInfos))
	{
		return dataPart;
	}
	else
	{
		recordsPool.Release(part);
		return NULL;
	}
}

//int chunkFormat(FastaDataChunk* &chunk, std::vector<Sequence*> &data, bool mHasQuality){
//	//format a whole chunk and return number of reads
//	int seq_count = 0;
//	int line_count = 0;
//	int pos_ = 0;
//
//	while(true){
//		//TODO rewrite to deal with part sequence
//		//string name = getLine(chunk, pos_);
//		//if(name.empty()) break;//dsrc guarantees that read are completed!
//		//std::cout << name << std::endl;
//
//		//string sequence = getSequence(chunk, pos_);
//		//std::cout << sequence << std::endl;
//
//		//data.push_back(new Read(name, sequence));
//		//seq_count++;
//
//	}
//
//	return seq_count;
//}


//FIXME:support enter = \n only
string getSequence(FastaDataChunk * &chunk, uint64 &pos)
{
	int start_pos = pos;
	char * data = (char *)chunk->data.Pointer();
	//cerr << "start pos: " << pos << endl << flush;
	//cerr << "chunk size: " << chunk->size << endl << flush;
	//cerr << "data[pos]: " << (int)data[pos] << endl << flush;
	string res="";

	while(pos < chunk->size - 1 )
	{
		if(data[pos] == '\n'){
			res += string(data+start_pos, pos-start_pos);
			pos++;
			start_pos = pos;
			if(data[pos] == '>')
				return res;		
		}
		else{
			pos++;
		}
	}

	//deal with last char	
	if(pos == chunk->size - 1)
	{
		if(data[pos] == '\n')
			res += string(data+start_pos, pos - start_pos);
		else
			res += string(data+start_pos, pos - start_pos + 1);

		return res;
	}	

		
	return "";
}

//only support uinx-like '\n'
string getLine(FastaDataChunk * &chunk, uint64 &pos)
{
	int start_pos = pos;
	char* data = (char *)chunk->data.Pointer();

	while(pos < chunk->size){
		if(data[pos] == '\n'){
			pos++;
			return string(data+start_pos, pos-start_pos - 1);
		}
		else{
			pos++;
		}
	}

	return "";
}

int chunkFormat(FastaChunk & fachunk, vector<Sketch::Reference> & refs)
{
	uint64 pos = 0;
	bool done = false;
	//cerr << "into chunkFormat" << endl;
	while(true){
	
		Sketch::Reference ref = getNextSeq(fachunk, done, pos);
		if(done) break;
		refs.push_back(ref);
	}

	ASSERT(refs.size() == fachunk.nseqs);

	return refs.size();
	
}

//design to filter out kmer apps
int chunkFormat(FastaChunk & fachunk, vector<Sketch::Reference> & refs, int kmerSize)
{
	uint64 pos = 0;
	bool done = false;
	//cerr << "into chunkFormat" << endl;
	uint64 short_count = 0; //counter for seqs shorter than kmer
	while(true){
	
		Sketch::Reference ref = getNextSeq(fachunk, done, pos);
		if(done) break;
		if(ref.seq.length() < kmerSize)
			short_count++;
		else
			refs.push_back(ref);
	}

	ASSERT(refs.size() + short_count == fachunk.nseqs);

	return refs.size();
	
}

Sketch::Reference getNextSeq(FastaChunk & fachunk, bool & done, uint64 & pos)
{
	Sketch::Reference ref;
	if(pos >= fachunk.chunk->size - 1){
		done = true;
		return ref;
	}

	char *data = (char *)fachunk.chunk->data.Pointer();	

	//while(data[pos] == '\n') pos++;//jump empty lines

	if(data[pos] != '>')
	{
		ref.seq = getSequence(fachunk.chunk, pos);	
		ref.length = ref.seq.size();
		ref.gid = fachunk.start;
		fachunk.start++;
	}else{
		string line = getLine(fachunk.chunk, pos);
		int str_pos = line.find_first_of(' ');
		ref.name = line.substr(1, str_pos-1);	//remove '>' and ' '
		if(str_pos < line.size()) 
			ref.comment = line.substr(str_pos + 1);
		//cerr << "name: " << ref.name << endl << flush;
		ref.seq = getSequence(fachunk.chunk, pos);	
		//cerr << "seq: " << ref.seq << endl << flush;
		ref.length = ref.seq.size();
		ref.gid = fachunk.start;
		fachunk.start++;
	}

	return ref;
}

} // namesapce fa

} // namespace mash
