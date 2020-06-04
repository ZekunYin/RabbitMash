/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00

  last modified by Zekun Yin 2020/5/18
*/
 
#include "FastxStream.h"
#include "../Sketch.h"
#include <iostream>
#include <string>

namespace mash
{

namespace fa
{

bool FastaFileReader::ReadNextChunk(FastaChunk* dataChunk_, SeqInfos& seqInfos)
{
	//std::cout << "==================Next Chunk =========================" << std::endl;
	FastaDataChunk *chunk_ = dataChunk_->chunk;
	if (Eof())
	{
		chunk_->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	uchar* data = chunk_->data.Pointer();
	const uint64 cbufSize = chunk_->data.Size();
	chunk_->size = 0;
	int64 toRead = cbufSize - bufferSize;// buffersize: size left from last chunk
	
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = this->Read(data + chunk_->size, toRead);
	//std::cout << "r is :" << r << std::endl;
	//std::cout << "toRead: " << toRead << std::endl;
	
	if (r > 0)
	{
		if (r == toRead)	// somewhere before end
		{
		    //uint64 chunkEnd = cbufSize - SwapBufferSize; // Swapbuffersize: 1 << 13
			////std::cout << "chunkend  cbufsize Swapbuffersize: " << chunkEnd <<" "<< cbufSize << " " << SwapBufferSize << std::endl;
			//chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
			//chunk_->size = chunkEnd - 1; //remove \n
			//if (usesCrlf)
			//	chunk_->size -= 1;

			//std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
			//bufferSize = cbufSize - chunkEnd;

			//dealing with halo region
			uint64 chunkEnd = FindCutPos(dataChunk_, data, cbufSize, mHalo, seqInfos);
			chunk_->size = chunkEnd;// - 1; //1 char back from last '>'

			//debug only
			//std::string content((char*)data, chunk_->size);
			//std::cout << "chunkEnd data: " << (char)data[chunkEnd-1] << std::endl;
			//std::cout << "chunkEnd: " << (uint64)chunkEnd << std::endl;
			//std::cout << "chunk_->size: " << chunk_->size << std::endl;
			//std::cout << "content_ori: " << data << std::endl;
			//std::cout << "content    : " << content << std::endl;
			//end debug
			
			if (usesCrlf)
				chunk_->size -= 1;
			//copy tail to swapBuffer
			//if(data[chunkEnd] == '\n') chunkEnd++;
			//TODO: dealing with halo region
			//get rid of '\n' in halo region FIXME: Crlf not supported
			int haloCount = 0;
			uint64 tmpPtr = chunkEnd - 1;
			if(mHalo > 0){
				while(true){
					if(data[tmpPtr] != '\n')
					{
						haloCount++;
						//std::cerr << (char)data[tmpPtr] << std::endl;
						if(haloCount == mHalo) break;
					}
					tmpPtr--;
				}
				//std::cerr << "haloCount: " << haloCount << std::endl;
			}else{
				tmpPtr = chunkEnd; //nocopy	
			}
			//std::copy(data + chunkEnd - mHalo, data + cbufSize, swapBuffer.Pointer());
			std::copy(data + tmpPtr, data + cbufSize, swapBuffer.Pointer());
			//bufferSize = cbufSize - chunkEnd + mHalo;
			bufferSize = cbufSize - tmpPtr;

		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
			if (usesCrlf)
				chunk_->size -= 1;

			//only for get seqsinfo
			uint64 chunkEnd = FindCutPos(dataChunk_, data, chunk_->size, mHalo, seqInfos);
			//debug only
			//std::string content((char*)data, chunk_->size);
			//std::cout << "tail content ori: " << data << std::endl;
			//std::cout << "tail content    : " << content << std::endl;
			//end debug

			eof = true;
		}
	}
	else
	{
		eof = true;
	}

	return true;
}

uint64 FastaFileReader::FindCutPos(FastaChunk* dataChunk_, uchar* data_, const uint64 size_, const uint64 halo_, SeqInfos& seqInfos)
{
	int count = 0;
	uint64 pos_ = 0;
	uint64 cut_ = 0; //cut_ point to next '>'
	uint64 lastSeq_ = 0; //-> the start of last sequences content
	uint64 lastName_ = 0; //-> the last '>'
	OneSeqInfo seqInfo;

	if(data_[0] == '>') //start with '>'
	{
		dataChunk_->start = this->totalSeqs;	
		while(pos_ < size_){
			if(data_[pos_] == '>'){
				lastName_ = pos_;
				if(FindEol(data_, pos_, size_)) //find name
				{
					++pos_;
					lastSeq_ = pos_;

					seqInfo.gid = this->totalSeqs;
					seqInfos.push_back(seqInfo);

					this->totalSeqs++;
				}else{
					cut_ = pos_; //find a cut: incomplete name
					//std::cerr << "cut char: " << (char)data_[cut_] << std::endl;
					break;
				}
			}else{
				++pos_;
			}
	
		}

	}else{ //start with {ACGT}
		dataChunk_->start = this->totalSeqs - 1;	
		while(pos_ < size_){
			if(data_[pos_] == '>'){
				lastName_ = pos_;
				if(FindEol(data_, pos_, size_)) //find name
				{
					++pos_;
					lastSeq_ = pos_;
					seqInfo.gid = this->totalSeqs;
					seqInfos.push_back(seqInfo);
					this->totalSeqs++;
				}else{
					cut_ = pos_; //find a cut -> '>'
					//std::cout << "cut char: " << (char)data_[cut_] << std::endl;
					break;
				}
			}else{
				++pos_;
			}
	
		}	
	}

	//no tail cut
	//make sure that cbufSize > name_len + halo
	if(cut_ == 0){
		uint64 lastSeqLen_ = size_ - lastSeq_;	
		if(lastSeqLen_ < halo_){
			cut_ = lastName_;	
			this->totalSeqs--;
		}
	}
	
	dataChunk_->nseqs = this->totalSeqs - dataChunk_->start;
	dataChunk_->end = this->totalSeqs - 1;

	//if(cut_ != 0) std::cout << "cut: " << cut_ << std::endl;

	return cut_ ? cut_ : size_;	
}

} // namespace fa

namespace fq
{

//bool IFastqStreamReader::ReadNextChunk(FastqDataChunk* chunk_)
bool FastqFileReader::ReadNextChunk(FastqDataChunk* chunk_)
{
	if (Eof())
	{
		chunk_->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	uchar* data = chunk_->data.Pointer();
	const uint64 cbufSize = chunk_->data.Size();
	chunk_->size = 0;
	int64 toRead = cbufSize - bufferSize;
	//---------------------------------
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = Read(data + chunk_->size, toRead);
	//std::cout << "r is :" << r << std::endl;
	
	if (r > 0)
	{
		if (r == toRead)	// somewhere before end
		{
		    uint64 chunkEnd = cbufSize - SwapBufferSize; // Swapbuffersize: 1 << 13
			//std::cout << "chunkend  cbufsize Swapbuffersize: " << chunkEnd <<" "<< cbufSize << " " << SwapBufferSize << std::endl;
			chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
			chunk_->size = chunkEnd - 1;
			if (usesCrlf)
				chunk_->size -= 1;

			std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
			bufferSize = cbufSize - chunkEnd;
		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
			if (usesCrlf)
				chunk_->size -= 1;

			eof = true;
		}
	}
	else
	{
		eof = true;
	}

	return true;
}
//bool FastqFileReader::ReadNextPairedChunk(FastqDataChunk* chunk_)
//{	//std::cout<<"ReadNextPairedChunk: enter" <<std::endl;
//	if (Eof())
//	{
//		chunk_->size = 0;
//		return false;
//	}
//
//	// flush the data from previous incomplete chunk
//	uchar* data = chunk_->data.Pointer();
//	const uint64 cbufSize = chunk_->data.Size();
//	chunk_->size = 0;
//	int64 toRead = cbufSize - bufferSize;
//	//---------------------------------
//	if (bufferSize > 0)
//	{
//		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
//		chunk_->size = bufferSize;
//		bufferSize = 0;
//	}
//
//	// read the next chunk
//	int64 r = Read(data + chunk_->size, toRead);
//	//std::cout << "r is :" << r << std::endl;
//	
//	if (r > 0)
//	{
//		if (r == toRead)	// somewhere before end
//		{
//		    //uint64 chunkEnd = cbufSize - SwapBufferSize; // Swapbuffersize: 1 << 13
//			//std::cout << "chunkend  : " << cbufSize-1 <<std::endl;
//			lastOneReadPos = GetPreviousRecordPos(data, cbufSize-1,cbufSize);
//			//std::cout << "one：" << lastOneReadPos<<std::endl;
//			lastTwoReadPos = GetPreviousRecordPos(data, lastOneReadPos-1,cbufSize);
//			//std::cout << "two :" << lastTwoReadPos<<std::endl;
//			chunk_->size = lastOneReadPos - 1;
//			if (usesCrlf) // TODO windows先不管
//				chunk_->size -= 1;
//
//			std::copy(data + lastTwoReadPos, data + cbufSize, swapBuffer.Pointer());
//			bufferSize = cbufSize - lastTwoReadPos;
//		}
//		else				// at the end of file
//		{
//			chunk_->size += r - 1;	// skip the last EOF symbol
//			if (usesCrlf)
//				chunk_->size -= 1;
//
//			eof = true;
//		}
//	}
//	else
//	{
//		eof = true;
//	}
//	//std::cout<<"ReadNextPairedChunk: success@!!!!!!!" <<std::endl;
//	return true;
//}
//
//需要向前找到倒数第2条read的开头，不仅仅是@号这么简单
uint64 FastqFileReader::GetPreviousRecordPos(uchar* data_, uint64 pos_,const uint64 size_)
{	int offset =2;

	SkipToSol(data_,pos_,size_);
	if(usesCrlf){
		offset=3;
	}
	while(data_[pos_+offset] !='@'){ //+2
		//std::cout<<"pos_"<<pos_<<std::endl;
		SkipToSol(data_,pos_,size_);
	}
	//标记一下，看是否是质量分
	uint64 pos0=pos_+offset;
	SkipToSol(data_,pos_,size_);
	if(data_[pos_+offset]== '+'){
		//说明上一个@号是质量分
		SkipToSol(data_,pos_,size_);
		SkipToSol(data_,pos_,size_);
		//此时应该是name
		if(data_[pos_+offset]!='@'){
			std::cout << "core dump is " << data_[pos_+offset] << std::endl;
			return pos_+offset;
		}else{
			return pos_+offset;
		}
	}
	else{
		return pos0;
	}
}




//uint64 IFastqStreamReader::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
uint64 FastqFileReader::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
{
	SkipToEol(data_, pos_, size_);
	++pos_;

	// find beginning of the next record
	while (data_[pos_] != '@')
	{
		SkipToEol(data_, pos_, size_);
		++pos_;
	}
	uint64 pos0 = pos_;
		
	SkipToEol(data_, pos_, size_);
	++pos_;

	if (data_[pos_] == '@')			// previous one was a quality field
		return pos_;
	
	SkipToEol(data_, pos_, size_);
	++pos_;
	if(data_[pos_] != '+')
		std::cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << std::endl;	
	ASSERT(data_[pos_] == '+');	// pos0 was the start of tag
	return pos0;
}

} // namespace fq

} // namespace mash

