//output the bin file in text formatting.
//addbyxxm

#include "CommandDumptri.h"
#include <iostream>
#include <fstream>
#include <string>

#include <sys/stat.h> //FIXME: port to windows
#include <sys/time.h>
#include <math.h>

#include "CommandDistance.h"
#include "CommandTriangle.h"
#include "Sketch.h"
#include "CommandDumpdist.h"

#define DIST 1

using namespace::std;

namespace mash{

CommandDumptri::CommandDumptri() : Command()
{
	name = "dumptri";
	summary = "Convert binary triangle results to human-readable texts.";
	description = "Convert binary results produced by \"triangle\" operation to human-readable texts using multiple threads.";
	argumentString = "<seq.msh> [<seq.msh>] <triangle.bin>";

    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <seq> specify paths to sequence files, one per line. The reference file is not affected.", ""));
	//the output format defined by CommandTriangle
    addOption("comment", Option(Option::Boolean, "C", "Output", "Use comment fields for sequence names instead of IDs.", ""));
    //addOption("edge", Option(Option::Boolean, "E", "Output", "Output edge list instead of Phylip matrix, with fields [seq1, seq2, dist, p-val, shared-hashes].", ""));

	addOption("output", Option(Option::String, "o", "Output", "output file", ""));
	useOption("threads");
	
	useOption("help");
	//addOption -o outputfile.
}

int CommandDumptri::run() const
{
	if(arguments.size() < 2 || options.at("help").active)
	{
		print();
		return 0;
	}

	vector<string> queryFiles;

	bool comment = options.at("comment").active;
	//bool edge = options.at("edge").active;
	bool list = options.at("list").active;
	bool edge = false;
	double distanceMax = 1.0;
	double pValueMax = 1.0;

	//string refMsh = arguments[0];
	string fileName = arguments.back();
	string oFileName = options.at("output").argument;

	for(int i = 0; i < arguments.size() - 1; i++){
		if(list)
		{
			splitFile(arguments[i], queryFiles);
		}
		else
		{
			queryFiles.push_back(arguments[i]);
		}
	
	}
	for(int i = 0; i < queryFiles.size(); i++)
	{
		if(!hasSuffix(queryFiles[i], ".msh")){
			cerr << queryFiles[i] << " is not msh format, please provide correct input" << endl;
			exit(1);
		}

	}

	if(oFileName == "")
		oFileName = fileName + ".tri";
	FILE * fout = fopen(oFileName.c_str(), "wb");
	if(fout == NULL){
		cerr << "can not open result file: " << oFileName << endl;
		exit(1);
	}
	else{
		cerr << "writting result to: " << oFileName << endl;
	}
	
	fstream resultFile(fileName.c_str(), ios::binary | ios::in);
	if(!resultFile.is_open()){
		cerr << "fail to open " << fileName << endl;
		return 1;
	}

	Sketch::Parameters parameters;
	//TODO: setup parameters
	int threads = options.at("threads").getArgumentAsNumber();

	Sketch querySketch;
	querySketch.initFromFiles(queryFiles, parameters);

	int64_t binSize = getFileSize(fileName.c_str());
	binSize -= sizeof(uint64_t) + 2 * sizeof(double);//edge header
	int64_t resSize = binSize / sizeof(CommandTriangle::Result);
	if(binSize % sizeof(CommandTriangle::Result) != 0)
	{
		cerr << "imcomplete binary file" << endl;
		exit(1);
	}
	
	cerr << "query sketches: " << querySketch.getReferenceCount() << endl;
	cerr << "binary file size: " << binSize << endl;
	cerr << "number of results: " << resSize << endl;

	//read header edge distance pvalue
	uint64_t header = 0;
	resultFile.read((char*)&header, sizeof(uint64_t));	
	resultFile.read((char*)&distanceMax, sizeof(double));	
	resultFile.read((char*)&pValueMax  , sizeof(double));	
	if(header == 0x1) edge = true;

	if(edge)
	{
		cerr << "distanceMax: " << distanceMax << endl;
		cerr << "pValueMax: "   << pValueMax   << endl;

	} else {
	
		int64_t correctResSize = querySketch.getReferenceCount() * (querySketch.getReferenceCount()-1) / 2;
		//cerr << "correctResSize of results: " << correctResSize << endl;

		if(resSize != correctResSize)
		{
			cerr << "unmatched msh file or bin file"  << endl;
			cerr << "please checkout whether the msh file and bin file is from the same data and parameters" << endl;
			
			exit(1);
		}
		//cerr << "msh and binary file are checked out !" << endl;

	}

	CommandTriangle::Result * buffer = new CommandTriangle::Result [binSize / sizeof(CommandTriangle::Result)];
	resultFile.read((char*)buffer, binSize);
	
	string oFilePrefix = fileName + ".tri";
	if(edge){
		#pragma omp parallel for default(shared) num_threads(threads)
		for(int i = 0; i < threads; i++)
		{
			fstream oFile(oFilePrefix + to_string(i), ios::out | ios::binary | ios:: trunc);
			int64_t start = i * ((resSize + threads -1) /threads);
			int64_t end = resSize < (i + 1) * ((resSize + threads - 1) / threads)
						? resSize : (i + 1) * ((resSize + threads - 1) / threads);
			
			for(int j = start; j < end; j++){
				if(!comment){
					string tmp = 
							querySketch.getReference(buffer[j].refID).name + "\t"
							+ querySketch.getReference(buffer[j].queryID).name + "\t"
							+ to_string(buffer[j].distance) + "\t"
							+ to_string(buffer[j].pValue) + "\t"
							+ to_string(buffer[j].numer) + "/"
							+ to_string(buffer[j].denom) + "\n";
					oFile.write(tmp.c_str(), tmp.size());
				}
				else{
					string tmp = 
							querySketch.getReference(buffer[j].refID).comment + "\t"
							+ querySketch.getReference(buffer[j].queryID).comment + "\t"
							+ to_string(buffer[j].distance) + "\t"
							+ to_string(buffer[j].pValue) + "\t"
							+ to_string(buffer[j].numer) + "/"
							+ to_string(buffer[j].denom) + "\n";
					oFile.write(tmp.c_str(), tmp.size());
				}
			}
		}
	}//end if(edge)

	else{//not edge triangle output
		
		//add the first two line of the triangle.
		fstream oFileHead(oFilePrefix + "head", ios::out | ios::binary | ios::trunc);
		string tmp1 = to_string(querySketch.getReferenceCount()) + "\n";
		oFileHead.write(tmp1.c_str(), tmp1.size());
		string tmp2 = (comment ? querySketch.getReference(0).comment : querySketch.getReference(0).name) + "\n";
		oFileHead.write(tmp2.c_str(), tmp2.size());
		oFileHead.close();

		int mean = (resSize + threads - 1) / threads;
		#pragma omp parallel for default(shared) num_threads(threads)
		for(int i = 0; i < threads; i++)
		{
			fstream oFile(oFilePrefix + to_string(i), ios::out | ios::binary | ios::trunc);
			int startLine = sqrt(i * mean * 2);
			int endLine = sqrt((i+1) * mean * 2);
			int start = (startLine+1)*(startLine+2)/2 - 1;
			int end = resSize < (endLine+1)*(endLine+2)/2 ? resSize-1 : (endLine+1)*(endLine+2)/2 - 1;
			
			//int curLine = startLine+1+ ceil((double)i/threads);
			int curLine = startLine + 1;
			string tmp;

			if(!comment)
				tmp = querySketch.getReference(buffer[start].refID).name;
			else
				tmp = querySketch.getReference(buffer[start].refID).comment;

			for(int j = start; j <= end; j++){
				if(!comment){
					tmp += "\t";
					tmp += to_string(buffer[j].distance);
					if(j+1 == curLine * (curLine+1) / 2){
						if(start == 0 || j != start){//erase multiThread redundant output
							tmp +="\n";
							oFile.write(tmp.c_str(), tmp.size());
						}
						curLine++;
						tmp = 
						querySketch.getReference(buffer[j+1].refID).name;
					}

				}
				else{
					tmp += "\t";
					tmp += to_string(buffer[j].distance);
					if(j+1 == curLine * (curLine+1) / 2){
						if(start == 0 || j != start){
							tmp +="\n";
							oFile.write(tmp.c_str(), tmp.size());
						}
						curLine++;
						tmp = 
						querySketch.getReference(buffer[j+1].refID).comment;
					}

				}
			}

//			tmp +="\n";
//			oFile.write(tmp.c_str(), tmp.size());

			oFile.close();

		}

	}//end else(not edge)
	
	double t1 = get_sec();
	int tmpSize = 1<<20;
	unsigned char *tmpBuffer = new unsigned char[tmpSize];

	if(!edge){
		FILE *tmpFileHead = fopen((oFilePrefix + "head").c_str(), "rb");
		if(tmpFileHead == NULL){
			cerr << "can not open" << (oFilePrefix + "head") << endl;
			exit(1);
		}
		int lengthHead;
		while(true){
			lengthHead = fread(tmpBuffer, sizeof(unsigned char), tmpSize, tmpFileHead);
			fwrite((void*)tmpBuffer, sizeof(unsigned char), lengthHead, fout);
			if(lengthHead < tmpSize) break;
		}
		remove((oFilePrefix + "head").c_str());
	}

	for(int i = 0; i < threads; i++){
		FILE *tmpFile = fopen((oFilePrefix + to_string(i)).c_str(), "rb");
		if(tmpFile == NULL){
			cerr << "can not open " << (oFilePrefix + to_string(i)) << endl;
			exit(1);
		}

		int n;
		while(true){
			n = fread(tmpBuffer, sizeof(unsigned char), tmpSize, tmpFile);
			fwrite((void*)tmpBuffer, sizeof(unsigned char), n, fout);
			if(n < tmpSize) break;
		}

		fclose(tmpFile);
		remove((oFilePrefix + to_string(i)).c_str());
	
	}

	double t2 = get_sec();

	cerr << "combine time is: " << t2 - t1 << endl;

	resultFile.close();
	fclose(fout);
	
	delete tmpBuffer;
	delete buffer;

	return 0;
}

}//namespace mash

