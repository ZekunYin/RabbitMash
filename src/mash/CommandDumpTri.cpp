//output the bin file in text formatting.
//addbyxxm

#include "CommandDumpTri.h"
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
//int64_t getFileSize_( const char * fileName)
//{
//	struct stat statbuf;
//	stat(fileName, &statbuf);
//	return (int64_t)statbuf.st_size;
//}
//
//double get_sec_(){
//	struct timeval tv;
//	gettimeofday(&tv, NULL);
//	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
//}

namespace mash{

CommandDumpTri::CommandDumpTri() : Command()
{
	name = "dumptri";
	summary = "Convert binary triangle results to human-readable texts.";
	description = "Convert binary results produced by \"triangle\" operation to human-readable texts using multiple threads.";
	argumentString = "<name.bin> <dist.bin>";
	//the output format defined by CommandTriangle
    addOption("comment", Option(Option::Boolean, "C", "Output", "Use comment fields for sequence names instead of IDs.", ""));
    addOption("edge", Option(Option::Boolean, "E", "Output", "Output edge list instead of Phylip matrix, with fields [seq1, seq2, dist, p-val, shared-hashes].", ""));

	addOption("output", Option(Option::String, "o", "Output", "output file", ""));
	useOption("threads");
	
	useOption("help");
	//addOption -o outputfile.
}

int CommandDumpTri::run() const
{
	if(arguments.size() < 2 || options.at("help").active)
	{
		print();
		return 0;
	}

	bool comment = options.at("comment").active;
	bool edge = options.at("edge").active;

	string refMsh = arguments[0];
	string fileName = arguments[1];
	string oFileName = options.at("output").argument;

	if(!hasSuffix(refMsh, ".msh")){
		cerr << refMsh << "is not msh format, please provide correct input" << endl;
		exit(1);
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

	Sketch refSketch;
	vector<string> refFiles;
	refFiles.push_back(refMsh);
	refSketch.initFromFiles(refFiles, parameters);

	int64_t binSize = getFileSize(fileName.c_str());
	int64_t resSize = binSize / sizeof(CommandTriangle::Result);
	if(binSize % sizeof(CommandTriangle::Result) != 0)
	{
		cerr << "imcomplete binary file" << endl;
		exit(1);
	}
	cerr << "ref sketches: " << refSketch.getReferenceCount() << endl;
	cerr << "binary file size: " << binSize << endl;
	cerr << "number of results: " << resSize << endl;

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
				if(edge)
				{
					if(!comment){
					//if(pair->pass)
						string tmp = 
								refSketch.getReference(buffer[j].refID).name + "\t"
								+ refSketch.getReference(buffer[j].queryID).name + "\t"
								+ to_string(buffer[j].distance) + "\t"
								+ to_string(buffer[j].pValue) + "\t"
								+ to_string(buffer[j].numer) + "/"
								+ to_string(buffer[j].denom) + "\n";
						oFile.write(tmp.c_str(), tmp.size());
					}
					else{
						string tmp = 
								refSketch.getReference(buffer[j].refID).comment + "\t"
								+ refSketch.getReference(buffer[j].queryID).comment + "\t"
								+ to_string(buffer[j].distance) + "\t"
								+ to_string(buffer[j].pValue) + "\t"
								+ to_string(buffer[j].numer) + "/"
								+ to_string(buffer[j].denom) + "\n";
						oFile.write(tmp.c_str(), tmp.size());
					}
				}
				else
				{
	
				}
			}
		}
	}//end if(edge)

	else{//not edge triangle output
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
				tmp = refSketch.getReference(buffer[start].refID).name;
			else
				tmp = refSketch.getReference(buffer[start].refID).comment;

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
						refSketch.getReference(buffer[j+1].refID).name;
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
						refSketch.getReference(buffer[j+1].refID).comment;
					}

				}
			}

//			tmp +="\n";
//			oFile.write(tmp.c_str(), tmp.size());

			oFile.close();

		}
		

	}
	
	double t1 = get_sec();
	int tmpSize = 1<<20;
	unsigned char *tmpBuffer = new unsigned char[tmpSize];

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

