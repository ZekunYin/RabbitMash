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

#define DIST 1

using namespace::std;
int64_t getFileSize_( const char * fileName)
{
	struct stat statbuf;
	stat(fileName, &statbuf);
	return (int64_t)statbuf.st_size;
}

double get_sec_(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

namespace mash{

CommandDumpTri::CommandDumpTri() : Command()
{
<<<<<<< HEAD
	name = "dumpTri";
	summary = "output the bin file in text formatting.";
	description = "in the last step of the \" triangle \" we get two binary out file name.bin and dist.bin but both of them are unreadable. The outputbin is used to output the bin file in text formatting.";
	argumentString = "<ref.msh> <triangle.bin>";

	//the output format defined by CommandTriangle
    addOption("comment", Option(Option::Boolean, "C", "Output", "Use comment fields for sequence names instead of IDs.", ""));
    addOption("edge", Option(Option::Boolean, "E", "Output", "Output edge list instead of Phylip matrix, with fields [seq1, seq2, dist, p-val, shared-hashes].", ""));
    //addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report in edge list. Implies -" + getOption("edge").identifier + ".", "1.0", 0., 1.));
    //addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report in edge list. Implies -" + getOption("edge").identifier + ".", "1.0", 0., 1.));

	addOption("output", Option(Option::String, "o", "Output", "output file", ""));
	useOption("threads");
=======
	name = "dumptri";
	summary = "Convert binary triangle results to human-readable texts.";
	description = "Convert binary results produced by \"triangle\" operation to human-readable texts using multiple threads.";
	argumentString = "<name.bin> <dist.bin>";
>>>>>>> 00a28fdc22fd8aaa2b48fd6200b837ad69d242a5
	
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
	int threads = options.at("threads").getArgumentAsNumber();

	Sketch refSketch;
	vector<string> refFiles;
	refFiles.push_back(refMsh);
	refSketch.initFromFiles(refFiles, parameters);

	int64_t binSize = getFileSize_(fileName.c_str());
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
	
	cerr << "end the resultFile.read " << endl;
	cerr << "the threads number is " << threads << endl;
	string oFilePrefix = fileName + ".tri";
	if(edge){
		#pragma omp parallel for default(shared) num_threads(threads)
		for(int i = 0; i < threads; i++)
		{
			fstream oFile(oFilePrefix + to_string(i), ios::out | ios::binary | ios:: trunc);
			int64_t start = i * ((resSize + threads -1) /threads);
			int64_t end = resSize < (i + 1) * ((resSize + threads - 1) / threads)
						? resSize : (i + 1) * ((resSize + threads - 1) / threads);
			
		cerr << "the start is: " << start << endl;
		cerr << "the end is: " << end << endl;
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
	}

	else{//not edge triangle output
		int mean = (resSize + threads - 1) / threads;
		#pragma omp parallel for default(shared) num_threads(threads)
		for(int i = 0; i < threads; i++)
		{
			fstream oFile(oFilePrefix + to_string(i), ios::out | ios::binary | ios::trunc);
			int startLine = sqrt(i * mean * 2);
			int endLine = sqrt((i+1) * mean * 2);
			int start = (startLine+1) * (startLine+2) / 2;
			int end = resSize < (endLine+1) * (endLine+2) / 2 ? resSize : (endLine+1) * (endLine+2) / 2;
			
			int curLine = startLine+1;
			string tmp;
			cout << "curLine is: " << curLine << endl;

			if(!comment)
				//cout << refSketch.getReference(curLine-1).name;
				//tmp = refSketch.getReference(buffer[start].refID).name;
				tmp = refSketch.getReference(curLine).name;
			else
				//cout << refSketch.getReference(curLine-1).comment;
				//tmp = refSketch.getReference(buffer[start].refID).comment;
				tmp = refSketch.getReference(curLine).comment;
			for(int j = start; j < end; j++){
				if(!comment){
					//cout << "\t" << to_string(buffer[j].distance);
					tmp += "\t";
					tmp += to_string(buffer[j].distance);
					if(j % (curLine * (curLine+1) / 2) == 0){
					//	cout << endl;
					//	curLine++;
					//	cout << refSketch.getReference(buffer[j].refID).name;
						tmp +="\n";
						oFile.write(tmp.c_str(), tmp.size());
						curLine++;
						tmp = 
						refSketch.getReference(curLine).name;
						//refSketch.getReference(buffer[j].refID+1).name;
					}

				}
				else{
					tmp += "\t";
					tmp += to_string(buffer[j].distance);
					if(j % (curLine * (curLine+1) / 2) == 0){
						tmp += "\n";
						oFile.write(tmp.c_str(), tmp.size());
						tmp = 
						refSketch.getReference(curLine).comment;
						//refSketch.getReference(buffer[j].refID).comment;
					}

				}
			}
			oFile.close();

		}
		

	}
	
	double t1 = get_sec_();
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

	double t2 = get_sec_();

	cerr << "combine time is: " << t2 - t1 << endl;

	resultFile.close();
	fclose(fout);
	
	delete tmpBuffer;
	delete buffer;

	return 0;
}

}//namespace mash


			



//	string nameFile = arguments[0];
//	string distFile = arguments[1];
//
//	fstream out1(nameFile.c_str(), out1.binary | out1.in);
//	fstream out2(distFile.c_str(), out2.binary | out2.in);
//
//	if(!out1.is_open()){
//		cerr << "fail to open " << nameFile << endl;
//		return 1;
//	}
//	if(!out2.is_open()){
//		cerr << "fail to open " << distFile << endl;
//		return 1;
//	}
//
//	double d;
//	char * name = new char[20];
//
//	int lineSize = 1;
//	int lineIndex = 0;
//	int readNum = 0;
//
//	out1.read(reinterpret_cast<char*>(name), 20 * sizeof(char));
//	for(int i = 0; i < 20; i++){
//		if(name[i] < 32) break;//invalid ASCII
//		cout << name[i];
//	}
//
//	while(1){
//		out2.read(reinterpret_cast<char*>(&d), sizeof d);
//		if(out2.eof()) break;
//		
//		readNum++;
//		if(d < 1.0){
//			if(lineIndex < lineSize){
//				cout << '\t' << d;
//				lineIndex++;
//			}
//			else{
//				cout << endl;
//				out1.read(reinterpret_cast<char*>(name), 20 * sizeof(char));
//				for(int i = 0; i < 20; i++){
//					if(name[i] < 32) break;//invalid ASCII
//					cout << name[i];
//				}
//				cout << '\t' << d;
//				lineIndex = 1;
//				lineSize++;
//			}
//		}
//
//		else{//d >= 1
//			for(int j = 0; j < d / 1; j++){
//				if(lineIndex < lineSize){
//					//cout << DIST << '\t';
//					cout << '\t' << DIST;
//					lineIndex++;
//				}
//				else{
//					cout << endl;
//					out1.read(reinterpret_cast<char*>(name), 20 * sizeof(char));
//					for(int i = 0; i < 20; i++){
//						if(name[i] < 32) break;//invalid ASCII
//						cout << name[i];
//					}
//					//cout << DIST << '\t';
//					cout << '\t' << DIST;
//
//					lineIndex = 1;
//					lineSize++;
//				}
//			}
//		} 
//			
//	}
//	cout << endl;//make sure the tail line endl in order to get the same md5-value with Mash.
//	cerr << "the readNum is: " << readNum << endl;
//	
//
//	delete name;
//	out1.close();
//	out2.close();
//
//	return 0;
//}
//
//} //namespace mash

