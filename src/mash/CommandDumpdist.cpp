//covert binary dist output to text format
//yzk 2020/5/28

#include "CommandDumpdist.h"
#include "CommandDistance.h"
#include "Sketch.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h> //FIXME: port to windows
#include <sys/time.h>

#include <stdint.h>
#include <cstdio>

using namespace::std;


namespace mash{

CommandDumpdist::CommandDumpdist() : Command()
{
	name = "dumpdist";
	summary = "Convert binary dist results to human-readable texts.";
	description = "Convert binary results produced by \"dist\" operation to human-readable texts using multiple threads.";
	argumentString = "<reference.msh> <query.msh> [<query.msh>] <dist.bin>";
	
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <query> specify paths to sequence files, one per line. The reference file is not affected.", ""));
    addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
    addOption("comment", Option(Option::Boolean, "C", "Output", "Show comment fields with reference/query names (denoted with ':').", "1.0", 0., 1.));
	addOption("output", Option(Option::String, "o", "Output", "output human-readable text file", ""));

	useOption("threads");
	useOption("help");
}

int64_t getFileSize( const char * fileName)
{
	struct stat statbuf;
	stat(fileName, &statbuf);
	return (int64_t)statbuf.st_size;
}


int CommandDumpdist::run() const
{
	if(arguments.size() < 3 || options.at("help").active)
	{
		print();
		return 0;
	}

	vector<string> queryFiles;
	vector<string> refFiles;

	string refMsh   = arguments.front();
	refFiles.push_back(refMsh);
	string fileName = arguments.back();

	string oFileName = options.at("output").argument;

    bool table = options.at("table").active;
    bool comment = options.at("comment").active;
    bool list = options.at("list").active;

	if(!hasSuffix(refMsh, ".msh")){
		cerr << refMsh << " is not msh format, please provide correct input" << endl;
		exit(1);
	}

	for(int i = 1; i < arguments.size() - 1; i++)
	{           
		if ( list )
            {
                splitFile(arguments[i], queryFiles);
            }
            else
            {
                queryFiles.push_back(arguments[i]);
            }

	}

	for( int i = 0; i < queryFiles.size(); i++ )
	{
		if(!hasSuffix(queryFiles[i], ".msh")){
			cerr << queryFiles[i] << " is not msh format, please provide correct input" << endl;
			exit(1);
		}
	}

	if(oFileName == "") oFileName = fileName + ".dist";
	FILE * fout = fopen(oFileName.c_str(), "wb");

	if(fout == NULL){
		cerr << "can not open result file: " << oFileName << endl;
		exit(1);
	}else{
		cerr << "writting result to " << oFileName << endl;
	}

	fstream resultFile(fileName.c_str(), ios::binary | ios::in);
	if(!resultFile.is_open()){
		cerr << "fail to open " << fileName << endl;
		return 1;
	}
	
	Sketch::Parameters parameters;
	//TODO:	setup parameters 
	int threads = options.at("threads").getArgumentAsNumber();

	Sketch querySketch;
	Sketch refSketch;


	//queryFiles.push_back(queryMsh);
	
	querySketch.initFromFiles(queryFiles, parameters);
	refSketch.initFromFiles(refFiles, parameters);

	int64_t binSize = getFileSize(fileName.c_str());
	int64_t resSize = binSize / sizeof(CommandDistance::Result);
	if(binSize % sizeof(CommandDistance::Result) != 0)
	{
		cerr << "imcomplete binary file" << endl;
		exit(1);
	}

	cerr << "ref sketches: "   << refSketch.getReferenceCount()   << endl;
	cerr << "query sketches: " << querySketch.getReferenceCount() << endl;
	cerr << "binary file size: " << binSize << endl;
	cerr << "number of results: " << resSize << endl;

	CommandDistance::Result *buffer = new CommandDistance::Result[binSize / sizeof(CommandDistance::Result)];

	resultFile.read((char*)buffer, binSize); //FIXME: portable but not memory efficient
	
	string oFilePrefix = fileName + ".dist";
#pragma omp parallel for default(shared) num_threads(threads)
	for(int i = 0; i < threads; i++)
	{
		fstream oFile(oFilePrefix + to_string(i), ios::out | ios::binary | ios::trunc);
		int64_t start = i * ((resSize + threads - 1) / threads);
		int64_t end = resSize < (i + 1) * ((resSize + threads - 1) / threads)
				    ? resSize : (i + 1) * ((resSize + threads - 1) / threads);

		for(int j = start; j < end; j++){
			//ostringstream tmp;
			//tmp
		    // << refSketch.getReference(buffer[j].refID).name << "\t" 
			// << querySketch.getReference(buffer[j].queryID).name << "\t" 
			// << buffer[j].distance << "\t" 
			// << buffer[j].pValue << "\t"
		    // << buffer[j].number << "/" << buffer[j].denom << endl;

			//oFile.write(tmp.str().c_str(), tmp.str().size());
			if( table && buffer[j].refID == 0 ){
				string tmp = querySketch.getReference(buffer[j].queryID).name;
				if(tmp == "")
					cerr << "WARNING: emptry string name. refID: " << buffer[j].refID 
					     << " queryID: " << buffer[j].queryID << endl;
				oFile.write(tmp.c_str(), tmp.size());
			}
			if( table ) 
			{
				string tmp = '\t' + to_string(buffer[j].distance);
				if( buffer[j].refID == (refSketch.getReferenceCount() - 1) ) 
					tmp += '\n';
				oFile.write(tmp.c_str(), tmp.size());
			} else {

				string tmp =
  	    		       refSketch.getReference(buffer[j].refID).name;
				if( comment ) tmp += ':' + refSketch.getReference(buffer[j].refID).comment;
				tmp += '\t' + querySketch.getReference(buffer[j].queryID).name;
				if( comment ) tmp += ':' + querySketch.getReference(buffer[j].queryID).comment;
				tmp += '\t' 
					+ to_string(buffer[j].distance) + '\t' 
					+ to_string(buffer[j].pValue  ) + '\t'
				    + to_string(buffer[j].number  ) + '/'
					+ to_string(buffer[j].denom   ) + '\n';// endl;
				oFile.write(tmp.c_str(), tmp.size());
			}
		}
		oFile.close();
	}
	//while( resultFile.peek() != EOF )
	//{
	//	resultFile.read((char*)&buffer, sizeof(CommandDistance::Result));
	//	cout 
	//	     << refSketch.getReference(buffer.refID).name << "\t" 
	//		 << querySketch.getReference(buffer.queryID).name << "\t" 
	//		 //<< buffer.refID   << "\t"
	//		 //<< buffer.queryID << "\t"
	//		 << buffer.distance << "\t" 
	//		 << buffer.pValue << "\t"
	//	     << buffer.number << "/" << buffer.denom << endl;
	//}

	//combine and remove tmp files
	double t1 = get_sec();

    if ( table )
    {
        string tmp =  "#query";

        for ( int i = 0; i < refSketch.getReferenceCount(); i++ )
        {
            tmp += '\t' + refSketch.getReference(i).name;
        }
		
		tmp += '\n';
		
		fwrite(tmp.c_str(), 1, tmp.size(), fout);
    }

	int tmpSize = 1<<20;
	unsigned char *tmpBuffer = new unsigned char[tmpSize];

	for(int i = 0; i < threads; i++)
	{
		FILE *tmpFile =  fopen((oFilePrefix + to_string(i)).c_str(), "rb");
		if(tmpFile == NULL)
		{
			cerr << "can not open " << (oFilePrefix + to_string(i)) << endl;
			exit(1);
		}
		int n;
		while(true)
		{
			n = fread(tmpBuffer, sizeof(unsigned char), tmpSize, tmpFile);
			fwrite((void*)tmpBuffer, sizeof(unsigned char), n, fout);
			if( n < tmpSize ) break;
		}

		fclose(tmpFile);
		remove((oFilePrefix + to_string(i)).c_str());
	}

	double t2 = get_sec();
	cerr << "combine time: " << t2 - t1 << endl;

	resultFile.close();
	fclose(fout);

	delete tmpBuffer;
	delete buffer;

	return 0;
}

} //namespace mash

