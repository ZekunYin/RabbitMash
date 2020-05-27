//covert binary dist output to text format
//yzk 2020/5/28

#include "CommandDumpdist.h"
#include "CommandDistance.h"
#include "Sketch.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace::std;

namespace mash{

CommandDumpdist::CommandDumpdist() : Command()
{
	name = "dumpdist";
	summary = "convert binary dist results to text format and print results to stdout";
	description = "convert binary dist results to text format.";
	argumentString = "<reference.msh> <query.msh> <dist.bin>";
	
	useOption("help");
	//addOption -o outputfile.
}

int CommandDumpdist::run() const
{
	if(arguments.size() < 3 || options.at("help").active)
	{
		print();
		return 0;
	}

	string refMsh   = arguments[0];
	string queryMsh = arguments[1];
	string fileName = arguments[2];

	fstream resultFile(fileName.c_str(), ios::binary | ios::in);
	if(!resultFile.is_open()){
		cerr << "fail to open " << fileName << endl;
		return 1;
	}
	
	Sketch::Parameters parameters;
	//TODO:	setup parameters 
	
	Sketch querySketch;
	Sketch refSketch;

	vector<string> queryFiles;
	vector<string> refFiles;

	queryFiles.push_back(queryMsh);
	refFiles.push_back(refMsh);
	
	querySketch.initFromFiles(queryFiles, parameters);
	refSketch.initFromFiles(refFiles, parameters);
	
	cerr << "ref sketches: "   << refSketch.getReferenceCount()   << endl;
	cerr << "query sketches: " << querySketch.getReferenceCount() << endl;

	CommandDistance::Result buffer;

	while( resultFile.peek() != EOF )
	{
		resultFile.read((char*)&buffer, sizeof(CommandDistance::Result));
		cout 
		     << refSketch.getReference(buffer.refID).name << "\t" 
			 << querySketch.getReference(buffer.queryID).name << "\t" 
			 //<< buffer.refID   << "\t"
			 //<< buffer.queryID << "\t"
			 << buffer.distance << "\t" 
			 << buffer.pValue << "\t"
		     << buffer.number << "/" << buffer.denom << endl;
	}

	resultFile.close();

	return 0;
}

} //namespace mash

