//output the bin file in text formatting.
//addbyxxm

#include "CommandOutputbin.h"
#include <iostream>
#include <fstream>
#include <string>

#define DIST 1

using namespace::std;

namespace mash{

CommandOutputbin::CommandOutputbin() : Command()
{
	name = "outputbin";
	summary = "output the bin file in text formatting.";
	description = "in the last step of the \" triangle \" we get two binary out file name.bin and dist.bin but both of them are unreadable. The outputbin is used to output the bin file in text formatting.";
	argumentString = "<name.bin> <dist.bin>";
	
	useOption("help");
	//addOption -o outputfile.
}

int CommandOutputbin::run() const
{
	if(arguments.size() < 2 || options.at("help").active)
	{
		print();
		return 0;
	}

	string nameFile = arguments[0];
	string distFile = arguments[1];

	fstream out1(nameFile.c_str(), out1.binary | out1.in);
	fstream out2(distFile.c_str(), out2.binary | out2.in);

	if(!out1.is_open()){
		cerr << "fail to open " << nameFile << endl;
		return 1;
	}
	if(!out2.is_open()){
		cerr << "fail to open " << distFile << endl;
		return 1;
	}

	double d;
	char * name = new char[20];

	int lineSize = 1;
	int lineIndex = 0;
	int readNum = 0;

	out1.read(reinterpret_cast<char*>(name), 20 * sizeof(char));
	for(int i = 0; i < 20; i++){
		if(name[i] < 32) break;//invalid ASCII
		cout << name[i];
	}

	while(1){
		out2.read(reinterpret_cast<char*>(&d), sizeof d);
		if(out2.eof()) break;
		
		readNum++;
		if(d < 1.0){
			if(lineIndex < lineSize){
				cout << '\t' << d;
				lineIndex++;
			}
			else{
				cout << endl;
				out1.read(reinterpret_cast<char*>(name), 20 * sizeof(char));
				for(int i = 0; i < 20; i++){
					if(name[i] < 32) break;//invalid ASCII
					cout << name[i];
				}
				cout << '\t' << d;
				lineIndex = 1;
				lineSize++;
			}
		}

		else{//d >= 1
			for(int j = 0; j < d / 1; j++){
				if(lineIndex < lineSize){
					//cout << DIST << '\t';
					cout << '\t' << DIST;
					lineIndex++;
				}
				else{
					cout << endl;
					out1.read(reinterpret_cast<char*>(name), 20 * sizeof(char));
					for(int i = 0; i < 20; i++){
						if(name[i] < 32) break;//invalid ASCII
						cout << name[i];
					}
					//cout << DIST << '\t';
					cout << '\t' << DIST;

					lineIndex = 1;
					lineSize++;
				}
			}
		} 
			
	}
	cout << endl;//make sure the tail line endl in order to get the same md5-value with Mash.
	cerr << "the readNum is: " << readNum << endl;
	

	delete name;
	out1.close();
	out2.close();

	return 0;
}

} //namespace mash

