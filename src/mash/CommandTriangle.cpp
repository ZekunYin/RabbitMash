// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandDistance.h"
#include "CommandTriangle.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
#include "sketchParameterSetup.h"
#include <math.h>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif


using namespace::std;

namespace mash {

CommandTriangle::CommandTriangle()
: Command()
{
    name = "triangle";
    summary = "Estimate a lower-triangular distance matrix.";
    description = "Estimate the distance of each input sequence to every other input sequence. Outputs a lower-triangular distance matrix in relaxed Phylip format. The input sequences can be fasta or fastq, gzipped or not, or Mash sketch files (.msh) with matching k-mer sizes. Input files can also be files of file names (see -l). If more than one input file is provided, whole files are compared by default (see -i).";
    argumentString = "<seq1> [<seq2>] ...";
    
    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <query> specify paths to sequence files, one per line. The reference file is not affected.", ""));
    addOption("comment", Option(Option::Boolean, "C", "Output", "Use comment fields for sequence names instead of IDs.", ""));
    addOption("edge", Option(Option::Boolean, "E", "Output", "Output edge list instead of Phylip matrix, with fields [seq1, seq2, dist, p-val, shared-hashes].", ""));
    addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report in edge list. Implies -" + getOption("edge").identifier + ".", "1.0", 0., 1.));
    addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report in edge list. Implies -" + getOption("edge").identifier + ".", "1.0", 0., 1.));
	addOption("outBin", Option(Option::String, "o", "Output", "output to the binary format with a higher speed.", "./rabbit-mash-output.bin"));
    //addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
    useSketchOptions();
}

int CommandTriangle::run() const
{
    if ( arguments.size() < 1 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int threads = options.at("threads").getArgumentAsNumber();
    bool list = options.at("list").active;
    //bool log = options.at("log").active;
    bool comment = options.at("comment").active;
    bool edge = options.at("edge").active;
    double pValueMax = options.at("pvalue").getArgumentAsNumber();
    double distanceMax = options.at("distance").getArgumentAsNumber();
    double pValuePeakToSet = 0;
	//bool outBin = options.at("outBin").active;

	//addbyxxm: test the option of string about bool return.
	//if(options.at("outBin").active){
	//	cerr << "the option of -o is true." << endl;
	//}
	//else 
	//	cerr << "the option of -o is true." << endl;
	//return 1;	
   
   	//if(outBin){
   	string ofile = options.at("outBin").argument;

	cerr << "writing result to " << ofile << endl;

	string nameFile = ofile + "Name.bin";
	string distFile = ofile + "Dist.bin";

	//fstream output(ofile, output.binary | output.out);
	fstream output1(nameFile, output1.binary | output1.out);
	fstream output2(distFile, output2.binary | output2.out);
	
	if(!output1.is_open()){
		cerr << "fail to open " << nameFile << endl;
		return 0;
	}

	if(!output2.is_open()){
		cerr << "fail to open " << distFile << endl;
		return 0;
	}
	//}

	



    if ( options.at("pvalue").active || options.at("distance").active )
    {
        edge = true;
    }
    
    Sketch::Parameters parameters;
    
    if ( sketchParameterSetup(parameters, *(Command *)this) )
    {
    	return 1;
    }
    
    if ( arguments.size() == 1 )
    {
    	parameters.concatenated = false;
    }
    
    Sketch sketch;
    
    uint64_t lengthMax;
    double randomChance;
    int kMin;
    string lengthMaxName;
    int warningCount = 0;
    
    vector<string> queryFiles;
    
    for ( int i = 0; i < arguments.size(); i++ )
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
    
    sketch.initFromFiles(queryFiles, parameters);
    
    double lengthThreshold = (parameters.warning * sketch.getKmerSpace()) / (1. - parameters.warning);
    
	for ( uint64_t i = 0; i < sketch.getReferenceCount(); i++ )
	{
		uint64_t length = sketch.getReference(i).length;
		
		if ( length > lengthThreshold )
		{
			if ( warningCount == 0 || length > lengthMax )
			{
				lengthMax = length;
				lengthMaxName = sketch.getReference(i).name;
				randomChance = sketch.getRandomKmerChance(i);
				kMin = sketch.getMinKmerSize(i);
			}
		
			warningCount++;
		}
	}
    
    if ( !edge )
    {
        cout << '\t' << sketch.getReferenceCount() << endl;
        cout << (comment ? sketch.getReference(0).comment : sketch.getReference(0).name) << endl;
    }
	
    
    ThreadPool<TriangleInput, TriangleOutput> threadPool(compare, threads);
   
   	int nameLength = 20;
   	char * output1Buffer = new char[nameLength];
   	double * output2Buffer = new double[sketch.getReferenceCount()];
	//add the first seq into the outputResult.
	//string name0 = sketch.getReference(0).name;
	//output1Buffer = &name0[0];
	//output1.write(reinterpret_cast<char*>(output1Buffer), (20) * sizeof(char));

	

	
	

    for ( uint64_t i = 1; i < sketch.getReferenceCount(); i++ )
    {
        threadPool.runWhenThreadAvailable(new TriangleInput(sketch, i, parameters, distanceMax, pValueMax));
        
        while ( threadPool.outputAvailable() )
        {
		//	if(outBin){	
            	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet, output1Buffer, output1, output2Buffer, output2);
		//	}
		//	else{
        //    	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet);
		//	}

        }
    }
    
    while ( threadPool.running() )
    {
	//	if(outBin){
            writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet, output1Buffer, output1, output2Buffer, output2);
        //	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet, outputBuffer, output);
	//	}
	//	else{
    //    	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet);
	//	}	


    }
    
    if ( !edge )
    {
        cerr << "Max p-value: " << pValuePeakToSet << endl;
    }
    
    if ( warningCount > 0 && ! parameters.reads )
    {
    	warnKmerSize(parameters, *this, lengthMax, lengthMaxName, randomChance, kMin, warningCount);
    }
   
   	delete output1Buffer;
   	delete output2Buffer;

	output1.close();	
	output2.close();	

    return 0;
}

void CommandTriangle::writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet, char * output1Buffer, fstream &output1File, double * output2Buffer, fstream & output2File) const
{
    const Sketch & sketch = output->sketch;
    const Sketch::Reference & ref = sketch.getReference(output->index);

	string name_ = ref.name;
	//cout << name_ << "===addbyxxm" <<endl;
	//exit(0);
	output1Buffer = &(name_[0]);
	output1File.write(reinterpret_cast<char*>(output1Buffer), (20) * sizeof(char));

    
    //if ( !edge )
    //{
    //    cout << (comment ? ref.comment : ref.name);
    //}
    
	double distCount = 0.0;
	int j = 0;
    for ( uint64_t i = 0; i < output->index; i++ )
    {
        const CommandDistance::CompareOutput::PairOutput * pair = &output->pairs[i];
        
        //if ( edge )
        //{
        //    if ( pair->pass )
        //    {
        //        const Sketch::Reference & qry = sketch.getReference(i);
        //        cout << (comment ? ref.comment : ref.name) << '\t'<< (comment ? qry.comment : qry.name) << '\t' << pair->distance << '\t' << pair->pValue << '\t' << pair->numer << '/' << pair->denom << endl;
        //    }
        //}
        //else
        //{
        //    cout << '\t' << pair->distance;
        //}

		//TODO
		if(pair->distance == 1.0){
			distCount += 1.0;
		}
		else{
			if(distCount == 0.0){
				output2Buffer[j] = pair->distance;
				j++;
			}
			else{
				output2Buffer[j] = distCount;
				j++;
				output2Buffer[j] = pair->distance;
				j++;
			}
			distCount = 0.0;
		}
			//if(pair->distance != 1.0){
			//	if(distCount == 0.0){
			//		outputBuffer[j] = pair->distance;
			//		j++;
			//	}
			//	else{
			//		outputBuffer[j] = distCount;
			//		j++;
			//		outputBuffer[j] = pair->distance;
			//		j++;
			//	}
			//	distCount = 0.0;
			//}
			//else{
			//	distCount += 1.0;
			//}
       	//outputBuffer[i] = pair->distance;
		//cout << outputBuffer[i] << endl;

        if ( pair->pValue > pValuePeakToSet )
        {
            pValuePeakToSet = pair->pValue;
        }
    }
	if(distCount != 0){//if all the dist is 1, save the distCount. 
		output2Buffer[j] = distCount;
		j++;
	}

	//outputFile.write(reinterpret_cast<char*>(outputBuffer), output->index * sizeof(double));

	//cout << "the output->index is: " << output->index << endl;
	//cout << "the j is: " << j << endl;

	output2File.write(reinterpret_cast<char*>(output2Buffer), (j) * sizeof(double));
    
    //if ( !edge )
    //{
    //    cout << endl;
    //}
    
    delete output;
}
void CommandTriangle::writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet) const
{
    const Sketch & sketch = output->sketch;
    const Sketch::Reference & ref = sketch.getReference(output->index);
    
    if ( !edge )
    {
        cout << (comment ? ref.comment : ref.name);
    }
    
    for ( uint64_t i = 0; i < output->index; i++ )
    {
        const CommandDistance::CompareOutput::PairOutput * pair = &output->pairs[i];
        
        if ( edge )
        {
            if ( pair->pass )
            {
                const Sketch::Reference & qry = sketch.getReference(i);
                cout << (comment ? ref.comment : ref.name) << '\t'<< (comment ? qry.comment : qry.name) << '\t' << pair->distance << '\t' << pair->pValue << '\t' << pair->numer << '/' << pair->denom << endl;
            }
        }
        else
        {
            cout << '\t' << pair->distance;
        }

		//TODO
       	//outputBuffer[i] = pair->distance;

        if ( pair->pValue > pValuePeakToSet )
        {
            pValuePeakToSet = pair->pValue;
        }
    }

	//outputFile.write(reinterpret_cast<char*>(outputBuffer), output->index * sizeof(double));
    
    if ( !edge )
    {
        cout << endl;
    }
    
    delete output;
}


CommandTriangle::TriangleOutput * compare(CommandTriangle::TriangleInput * input)
{
    const Sketch & sketch = input->sketch;
    
    CommandTriangle::TriangleOutput * output = new CommandTriangle::TriangleOutput(input->sketch, input->index);
    
    uint64_t sketchSize = sketch.getMinHashesPerWindow();
    
    for ( uint64_t i = 0; i < input->index; i++ )
    {
        compareSketches(&output->pairs[i], sketch.getReference(input->index), sketch.getReference(i), sketchSize, sketch.getKmerSize(), sketch.getKmerSpace(), input->maxDistance, input->maxPValue);
    }
    
    return output;
}

} // namespace mash
