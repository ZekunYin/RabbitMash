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
	addOption("outBin", Option(Option::String, "o", "Output", "output results to binary phylip format for better performance. If -o is not specified, results will be redirected to stdout. This will be slow.", ""));
	//addOption("outBin", Option(Option::String, "o", "Output", "output to the binary format with a higher speed.", "./rabbit-mash-output.bin"));
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
	bool outBin; 

   
   	string oFileName = options.at("outBin").argument;

	if(oFileName != ""){
		outBin = true;
		cerr << "writing result to " << oFileName << endl;
	}
	else{
		outBin = false;
	}

	ofstream oFile;

	if(outBin)
	{
		oFile.open(oFileName, ios::out | ios::binary | ios::trunc);
		if(!oFile.is_open()){
			cerr << "Can't opene output file: " << oFileName << endl;
			exit(1);
		}
	}
		

    if ( options.at("pvalue").active || options.at("distance").active )
    {
        edge = true;
    }
   
	//write header to binary
	if( edge )
	{
		uint64_t isEdge = 0x1;
		oFile.write((char *)&isEdge, sizeof(uint64_t));
		oFile.write((char *)&distanceMax, sizeof(double));
		oFile.write((char *)&pValueMax, sizeof(double));
	} else {
		uint64_t isEdge = 0x0;
		oFile.write((char *)&isEdge, sizeof(uint64_t));
		oFile.write((char *)&distanceMax, sizeof(double));
		oFile.write((char *)&pValueMax, sizeof(double));
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
    
    if ( !edge && !outBin)
    {
        cout << '\t' << sketch.getReferenceCount() << endl;
        cout << (comment ? sketch.getReference(0).comment : sketch.getReference(0).name) << endl;
    }
	
    
    ThreadPool<TriangleInput, TriangleOutput> threadPool(compare, threads);
	

    for ( uint64_t i = 1; i < sketch.getReferenceCount(); i++ )
    {
        threadPool.runWhenThreadAvailable(new TriangleInput(sketch, i, parameters, distanceMax, pValueMax));
        
        while ( threadPool.outputAvailable() )
        {
			if(outBin){	
            	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet, oFile);
			}
			else{
            	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet);
			}

        }
    }
    
    while ( threadPool.running() )
    {
		if(outBin){
        	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet, oFile);
		}
		else{
        	writeOutput(threadPool.popOutputWhenAvailable(), comment, edge, pValuePeakToSet);
		}	


    }
    
    if ( !edge )
    {
        cerr << "Max p-value: " << pValuePeakToSet << endl;
    }
    
    if ( warningCount > 0 && ! parameters.reads )
    {
    	warnKmerSize(parameters, *this, lengthMax, lengthMaxName, randomChance, kMin, warningCount);
    }

  	if(outBin){ 
		oFile.close();
	}
    return 0;
}


void CommandTriangle::writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet, ofstream &oFile) const
{
    const Sketch & sketch = output->sketch;
    const Sketch::Reference & ref = sketch.getReference(output->index);

	//CommandDistance::CompareOutput::PairOutput * bufferPair = new CommandDistance::CompareOutput::PairOutput[output->index];
	Result * bufferPair = new Result[output->index];
   	uint64_t passCount = 0; 
    for ( uint64_t i = 0; i < output->index; i++ )
    {
        const CommandDistance::CompareOutput::PairOutput * pair = &output->pairs[i];
		if(pair->pass){
			bufferPair[passCount].refID = output->index;
			bufferPair[passCount].queryID = i;
			bufferPair[passCount].numer = pair->numer;
			bufferPair[passCount].denom = pair->denom;
			bufferPair[passCount].distance = pair->distance;
			bufferPair[passCount].pValue = pair->pValue;
			passCount++;
		}
        if ( pair->pValue > pValuePeakToSet )
        {
            pValuePeakToSet = pair->pValue;
        }
    }

	//oFile.write((char*)bufferPair, output->index * sizeof(CommandDistance::CompareOutput::PairOutput));
	oFile.write((char*)bufferPair, passCount * sizeof(Result));
	oFile.flush();
	delete bufferPair;
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
