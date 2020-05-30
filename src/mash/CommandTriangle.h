// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandTriangle
#define INCLUDED_CommandTriangle

#include "Command.h"
#include "CommandDistance.h"
#include "Sketch.h"

#include <fstream>

namespace mash {

class CommandTriangle : public Command
{
public:
	
	struct Result{
		int refID;
		int queryID;
		double distance;
		double pValue;
		int numer = 0;
		int denom = 0;
	};
    
    struct TriangleInput
    {
        TriangleInput(const Sketch & sketchNew, uint64_t indexNew, const Sketch::Parameters & parametersNew, double maxDistanceNew, double maxPValueNew)
            :
            sketch(sketchNew),
            index(indexNew),
            parameters(parametersNew),
            maxDistance(maxDistanceNew),
            maxPValue(maxPValueNew)
            {}
        
        const Sketch & sketch;
        uint64_t index;
        const Sketch::Parameters & parameters;
        double maxDistance;
        double maxPValue;
    };
    
    struct TriangleOutput
    {
        TriangleOutput(const Sketch & sketchNew, uint64_t indexNew)
            :
            sketch(sketchNew),
            index(indexNew)
        {
            pairs = new CommandDistance::CompareOutput::PairOutput[index];
        }
        
        ~TriangleOutput()
        {
            delete [] pairs;
        }
        
        const Sketch & sketch;
        uint64_t index;
        
        CommandDistance::CompareOutput::PairOutput * pairs;
    };
    
    CommandTriangle();
    
    int run() const; // override
    
private:
    
    double pValueMax;
    bool comment;
    
//void writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet, char * output1Buffer, std::fstream &output1File, double * output2Buffer, std::fstream & output2File) const;
    //void writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet, double * outputBuffer, std::fstream & outputFile) const;
void writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet, std::ofstream &oFile) const;
void writeOutput(TriangleOutput * output, bool comment, bool edge, double & pValuePeakToSet) const;
};

CommandTriangle::TriangleOutput * compare(CommandTriangle::TriangleInput * input);

} // namespace mash

#endif
