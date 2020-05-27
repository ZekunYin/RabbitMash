// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef INCLUDED_CommandDistance
#define INCLUDED_CommandDistance

#include "Command.h"
#include "Sketch.h"
#include <fstream>

namespace mash {

class CommandDistance : public Command
{
public:

   	struct Result{
		int refID;
		int queryID;
		double distance;
		double pValue;
		int number = 0;
		int denom = 0;
	};

    struct CompareInput
    {
        CompareInput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, uint64_t pairCountNew, const Sketch::Parameters & parametersNew, double maxDistanceNew, double maxPValueNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            pairCount(pairCountNew),
            parameters(parametersNew),
            maxDistance(maxDistanceNew),
            maxPValue(maxPValueNew)
            {}
        
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
        
        uint64_t indexRef;
        uint64_t indexQuery;
        uint64_t pairCount;
        
        const Sketch::Parameters & parameters;
        double maxDistance;
        double maxPValue;
    };
    
    struct CompareOutput
    {
        CompareOutput(const Sketch & sketchRefNew, const Sketch & sketchQueryNew, uint64_t indexRefNew, uint64_t indexQueryNew, uint64_t pairCountNew)
            :
            sketchRef(sketchRefNew),
            sketchQuery(sketchQueryNew),
            indexRef(indexRefNew),
            indexQuery(indexQueryNew),
            pairCount(pairCountNew)
        {
            pairs = new PairOutput[pairCount];
        }
        
        ~CompareOutput()
        {
            delete [] pairs;
        }
        
        struct PairOutput
        {
            uint64_t numer;
            uint64_t denom;
            double distance;
            double pValue;
            bool pass;
        };
        
        const Sketch & sketchRef;
        const Sketch & sketchQuery;
        
        uint64_t indexRef;
        uint64_t indexQuery;
        uint64_t pairCount;
        
        PairOutput * pairs;
    };
    
    CommandDistance();
    
    int run() const; // override
    
private:
    
    void writeOutput(CompareOutput * output, bool table, bool comment) const;
    void writeOutput(CompareOutput * output, bool table, bool comment, std::ofstream &) const;
};

CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * input);
void compareSketches(CommandDistance::CompareOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxDistance, double maxPValue);
double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize);

//#if defined (__ICC) || defined (__INTEL_COMPILER)
	#if defined __AVX512F__ && defined __AVX512CD__
uint32_t u32_intersect_vector_avx512(const uint32_t *list1,  uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t *i_a, uint64_t *i_b);
uint64_t u64_intersect_vector_avx512(const uint64_t *list1,  uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3, uint64_t *i_a, uint64_t *i_b);
	#else 
		#if defined __AVX2__
			//implement by avx2

            size_t u32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t* i_a, uint64_t* i_b);

            size_t u64_intersect_vector_avx2(const uint64_t *list1, uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3, uint64_t* i_a, uint64_t* i_b);
		#else
			#if defined __SSE4_1__
			//implement by sse
            uint64_t u64_intersection_vector_sse(const uint64_t *list1, uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3, uint64_t *i_a, uint64_t *i_b);


            uint64_t u32_intersection_vector_sse(const uint32_t *list1, uint64_t size1, const uint32_t *list2, uint64_t size2, uint64_t size3, uint64_t *i_a, uint64_t *i_b);

			#else
			//implement without optimization
			#endif
		#endif
	#endif

//#endif

} // namespace mash

#endif
