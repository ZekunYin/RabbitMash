// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandDistance.h"
#include "Sketch.h"
#include <iostream>
#include <zlib.h>
#include "ThreadPool.h"
#include "sketchParameterSetup.h"
#include <math.h>

#if defined __AVX512F__ && defined __AVX512CD__
#include <immintrin.h>
#else 
#if defined __AVX2__
#include <immintrin.h>
#else
#if defined __SSE4_1__
#include <nmmintrin.h>
#else

#endif
#endif
#endif

#ifdef USE_BOOST
#include <boost/math/distributions/binomial.hpp>
using namespace::boost::math;
#else
#include <gsl/gsl_cdf.h>
#endif

using namespace::std;

namespace mash {

    CommandDistance::CommandDistance()
        : Command()
    {
        name = "dist";
        summary = "Estimate the distance of query sequences to references.";
        description = "Estimate the distance of each query sequence to the reference. Both the reference and queries can be fasta or fastq, gzipped or not, or Mash sketch files (.msh) with matching k-mer sizes. Query files can also be files of file names (see -l). Whole files are compared by default (see -i). The output fields are [reference-ID, query-ID, distance, p-value, shared-hashes].";
        argumentString = "<reference> <query> [<query>] ...";

        useOption("help");
        addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <query> specify paths to sequence files, one per line. The reference file is not affected.", ""));
        addOption("table", Option(Option::Boolean, "t", "Output", "Table output (will not report p-values, but fields will be blank if they do not meet the p-value threshold).", ""));
        //addOption("log", Option(Option::Boolean, "L", "Output", "Log scale distances and divide by k-mer size to provide a better analog to phylogenetic distance. The special case of zero shared min-hashes will result in a distance of 1.", ""));
        addOption("pvalue", Option(Option::Number, "v", "Output", "Maximum p-value to report.", "1.0", 0., 1.));
        addOption("distance", Option(Option::Number, "d", "Output", "Maximum distance to report.", "1.0", 0., 1.));
        addOption("comment", Option(Option::Boolean, "C", "Output", "Show comment fields with reference/query names (denoted with ':').", "1.0", 0., 1.));
        addOption("binOutput", Option(Option::String, "o", "Output", "Output file name in binary format", ""));
        useSketchOptions();
    }

    int CommandDistance::run() const
    {
        if ( arguments.size() < 2 || options.at("help").active )
        {
            print();
            return 0;
        }
#if defined __AVX512F__ && defined __AVX512CD__
        cerr << "Using AVX512 instructions" << endl;
#else 
#if defined __AVX2__
        cerr << "Using AVX2 instructions" << endl;
        //cerr << "Not implemented yet! Please use:" << endl;
        //cerr << "./configure --disable-simd" << endl;
        //TODO: implement by avx2-----DONE by zjx
#else
#if defined __SSE4_1__
        cerr << "Using SSE4 instructions" << endl;
        //cerr << "Not implemented yet! Please use:" << endl;
        //cerr << "./configure --disable-simd" << endl;
        //TODO:implement by sse----DONE by zjx
#else
        cerr << "No SIMD instructions used" << endl;
        //implement without optimization
#endif
#endif
#endif


        int threads = options.at("threads").getArgumentAsNumber();
        bool list = options.at("list").active;
        bool table = options.at("table").active;
        bool comment = options.at("comment").active;
        //bool log = options.at("log").active;
        double pValueMax = options.at("pvalue").getArgumentAsNumber();
        double distanceMax = options.at("distance").getArgumentAsNumber();
		bool binOut = true;

		string oFileName = options.at("binOutput").argument;

		if (oFileName == "") binOut = false;
		if(binOut)
			cerr << "Results will be written to " << oFileName << endl;

        Sketch::Parameters parameters;

        if ( sketchParameterSetup(parameters, *(Command *)this) )
        {
            return 1;
        }
        //parameters.use64 = false;

        Sketch sketchRef;

        uint64_t lengthMax;
        double randomChance;
        int kMin;
        string lengthMaxName;
        int warningCount = 0;

        const string & fileReference = arguments[0];

		//string oFileName = "/home/ssd/dist_output.bin";
		
		ofstream oFile;
		if(binOut)
		{
			oFile.open(oFileName, ios::out | ios::binary | ios::trunc);
			if(!oFile.is_open()){
				cerr << "Cann't open output file: " << oFileName << endl;
				exit(1);
			}
		}

        bool isSketch = hasSuffix(fileReference, suffixSketch);

        if ( isSketch )
        {
            if ( options.at("kmer").active )
            {
                cerr << "ERROR: The option -" << options.at("kmer").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
                return 1;
            }

            if ( options.at("noncanonical").active )
            {
                cerr << "ERROR: The option -" << options.at("noncanonical").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
                return 1;
            }

            if ( options.at("protein").active )
            {
                cerr << "ERROR: The option -" << options.at("protein").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
                return 1;
            }

            if ( options.at("alphabet").active )
            {
                cerr << "ERROR: The option -" << options.at("alphabet").identifier << " cannot be used when a sketch is provided; it is inherited from the sketch." << endl;
                return 1;
            }
        }
        else
        {
            cerr << "Sketching " << fileReference << " (provide sketch file made with \"mash sketch\" to skip)...";
        }

        vector<string> refArgVector;
        refArgVector.push_back(fileReference);

        //cerr << "Sketch for " << fileReference << " not found or out of date; creating..." << endl;

        sketchRef.initFromFiles(refArgVector, parameters);

        double lengthThreshold = (parameters.warning * sketchRef.getKmerSpace()) / (1. - parameters.warning);

        if ( isSketch )
        {
            if ( options.at("sketchSize").active )
            {
                if ( parameters.reads && parameters.minHashesPerWindow != sketchRef.getMinHashesPerWindow() )
                {
                    cerr << "ERROR: The sketch size must match the reference when using a bloom filter (leave this option out to inherit from the reference sketch)." << endl;
                    return 1;
                }
            }

            parameters.minHashesPerWindow = sketchRef.getMinHashesPerWindow();
            parameters.kmerSize = sketchRef.getKmerSize();
            parameters.noncanonical = sketchRef.getNoncanonical();
            parameters.preserveCase = sketchRef.getPreserveCase();
            parameters.seed = sketchRef.getHashSeed();

            string alphabet;
            sketchRef.getAlphabetAsString(alphabet);
            setAlphabetFromString(parameters, alphabet.c_str());
        }
        else
        {
            for ( uint64_t i = 0; i < sketchRef.getReferenceCount(); i++ )
            {
                uint64_t length = sketchRef.getReference(i).length;

                if ( length > lengthThreshold )
                {
                    if ( warningCount == 0 || length > lengthMax )
                    {
                        lengthMax = length;
                        lengthMaxName = sketchRef.getReference(i).name;
                        randomChance = sketchRef.getRandomKmerChance(i);
                        kMin = sketchRef.getMinKmerSize(i);
                    }

                    warningCount++;
                }
            }

            cerr << "done.\n";
        }

        if ( table )
        {
            cout << "#query";

            for ( int i = 0; i < sketchRef.getReferenceCount(); i++ )
            {
                cout << '\t' << sketchRef.getReference(i).name;
            }

            cout << endl;
        }

        ThreadPool<CompareInput, CompareOutput> threadPool(compare, threads);

        vector<string> queryFiles;

        for ( int i = 1; i < arguments.size(); i++ )
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

        Sketch sketchQuery;

        sketchQuery.initFromFiles(queryFiles, parameters, 0, true);

        uint64_t pairCount = sketchRef.getReferenceCount() * sketchQuery.getReferenceCount();
        uint64_t pairsPerThread = pairCount / parameters.parallelism;

        if ( pairsPerThread == 0 )
        {
            pairsPerThread = 1;
        }

        static uint64_t maxPairsPerThread = 0x1000;

        if ( pairsPerThread > maxPairsPerThread )
        {
            pairsPerThread = maxPairsPerThread;
        }

        uint64_t iFloor = pairsPerThread / sketchRef.getReferenceCount();
        uint64_t iMod = pairsPerThread % sketchRef.getReferenceCount();

        for ( uint64_t i = 0, j = 0; i < sketchQuery.getReferenceCount(); i += iFloor, j += iMod )
        {
            if ( j >= sketchRef.getReferenceCount() )
            {
                if ( i == sketchQuery.getReferenceCount() - 1 )
                {
                    break;
                }

                i++;
                j -= sketchRef.getReferenceCount();
            }

            threadPool.runWhenThreadAvailable(new CompareInput(sketchRef, sketchQuery, j, i, pairsPerThread, parameters, distanceMax, pValueMax));

            while ( threadPool.outputAvailable() )
            {
				if(binOut)
                	writeOutput(threadPool.popOutputWhenAvailable(), table, comment, oFile);
				else
                	writeOutput(threadPool.popOutputWhenAvailable(), table, comment);
            }
        }

        while ( threadPool.running() )
        {
			if(binOut)
            	writeOutput(threadPool.popOutputWhenAvailable(), table, comment, oFile);
			else
            	writeOutput(threadPool.popOutputWhenAvailable(), table, comment);
        }

        if ( warningCount > 0 && ! parameters.reads )
        {
            warnKmerSize(parameters, *this, lengthMax, lengthMaxName, randomChance, kMin, warningCount);
        }

		if(binOut) oFile.close();
        return 0;
    }

    void CommandDistance::writeOutput(CompareOutput * output, bool table, bool comment, ofstream &oFile) const
    {
        uint64_t i = output->indexQuery;
        uint64_t j = output->indexRef;
		Result *buffer = new Result[output->pairCount];
		uint64_t passCount = 0;
		uint64_t k;
        for ( k = 0; k < output->pairCount && i < output->sketchQuery.getReferenceCount(); k++ )
        {
            const CompareOutput::PairOutput * pair = &output->pairs[k];

        //    if ( table && j == 0 )
        //    {
        //        cout << output->sketchQuery.getReference(i).name;
        //    }

            //if ( table )
            //{

            //    if ( pair->pass )
            //    {
			//		buffer[k].refID = j;
			//		buffer[k].queryID = i;
			//		buffer[k].distance = pair->distance;
            //    }
			//	j++;
            //}
            //else 
            if ( pair->pass )
            {
				buffer[passCount].refID = j;
				buffer[passCount].queryID = i;
				buffer[passCount].distance = pair->distance;
				buffer[passCount].pValue = pair->pValue;
				buffer[passCount].number = pair->numer;
				buffer[passCount].denom = pair->denom;
			
				passCount++;

			}

            j++;

			if(j == output->sketchRef.getReferenceCount() )
			{
				j = 0;
				i++;
			}
				
        }

	
		oFile.write((char*)buffer, passCount * sizeof(Result));
		oFile.flush();
		delete buffer;
        delete output;
    }

    CommandDistance::CompareOutput * compare(CommandDistance::CompareInput * input)
    {
        const Sketch & sketchRef = input->sketchRef;
        const Sketch & sketchQuery = input->sketchQuery;

        CommandDistance::CompareOutput * output = new CommandDistance::CompareOutput(input->sketchRef, input->sketchQuery, input->indexRef, input->indexQuery, input->pairCount);

        uint64_t sketchSize = sketchQuery.getMinHashesPerWindow() < sketchRef.getMinHashesPerWindow() ?
            sketchQuery.getMinHashesPerWindow() :
            sketchRef.getMinHashesPerWindow();

        uint64_t i = input->indexQuery;
        uint64_t j = input->indexRef;

        for ( uint64_t k = 0; k < input->pairCount && i < sketchQuery.getReferenceCount(); k++ )
        {
            compareSketches(&output->pairs[k], sketchRef.getReference(j), sketchQuery.getReference(i), sketchSize, sketchRef.getKmerSize(), sketchRef.getKmerSpace(), input->maxDistance, input->maxPValue);

            j++;

            if ( j == sketchRef.getReferenceCount() )
            {
                j = 0;
                i++;
            }
        }

        return output;
    }

	void CommandDistance::writeOutput(CompareOutput * output, bool table, bool comment) const
	{
	    uint64_t i = output->indexQuery;
	    uint64_t j = output->indexRef;
	    
	    for ( uint64_t k = 0; k < output->pairCount && i < output->sketchQuery.getReferenceCount(); k++ )
	    {
	        const CompareOutput::PairOutput * pair = &output->pairs[k];
	        
	        if ( table && j == 0 )
	        {
	            cout << output->sketchQuery.getReference(i).name;
	        }
	        
	        if ( table )
	        {
	            cout << '\t';
	    
	            if ( pair->pass )
	            {
	                cout << pair->distance;
	            }
	        }
	        else if ( pair->pass )
	        {
	            cout << output->sketchRef.getReference(j).name;
	            
	            if ( comment )
	            {
	                cout << ':' << output->sketchRef.getReference(j).comment;
	            }
	            
	            cout << '\t' << output->sketchQuery.getReference(i).name;
	            
	            if ( comment )
	            {
	                cout << ':' << output->sketchQuery.getReference(i).comment;
	            }
	            
	            cout << '\t' << pair->distance << '\t' << pair->pValue << '\t' << pair->numer << '/' << pair->denom << endl;
	        }
	    
	        j++;
	        
	        if ( j == output->sketchRef.getReferenceCount() )
	        {
	            if ( table )
	            {
	                cout << endl;
	            }
	            
	            j = 0;
	            i++;
	        }
	    }
	    
	    delete output;
	}
    void compareSketches(CommandDistance::CompareOutput::PairOutput * output, const Sketch::Reference & refRef, const Sketch::Reference & refQry, uint64_t sketchSize, int kmerSize, double kmerSpace, double maxDistance, double maxPValue)
    {
        uint64_t i = 0;
        uint64_t j = 0;
        uint64_t common = 0;
        uint64_t denom = 0;
        const HashList & hashesSortedRef = refRef.hashesSorted;
        const HashList & hashesSortedQry = refQry.hashesSorted;

        output->pass = false;

        //#if defined (__ICC) || defined (__INTEL_COMPILER)
#if defined __AVX512F__ && defined __AVX512CD__
        //if(parameters.use64)
        if(hashesSortedRef.get64())
        {
            //cerr << "implement the 64bit avx512 addbyxxm " << endl;
            //exit(0);
            common = u64_intersect_vector_avx512((uint64_t*)hashesSortedRef.hashes64.data(), hashesSortedRef.size(), (uint64_t*)hashesSortedQry.hashes64.data(), hashesSortedQry.size(), sketchSize, &i, &j);

            denom = i + j - common;
            //cout << "denom: " << denom << endl;
            //cout << "common: " << common << endl;
            //cout << "i: " << i<< endl;
            //cout << "j: " << j<< endl;
            //cout << "hash i: " << (uint64_t)hashesSortedRef.at(i).hash64 << endl;
            //cout << "hash j: " << (uint64_t)hashesSortedQry.at(j).hash64 << endl;
        }
        else //if(hashesSortedRef.get32())
        {
            //cerr << "implement the 32bit avx512 addbyxxm " << endl;
            //exit(0);
            common = u32_intersect_vector_avx512((uint32_t*)hashesSortedRef.hashes32.data(), hashesSortedRef.size(), (uint32_t*)hashesSortedQry.hashes32.data(), hashesSortedQry.size(), sketchSize, &i, &j);
            denom = i + j - common;
        }
#else
#ifdef __AVX2__
        // implement by avx2
    
        if(hashesSortedRef.get64())
        {
            common = u64_intersect_vector_avx2((uint64_t*)hashesSortedRef.hashes64.data(), hashesSortedRef.size(), (uint64_t*)hashesSortedQry.hashes64.data(), hashesSortedQry.size(), sketchSize, &i, &j);
            //add later
            denom = i + j - common;
        }
        else //if(hashesSortedRef.get32())
        {
            //cout << "hi before call funtion" << endl;
            common = u32_intersect_vector_avx2((uint32_t*)hashesSortedRef.hashes32.data(), hashesSortedRef.size(), (uint32_t*)hashesSortedQry.hashes32.data(), hashesSortedQry.size(), sketchSize, &i, &j);
            denom = i + j - common;
        }
#else
#ifdef __SSE4_1__
        // implement by sse

        if(hashesSortedRef.get64())
        {
            common = u64_intersection_vector_sse((uint64_t*)hashesSortedRef.hashes64.data(), hashesSortedRef.size(), (uint64_t*)hashesSortedQry.hashes64.data(), hashesSortedQry.size(), sketchSize, &i, &j);
            //add later
            denom = i + j - common;
        }
        else //if(hashesSortedRef.get32())
        {
            //cout << "hi before call funtion" << endl;
            common = u32_intersection_vector_sse((uint32_t*)hashesSortedRef.hashes32.data(), hashesSortedRef.size(), (uint32_t*)hashesSortedQry.hashes32.data(), hashesSortedQry.size(), sketchSize, &i, &j);
            denom = i + j - common;
        }
#else
        //implement without optimization	
        while ( denom < sketchSize && i < hashesSortedRef.size() && j < hashesSortedQry.size() )
        {
            if ( hashLessThan(hashesSortedRef.at(i), hashesSortedQry.at(j), hashesSortedRef.get64()) )
            {
                i++;
            }
            else if ( hashLessThan(hashesSortedQry.at(j), hashesSortedRef.at(i), hashesSortedRef.get64()) )
            {
                j++;
            }
            else
            {
                //		cout << "res: " << (uint64_t)hashesSortedRef.at(i).hash64 << endl;
                i++;
                j++;
                common++;
            }

            denom++;
        }
#endif
#endif
#endif

        //#endif


        if ( denom < sketchSize )
        {
            // complete the union operation if possible

            if ( i < hashesSortedRef.size() )
            {
                denom += hashesSortedRef.size() - i;
            }

            if ( j < hashesSortedQry.size() )
            {
                denom += hashesSortedQry.size() - j;
            }

            if ( denom > sketchSize )
            {
                denom = sketchSize;
            }
        }

        double distance;
        double jaccard = double(common) / denom;

        if ( common == denom ) // avoid -0
        {
            distance = 0;
        }
        else if ( common == 0 ) // avoid inf
        {
            distance = 1.;
        }
        else
        {
            //distance = log(double(common + 1) / (denom + 1)) / log(1. / (denom + 1));
            distance = -log(2 * jaccard / (1. + jaccard)) / kmerSize;

            if ( distance > 1 )
            {
                distance = 1;
            }
        }

        if ( maxDistance >= 0 && distance > maxDistance )
        {
            return;
        }

        output->numer = common;
        output->denom = denom;
        output->distance = distance;
        output->pValue = pValue(common, refRef.length, refQry.length, kmerSpace, denom);

        if ( maxPValue >= 0 && output->pValue > maxPValue )
        {
            return;
        }

        output->pass = true;


    }

    double pValue(uint64_t x, uint64_t lengthRef, uint64_t lengthQuery, double kmerSpace, uint64_t sketchSize)
    {
        if ( x == 0 )
        {
            return 1.;
        }

        double pX = 1. / (1. + kmerSpace / lengthRef);
        double pY = 1. / (1. + kmerSpace / lengthQuery);

        double r = pX * pY / (pX + pY - pX * pY);

        //double M = (double)kmerSpace * (pX + pY) / (1. + r);

        //return gsl_cdf_hypergeometric_Q(x - 1, r * M, M - r * M, sketchSize);

#ifdef USE_BOOST
        return cdf(complement(binomial(sketchSize, r), x - 1));
#else
        return gsl_cdf_binomial_Q(x - 1, r, sketchSize);
#endif
    }

    //#if defined (__ICC) || defined (__INTEL_COMPILER)

    uint64_t u64_intersect_scalar_stop(const uint64_t *list1, uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3,
            uint64_t *i_a, uint64_t *i_b){
        uint64_t counter=0;
        const uint64_t *end1 = list1+size1, *end2 = list2+size2;
        *i_a = 0;
        *i_b = 0;
        //uint64_t stop = 0;
        // hard to get only the loop instructions, now only a tiny check at the top wrong
#if IACA_INTERSECT_SCALAR
        IACA_START
#endif
            while(list1 != end1 && list2 != end2 ){
                if(*list1 < *list2){
                    list1++;
                    (*i_a)++;
                    size3--;
                }else if(*list1 > *list2){
                    list2++; 
                    (*i_b)++;
                    size3--;
                }else{
                    //result[counter++] = *list1;
                    counter++;
                    list1++; list2++; 
                    (*i_a)++;
                    (*i_b)++;
                    size3--;
                }
                if(size3 == 0) break;
            }
#if IACA_INTERSECT_SCALAR
        IACA_END
#endif
            return counter;
    }

    uint64_t u32_intersect_scalar_stop(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3,
            uint64_t *i_a, uint64_t *i_b){
        uint64_t counter=0;
        const uint32_t *end1 = list1+size1, *end2 = list2+size2;
        *i_a = 0;
        *i_b = 0;
        //uint64_t stop = 0;
        // hard to get only the loop instructions, now only a tiny check at the top wrong
#if IACA_INTERSECT_SCALAR
        IACA_START
#endif
            while(list1 != end1 && list2 != end2 ){
                if(*list1 < *list2){
                    list1++;
                    (*i_a)++;
                    size3--;
                }else if(*list1 > *list2){
                    list2++; 
                    (*i_b)++;
                    size3--;
                }else{
                    //result[counter++] = *list1;
                    counter++;
                    list1++; list2++; 
                    (*i_a)++;
                    (*i_b)++;
                    size3--;
                }
                if(size3 == 0) break;
            }
#if IACA_INTERSECT_SCALAR
        IACA_END
#endif
            return counter;
    }
#if defined __AVX512F__ && defined __AVX512CD__

    static /*constexpr*/ std::array<uint64_t,8*7> u64_prepare_shuffle_vectors(){
        std::array<uint64_t,8*7> arr = {};
        uint64_t start=1;
        for(uint64_t i=0; i<7; ++i){
            uint64_t counter = start;
            for(uint64_t j=0; j<8; ++j){
                arr[i*8 + j] = counter % 8;
                ++counter;
            }
            ++start;
        }
        return arr;
    }
    static const /*constexpr*/ auto u64_shuffle_vectors_arr = u64_prepare_shuffle_vectors();

    static const /*constexpr*/ __m512i *u64_shuffle_vectors = (__m512i*)u64_shuffle_vectors_arr.data();
    static void inline
        inspect(__m512i v){
            uint64_t f[8] __attribute__((aligned(64)));
            _mm512_store_epi64(f,v);
            printf("[%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld]\n", f[0], f[1], f[2], f[3],f[4], f[5], f[6], f[7]);
        }
    uint64_t u64_intersect_vector_avx512(const uint64_t *list1,  uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3, uint64_t *i_a, uint64_t *i_b)
    {
        //assert(size3 <= size1 + size2);
        uint64_t count=0;
        *i_a = 0;
        *i_b = 0;
        uint64_t st_a = (size1 / 8) * 8;
        uint64_t st_b = (size2 / 8) * 8;
        //	uint64_t stop = (size3 / 16) * 16;

        uint64_t i_a_s, i_b_s;

        if(size3 <= 8){
            count += u64_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
            return count;
        }

        uint64_t stop = size3 - 8;
        //cout << "stop: " << stop <<  endl;
        //__m512i sv0  = u64_shuffle_vectors[ 0];//_mm512_set_epi32(0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);
        //__m512i sv1  = u64_shuffle_vectors[ 1];//_mm512_set_epi32(1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2);
        //__m512i sv2  = u64_shuffle_vectors[ 2];//_mm512_set_epi32(2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3);
        //__m512i sv3  = u64_shuffle_vectors[ 3];//_mm512_set_epi32(3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4);
        //__m512i sv4  = u64_shuffle_vectors[ 4];//_mm512_set_epi32(4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5);
        //__m512i sv5  = u64_shuffle_vectors[ 5];//_mm512_set_epi32(5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6);
        //__m512i sv6  = u64_shuffle_vectors[ 6];//_mm512_set_epi32(6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7);

        __m512i sv0  = _mm512_set_epi64(0,7,6,5,4,3,2,1);
        __m512i sv1  = _mm512_set_epi64(1,0,7,6,5,4,3,2);
        __m512i sv2  = _mm512_set_epi64(2,1,0,7,6,5,4,3);
        __m512i sv3  = _mm512_set_epi64(3,2,1,0,7,6,5,4);
        __m512i sv4  = _mm512_set_epi64(4,3,2,1,0,7,6,5);
        __m512i sv5  = _mm512_set_epi64(5,4,3,2,1,0,7,6);
        __m512i sv6  = _mm512_set_epi64(6,5,4,3,2,1,0,7);

        //__m512i vzero = _mm512_setzero_epi32();
        while(*i_a < st_a && *i_b < st_b){
            //__m512i v_a = _mm512_loadu_epi64((__m512i*)&list1[*i_a]);
            //__m512i v_b = _mm512_loadu_epi64((__m512i*)&list2[*i_b]);
            __m512i v_a = _mm512_loadu_si512((__m512i*)&list1[*i_a]);
            __m512i v_b = _mm512_loadu_si512((__m512i*)&list2[*i_b]);

            uint64_t a_max = list1[*i_a+7];
            uint64_t b_max = list2[*i_b+7];
            //cout << "a " << *i_a << ": " << list1[*i_a+7] << endl;
            //cout << "b " << *i_b << ": " << list2[*i_b+7] << endl;
            //if(a_max <= b_max)
            //	cout << "choose a" << endl;
            //else
            //	cout << "choose b" << endl;
            //cout << endl;

            *i_a += (a_max <= b_max) * 8;
            *i_b += (a_max >= b_max) * 8;

            __mmask16 cmp0 = _mm512_cmpeq_epu64_mask(v_a, v_b);
            __m512i rot0 = _mm512_permutexvar_epi64(sv0, v_b);
            __mmask16 cmp1 = _mm512_cmpeq_epu64_mask(v_a, rot0);
            __m512i rot1 = _mm512_permutexvar_epi64(sv1, v_b);
            __mmask16 cmp2 = _mm512_cmpeq_epu64_mask(v_a, rot1);
            __m512i rot2 = _mm512_permutexvar_epi64(sv2, v_b);
            __mmask16 cmp3 = _mm512_cmpeq_epu64_mask(v_a, rot2);
            cmp0 = _mm512_kor(_mm512_kor(cmp0, cmp1), _mm512_kor(cmp2, cmp3));

            __m512i rot3 = _mm512_permutexvar_epi64(sv3, v_b);
            __mmask16 cmp4 = _mm512_cmpeq_epu64_mask(v_a, rot3);
            __m512i rot4 = _mm512_permutexvar_epi64(sv4, v_b);
            __mmask16 cmp5 = _mm512_cmpeq_epu64_mask(v_a, rot4);
            __m512i rot5 = _mm512_permutexvar_epi64(sv5, v_b);
            __mmask16 cmp6 = _mm512_cmpeq_epu64_mask(v_a, rot5);
            __m512i rot6 = _mm512_permutexvar_epi64(sv6, v_b);
            __mmask16 cmp7 = _mm512_cmpeq_epu64_mask(v_a, rot6);
            cmp4 = _mm512_kor(_mm512_kor(cmp4, cmp5), _mm512_kor(cmp6, cmp7));


            cmp0 = _mm512_kor(cmp0, cmp4);

            //_mm512_mask_compressstoreu_epi64(&result[count], cmp0, v_a);
            //__m512i vres = _mm512_mask_compress_epi64(_mm512_setzero_epi32(), cmp0, v_a);
            //if(cmp0 > 0)
            //	inspect(vres);
            count += _mm_popcnt_u64(cmp0);
            if(*i_a + *i_b - count >= stop){
                count -= _mm_popcnt_u64(cmp0);
                *i_a -= (a_max <= b_max) * 8;
                *i_b -= (a_max >= b_max) * 8;
                break;
            }

        }
        //cout << "avx512 i_a: " << *i_a << endl;
        //cout << "avx512 i_b: " << *i_b << endl;
        //cout << "avx512 count " << count << endl;
        //cout << "avx512 list[i_a]: " << list1[*i_a] << endl;
        //cout << "avx512 list[i_b]: " << list2[*i_b] << endl;
        // intersect the tail using scalar intersection
        //count += u64_intersect_scalar(list1+i_a, size1-i_a, list2+i_b, size2-i_b, result+count);
        //if(size3 - *i_a - *i_b == 0){
        //	return count;
        //}
        //else{
        count += u64_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3 - (*i_a+*i_b - count), &i_a_s, &i_b_s);

        *i_a += i_a_s;
        *i_b += i_b_s;
        //}
        return count;

    }


    static /*constexpr*/ std::array<uint32_t,16*16> u32_prepare_shuffle_vectors(){
        std::array<uint32_t,16*16> arr_unalign = {};
        //std::array<uint64_t,8*7> arr = {};
        uint32_t *arr = (uint32_t *)(((long)arr_unalign.data() + 64) & (~63));
        //__m512i *temp;
        uint64_t start=1;
        for(uint64_t i=0; i<15; ++i){
            uint64_t counter = start;
            for(uint64_t j=0; j<16; ++j){
                arr[i*16 + j] = counter % 16;
                ++counter;
            }
            ++start;
        }
        return arr_unalign;
    }
    static const /*constexpr*/ auto u32_shuffle_vectors_arr = u32_prepare_shuffle_vectors();
    static const /*constexpr*/ __m512i *u32_shuffle_vectors_unalign = (__m512i*)u32_shuffle_vectors_arr.data();
    //static const /*constexpr*/ __m512i *u64_shuffle_vectors = (__m512i*)u64_shuffle_vectors_arr.data();
    static const __m512i *u32_shuffle_vectors = (__m512i *)(((long)u32_shuffle_vectors_unalign + 64) & (~63));
    //size3 is the stop threshold of the sum of size1&size2 
    uint32_t u32_intersect_vector_avx512(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t *i_a, uint64_t *i_b){
        //assert(size3 <= size1 + size2);
        uint64_t count=0;
        *i_a = 0;
        *i_b = 0;
        uint64_t st_a = (size1 / 16) * 16;
        uint64_t st_b = (size2 / 16) * 16;
        //	uint64_t stop = (size3 / 16) * 16;

        uint64_t i_a_s, i_b_s;

        if(size3 <= 16){
            count += u32_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
            return count;
        }

        uint64_t stop = size3 - 16;
        //cout << "stop: " << stop <<  endl;
        __m512i sv0   = _mm512_set_epi32(0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1); //u32_shuffle_vectors[0 ];
        __m512i sv1   = _mm512_set_epi32(1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2); //u32_shuffle_vectors[1 ];
        __m512i sv2   = _mm512_set_epi32(2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3); //u32_shuffle_vectors[2 ];
        __m512i sv3   = _mm512_set_epi32(3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4); //u32_shuffle_vectors[3 ];
        __m512i sv4   = _mm512_set_epi32(4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5); //u32_shuffle_vectors[4 ];
        __m512i sv5   = _mm512_set_epi32(5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6); //u32_shuffle_vectors[5 ];
        __m512i sv6   = _mm512_set_epi32(6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7); //u32_shuffle_vectors[6 ];
        __m512i sv7   = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8); //u32_shuffle_vectors[7 ];
        __m512i sv8   = _mm512_set_epi32(8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9); //u32_shuffle_vectors[8 ];
        __m512i sv9   = _mm512_set_epi32(9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10); //u32_shuffle_vectors[9 ];
        __m512i sv10  = _mm512_set_epi32(10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11); //u32_shuffle_vectors[10];
        __m512i sv11  = _mm512_set_epi32(11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12); //u32_shuffle_vectors[11];
        __m512i sv12  = _mm512_set_epi32(12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13); //u32_shuffle_vectors[12];
        __m512i sv13  = _mm512_set_epi32(13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14); //u32_shuffle_vectors[13];
        __m512i sv14  = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15); //u32_shuffle_vectors[14];
        //__m512i vzero = _mm512_setzero_epi32();
        while(*i_a < st_a && *i_b < st_b){



            uint32_t a_max = list1[*i_a+15];
            uint32_t b_max = list2[*i_b+15];

            //cout << a_max << endl;
            //cout << b_max << endl;

            __m512i v_a = _mm512_loadu_si512((__m512i*)&(list1[*i_a]));
            __m512i v_b = _mm512_loadu_si512((__m512i*)&(list2[*i_b]));
            //cout << "a " << *i_a << ": " << list1[*i_a+7] << endl;
            //cout << "b " << *i_b << ": " << list2[*i_b+7] << endl;
            //if(a_max <= b_max)
            //	cout << "choose a" << endl;
            //else
            //	cout << "choose b" << endl;
            //cout << endl;

            *i_a += (a_max <= b_max) * 16;
            *i_b += (a_max >= b_max) * 16;

            __mmask16 cmp0 = _mm512_cmpeq_epu32_mask(v_a, v_b);
            __m512i rot0 = _mm512_permutexvar_epi32(sv0, v_b);
            __mmask16 cmp1 = _mm512_cmpeq_epu32_mask(v_a, rot0);
            __m512i rot1 = _mm512_permutexvar_epi32(sv1, v_b);
            __mmask16 cmp2 = _mm512_cmpeq_epu32_mask(v_a, rot1);
            __m512i rot2 = _mm512_permutexvar_epi32(sv2, v_b);
            __mmask16 cmp3 = _mm512_cmpeq_epu32_mask(v_a, rot2);
            cmp0 = _mm512_kor(_mm512_kor(cmp0, cmp1), _mm512_kor(cmp2, cmp3));

            __m512i rot3 = _mm512_permutexvar_epi32(sv3, v_b);
            __mmask16 cmp4 = _mm512_cmpeq_epu32_mask(v_a, rot3);
            __m512i rot4 = _mm512_permutexvar_epi32(sv4, v_b);
            __mmask16 cmp5 = _mm512_cmpeq_epu32_mask(v_a, rot4);
            __m512i rot5 = _mm512_permutexvar_epi32(sv5, v_b);
            __mmask16 cmp6 = _mm512_cmpeq_epu32_mask(v_a, rot5);
            __m512i rot6 = _mm512_permutexvar_epi32(sv6, v_b);
            __mmask16 cmp7 = _mm512_cmpeq_epu32_mask(v_a, rot6);
            cmp4 = _mm512_kor(_mm512_kor(cmp4, cmp5), _mm512_kor(cmp6, cmp7));

            __m512i rot7 = _mm512_permutexvar_epi32(sv7, v_b);
            cmp1 = _mm512_cmpeq_epu32_mask(v_a, rot7);
            __m512i rot8 = _mm512_permutexvar_epi32(sv8, v_b);
            cmp2 = _mm512_cmpeq_epu32_mask(v_a, rot8);
            __m512i rot9 = _mm512_permutexvar_epi32(sv9, v_b);
            cmp3 = _mm512_cmpeq_epu32_mask(v_a, rot9);
            __m512i rot10 = _mm512_permutexvar_epi32(sv10, v_b);
            cmp5 = _mm512_cmpeq_epu32_mask(v_a, rot10);
            cmp1 = _mm512_kor(_mm512_kor(cmp1, cmp2), _mm512_kor(cmp3, cmp5));

            __m512i rot11 = _mm512_permutexvar_epi32(sv11, v_b);
            cmp2 = _mm512_cmpeq_epu32_mask(v_a, rot11);
            __m512i rot12 = _mm512_permutexvar_epi32(sv12, v_b);
            cmp3 = _mm512_cmpeq_epu32_mask(v_a, rot12);
            __m512i rot13 = _mm512_permutexvar_epi32(sv13, v_b);
            cmp5 = _mm512_cmpeq_epu32_mask(v_a, rot13);
            __m512i rot14 = _mm512_permutexvar_epi32(sv14, v_b);
            cmp6 = _mm512_cmpeq_epu32_mask(v_a, rot14);
            cmp2 = _mm512_kor(_mm512_kor(cmp2, cmp3), _mm512_kor(cmp5, cmp6));


            cmp0 = _mm512_kor(_mm512_kor(cmp0, cmp4), _mm512_kor(cmp1, cmp2));




            //_mm512_mask_compressstoreu_epi64(&result[count], cmp0, v_a);
            //__m512i vres = _mm512_mask_compress_epi64(_mm512_setzero_epi32(), cmp0, v_a);
            //if(cmp0 > 0)
            //	inspect(vres);
            count += _mm_popcnt_u32(cmp0);
            if(*i_a + *i_b - count >= stop){
                count -= _mm_popcnt_u32(cmp0);
                *i_a -= (a_max <= b_max) * 16;
                *i_b -= (a_max >= b_max) * 16;
                break;
            }

        }
        //cout << "avx512 i_a: " << *i_a << endl;
        //cout << "avx512 i_b: " << *i_b << endl;
        //cout << "avx512 count " << count << endl;
        //cout << "avx512 list[i_a]: " << list1[*i_a] << endl;
        //cout << "avx512 list[i_b]: " << list2[*i_b] << endl;
        // intersect the tail using scalar intersection
        //count += u64_intersect_scalar(list1+i_a, size1-i_a, list2+i_b, size2-i_b, result+count);
        //if(size3 - *i_a - *i_b == 0){
        //	return count;
        //}
        //else{
        count += u32_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3 - (*i_a+*i_b - count), &i_a_s, &i_b_s);

        *i_a += i_a_s;
        *i_b += i_b_s;
        //}
        return count;
    }



#else
#ifdef __AVX2__
    
    size_t u32_intersect_vector_avx2(const uint32_t *list1, uint32_t size1, const uint32_t *list2, uint32_t size2, uint32_t size3, uint64_t* i_a, uint64_t* i_b){
        //assert(size3 <= size1 + size2);
        uint64_t count=0;
        *i_a = 0;
        *i_b = 0;
        uint64_t st_a = (size1 / 8) * 8;
        uint64_t st_b = (size2 / 8) * 8;
        //	uint64_t stop = (size3 / 16) * 16;

        uint64_t i_a_s, i_b_s;

        if(size3 <= 16){
            count += u32_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
            return count;
        }

        uint64_t stop = size3 - 16;
        while(*i_a < st_a && *i_b < st_b){



            uint32_t a_max = list1[*i_a+7];
            uint32_t b_max = list2[*i_b+7];

            //cout << a_max << endl;
            //cout << b_max << endl;

            __m256i v_a = _mm256_loadu_si256((__m256i*)&(list1[*i_a]));
            __m256i v_b = _mm256_loadu_si256((__m256i*)&(list2[*i_b]));
            //cout << "a " << *i_a << ": " << list1[*i_a+7] << endl;
            //cout << "b " << *i_b << ": " << list2[*i_b+7] << endl;
            //if(a_max <= b_max)
            //	cout << "choose a" << endl;
            //else
            //	cout << "choose b" << endl;
            //cout << endl;

            *i_a += (a_max <= b_max) * 8;
            *i_b += (a_max >= b_max) * 8;


            /*constexpr*/ const int32_t cyclic_shift = _MM_SHUFFLE(0,3,2,1); //rotating right
            /*constexpr*/ const int32_t cyclic_shift2= _MM_SHUFFLE(2,1,0,3); //rotating left
            /*constexpr*/ const int32_t cyclic_shift3= _MM_SHUFFLE(1,0,3,2); //between
            __m256i cmp_mask1 = _mm256_cmpeq_epi32(v_a, v_b);
            __m256 rot1 = _mm256_permute_ps((__m256)v_b, cyclic_shift);
            __m256i cmp_mask2 = _mm256_cmpeq_epi32(v_a, (__m256i)rot1);
            __m256 rot2 = _mm256_permute_ps((__m256)v_b, cyclic_shift3);
            __m256i cmp_mask3 = _mm256_cmpeq_epi32(v_a, (__m256i)rot2);
            __m256 rot3 = _mm256_permute_ps((__m256)v_b, cyclic_shift2);
            __m256i cmp_mask4 = _mm256_cmpeq_epi32(v_a, (__m256i)rot3);

            __m256 rot4 = _mm256_permute2f128_ps((__m256)v_b, (__m256)v_b, 1);

            __m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a, (__m256i)rot4);
            __m256 rot5 = _mm256_permute_ps(rot4, cyclic_shift);
            __m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a, (__m256i)rot5);
            __m256 rot6 = _mm256_permute_ps(rot4, cyclic_shift3);
            __m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a, (__m256i)rot6);
            __m256 rot7 = _mm256_permute_ps(rot4, cyclic_shift2);
            __m256i cmp_mask8 = _mm256_cmpeq_epi32(v_a, (__m256i)rot7);

            __m256i cmp_mask = _mm256_or_si256(
                    _mm256_or_si256(
                        _mm256_or_si256(cmp_mask1, cmp_mask2),
                        _mm256_or_si256(cmp_mask3, cmp_mask4)
                        ),
                    _mm256_or_si256(
                        _mm256_or_si256(cmp_mask5, cmp_mask6),
                        _mm256_or_si256(cmp_mask7, cmp_mask8)
                        )
                    );
            int32_t mask = _mm256_movemask_ps((__m256)cmp_mask);

            //__m256i idx = _mm256_load_si256((const __m256i*)&shuffle_mask_avx[mask*8]);
            //__m256i p = _mm256_permutevar8x32_epi32(v_a, idx);
            //_mm256_storeu_si256((__m256i*)&result[count], p);
            //_mm256_storeu_si256((__m256i*)&result[count], v_a);
            //_mm256_storeu_si256((__m256i*)&result[count+1], v_b);
            //_mm256_storeu_si256((__m256i*)&result[count+2], v_a2);
            //_mm256_storeu_si256((__m256i*)&result[count+3], v_b2);
            count += _mm_popcnt_u32(mask);

            //_mm512_mask_compressstoreu_epi64(&result[count], cmp0, v_a);
            //__m512i vres = _mm512_mask_compress_epi64(_mm512_setzero_epi32(), cmp0, v_a);
            //if(cmp0 > 0)
            //	inspect(vres);

            if(*i_a + *i_b - count >= stop){
                //count -= _mm_popcnt_u32(cmp0);
                //*i_a -= (a_max <= b_max) * 16;
                //*i_b -= (a_max >= b_max) * 16;
                break;
            }

        }
        //cout << "avx512 i_a: " << *i_a << endl;
        //cout << "avx512 i_b: " << *i_b << endl;
        //cout << "avx512 count " << count << endl;
        //cout << "avx512 list[i_a]: " << list1[*i_a] << endl;
        //cout << "avx512 list[i_b]: " << list2[*i_b] << endl;
        // intersect the tail using scalar intersection
        //count += u64_intersect_scalar(list1+i_a, size1-i_a, list2+i_b, size2-i_b, result+count);
        //if(size3 - *i_a - *i_b == 0){
        //	return count;
        //}
        //else{
        count += u32_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3 - (*i_a+*i_b - count), &i_a_s, &i_b_s);

        *i_a += i_a_s;
        *i_b += i_b_s;
        //}
        return count;
    }

    size_t u64_intersect_vector_avx2(const uint64_t *list1, uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3, uint64_t* i_a, uint64_t* i_b){
    		//assert(size3 <= size1 + size2);
    		uint64_t count=0;
    		*i_a = 0;
    		*i_b = 0;
    		uint64_t st_a = (size1 / 4) * 4;
    		uint64_t st_b = (size2 / 4) * 4;
    	//	uint64_t stop = (size3 / 16) * 16;
    
    		uint64_t i_a_s, i_b_s;
    
    		if(size3 <= 8){
    			count += u64_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
    			return count;
    		}
    		
    		uint64_t stop = size3 - 8;
            //cout << "stop:" << stop << endl;
    		while(*i_a < st_a && *i_b < st_b){
    				
    
    
    				uint64_t a_max = list1[*i_a+3];
    				uint64_t b_max = list2[*i_b+3];
    			
    				//cout << a_max << endl;
    				//cout << b_max << endl;
    
    				__m256i v_a = _mm256_loadu_si256((__m256i*)&(list1[*i_a]));
    				__m256i v_b = _mm256_loadu_si256((__m256i*)&(list2[*i_b]));
    				//cout << "a " << *i_a << ": " << list1[*i_a+7] << endl;
    				//cout << "b " << *i_b << ": " << list2[*i_b+7] << endl;
    				//if(a_max <= b_max)
    				//	cout << "choose a" << endl;
    				//else
    				//	cout << "choose b" << endl;
    				//cout << endl;
    
    				*i_a += (a_max <= b_max) * 4;
    				*i_b += (a_max >= b_max) * 4;
    
    
    				/*constexpr*/// const int64_t cyclic_shift = _MM_SHUFFLE(0,3,2,1); //rotating right
    				/*constexpr*/// const int64_t cyclic_shift2= _MM_SHUFFLE(2,1,0,3); //rotating left
    				/*constexpr*/// const int64_t cyclic_shift3= _MM_SHUFFLE(1,0,3,2); //between
    				__m256i cmp_mask1 = _mm256_cmpeq_epi64(v_a, v_b);
    				__m256i rot1 = _mm256_permute4x64_epi64(v_b, 57);//00111001
    				__m256i cmp_mask2 = _mm256_cmpeq_epi64(v_a, (__m256i)rot1);
    				__m256i rot2 = _mm256_permute4x64_epi64(v_b, 78);//01001110
    				__m256i cmp_mask3 = _mm256_cmpeq_epi64(v_a, (__m256i)rot2);
    				__m256i rot3 = _mm256_permute4x64_epi64(v_b, 147);//10010011
    				__m256i cmp_mask4 = _mm256_cmpeq_epi64(v_a, (__m256i)rot3);
                    //inspect(v_a);
                    //inspect(v_b);
                 
                    //inspect(rot1);
                    //inspect(rot2);
                    //inspect(rot3);
                    //exit(0);
    				//__m256 rot4 = _mm256_permute2f128_ps((__m256)v_b, (__m256)v_b, 1);
    
    				//__m256i cmp_mask5 = _mm256_cmpeq_epi32(v_a, (__m256i)rot4);
    				//__m256 rot5 = _mm256_permute_ps(rot4, cyclic_shift);
    				//__m256i cmp_mask6 = _mm256_cmpeq_epi32(v_a, (__m256i)rot5);
    				//__m256 rot6 = _mm256_permute_ps(rot4, cyclic_shift3);
    				//__m256i cmp_mask7 = _mm256_cmpeq_epi32(v_a, (__m256i)rot6);
    				//__m256 rot7 = _mm256_permute_ps(rot4, cyclic_shift2);
    				//__m256i cmp_mask8 = _mm256_cmpeq_epi32(v_a, (__m256i)rot7);
    
    				__m256i cmp_mask = _mm256_or_si256(
    								//_mm256_or_si256(
    										_mm256_or_si256(cmp_mask1, cmp_mask2),
    										_mm256_or_si256(cmp_mask3, cmp_mask4)
    							//			),
    							//	_mm256_or_si256(
    							//			_mm256_or_si256(cmp_mask5, cmp_mask6),
    							//			_mm256_or_si256(cmp_mask7, cmp_mask8)
    								//		)
    								);
    				//inspect(cmp_mask);
                    int64_t mask = _mm256_movemask_pd((__m256d)cmp_mask);
                    //cout << "mask:"<< mask << endl;
    		        //printf("[%x,%x,%x,%x]\n", f[0], f[1], f[2], f[3]);
    				
                    //__m256i idx = _mm256_load_si256((const __m256i*)&shuffle_mask_avx[mask*8]);
    				//__m256i p = _mm256_permutevar8x32_epi32(v_a, idx);
    				//_mm256_storeu_si256((__m256i*)&result[count], p);
    				//_mm256_storeu_si256((__m256i*)&result[count], v_a);
    				//_mm256_storeu_si256((__m256i*)&result[count+1], v_b);
    				//_mm256_storeu_si256((__m256i*)&result[count+2], v_a2);
    				//_mm256_storeu_si256((__m256i*)&result[count+3], v_b2);
    				count += _mm_popcnt_u64(mask);
                    //cout << count <<endl;
    				//_mm512_mask_compressstoreu_epi64(&result[count], cmp0, v_a);
    				//__m512i vres = _mm512_mask_compress_epi64(_mm512_setzero_epi32(), cmp0, v_a);
    				//if(cmp0 > 0)
    				//	inspect(vres);
    				if(*i_a + *i_b - count >= stop){
    					//count -= _mm_popcnt_u32(cmp0);
    					//*i_a -= (a_max <= b_max) * 16;
    					//*i_b -= (a_max >= b_max) * 16;
    					break;
    				}
    
    		}
    		//cout << "avx512 i_a: " << *i_a << endl;
    		//cout << "avx512 i_b: " << *i_b << endl;
    		//cout << "avx512 count " << count << endl;
    		//cout << "avx512 list[i_a]: " << list1[*i_a] << endl;
    		//cout << "avx512 list[i_b]: " << list2[*i_b] << endl;
    		// intersect the tail using scalar intersection
    		//count += u64_intersect_scalar(list1+i_a, size1-i_a, list2+i_b, size2-i_b, result+count);
    		//if(size3 - *i_a - *i_b == 0){
    		//	return count;
    		//}
    		//else{
    			//cout << "ia:" << *i_a << endl;
    			//cout << "ib:" << *i_b << endl;
    			//cout << "count:" << count << endl;
    
                count += u64_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3 - (*i_a+*i_b - count), &i_a_s, &i_b_s);
    
    			*i_a += i_a_s;
    			*i_b += i_b_s;
    		//}
    		return count;
    }
#else
#ifdef __SSE4_1__
    // implement by sse

uint64_t u32_intersection_vector_sse(const uint32_t *list1, uint64_t size1, const uint32_t *list2, uint64_t size2, uint64_t size3, uint64_t *i_a, uint64_t *i_b){
    uint64_t count = 0;
    *i_a = 0;
    *i_b = 0;
    uint64_t st_a = (size1 / 4) * 4;
    uint64_t st_b = (size2 / 4) * 4;

    uint64_t i_a_s, i_b_s;

    if(size3 <= 8){
        count += u32_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
        return count;
    }
    uint64_t stop = size3 - 8;
    while(*i_a < st_a && *i_b < st_b){
        
        uint32_t a_max = list1[*i_a + 3];
        uint32_t b_max = list2[*i_b + 3];
        
        __m128i v_a = _mm_loadu_si128((__m128i*)&(list1[*i_a]));
        __m128i v_b = _mm_loadu_si128((__m128i*)&(list2[*i_b]));

        *i_a += (a_max <= b_max) * 4;
        *i_b += (a_max >= b_max) * 4;

        __m128i cmp1 = _mm_cmpeq_epi32(v_a, v_b);
        __m128i rot1 = _mm_shuffle_epi32(v_b, 57);//00111001
        __m128i cmp2 = _mm_cmpeq_epi32(v_a, rot1);
        __m128i rot2 = _mm_shuffle_epi32(v_b, 78);//01001110
        __m128i cmp3 = _mm_cmpeq_epi32(v_a, rot2);
        __m128i rot3 = _mm_shuffle_epi32(v_b, 147);//10010011
        __m128i cmp4 = _mm_cmpeq_epi32(v_a, rot3);

        __m128i cmp = _mm_or_si128(_mm_or_si128(cmp1, cmp2),_mm_or_si128(cmp3, cmp4));
        int32_t mask = _mm_movemask_ps((__m128)cmp);
        count += _mm_popcnt_u32(mask);

        if(*i_a + *i_b - count >= stop) break;
    }
    count += u32_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3-(*i_a+*i_b-count), &i_a_s, &i_b_s);
    
    *i_a += i_a_s;
    *i_b += i_b_s;

    return count;
}


uint64_t u64_intersection_vector_sse(const uint64_t *list1, uint64_t size1, const uint64_t *list2, uint64_t size2, uint64_t size3, uint64_t *i_a, uint64_t *i_b){
    uint64_t count = 0;
    *i_a = 0;
    *i_b = 0;
    uint64_t st_a = (size1 / 2) * 2;
    uint64_t st_b = (size2 / 2) * 2;

    uint64_t i_a_s, i_b_s;

    if(size3 <= 4){
        count += u64_intersect_scalar_stop(list1, size1, list2, size2, size3, i_a, i_b);
        return count;
    }
    uint64_t stop = size3 - 4;
    while(*i_a < st_a && *i_b < st_b){
        
        uint64_t a_max = list1[*i_a + 1];
        uint64_t b_max = list2[*i_b + 1];
        
        __m128i v_a = _mm_loadu_si128((__m128i*)&(list1[*i_a]));
        __m128i v_b = _mm_loadu_si128((__m128i*)&(list2[*i_b]));

        *i_a += (a_max <= b_max) * 2;
        *i_b += (a_max >= b_max) * 2;

        __m128i cmp1 = _mm_cmpeq_epi64(v_a, v_b);
        //inspect((__m128d)v_b);
        __m128d rot1 = _mm_shuffle_pd((__m128d)v_b, (__m128d)v_b, 1);//00000001
        //inspect(rot1);
        __m128i cmp2 = _mm_cmpeq_epi64(v_a, (__m128i)rot1);

        cmp1 = _mm_or_si128(cmp1, cmp2);
        int64_t mask = _mm_movemask_pd((__m128d)cmp1);
        count += _mm_popcnt_u64(mask);

        if(*i_a + *i_b - count >= stop){ break;}
    }
    count += u64_intersect_scalar_stop(list1+*i_a, size1-*i_a, list2+*i_b, size2-*i_b, size3-(*i_a+*i_b-count), &i_a_s, &i_b_s);
    
    *i_a += i_a_s;
    *i_b += i_b_s;

    return count;
}
#endif
#endif
#endif

    //#endif

} // namespace mash
