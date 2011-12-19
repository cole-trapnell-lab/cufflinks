//
//  sampling.cpp
//  cufflinks
//
//  Created by Cole Trapnell on 12/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "sampling.h"
#include <limits>

using namespace std;

void generate_importance_samples(multinormal_generator<double>& generator,
                                 std::vector<boost::numeric::ublas::vector<double> >& samples, 
                                 int num_samples,
                                 bool no_zeros)
{
	for (int i = 0; i < num_samples; ++i)
	{
		boost::numeric::ublas::vector<double> r = generator.next_rand();
        
		boost::numeric::ublas::vector<double> scaled_sample = r;
		
		for (size_t j = 0; j < scaled_sample.size(); ++j) {
            //			if (scaled_sample(j) < 0)
            //				scaled_sample(j) = 1e-10;
            if (scaled_sample(j) < 0)
                scaled_sample(j) = -scaled_sample(j);
		}
		
		double m = sum(scaled_sample);
		if (m && !isnan(m))
		{
			for (size_t j = 0; j < scaled_sample.size(); ++j) 
            {
				scaled_sample(j) = scaled_sample(j) / m;
			}
			if (no_zeros)
            {
                bool has_zero = false;
                for (size_t j = 0; j < scaled_sample.size(); ++j)
                {
                    if (scaled_sample[j] == 0)
                    {
                        has_zero = true;
                        break;
                    }
                }
                
                if (has_zero)
                    continue;
            }
			samples.push_back(scaled_sample);
		}
		else
		{
            samples.push_back(boost::numeric::ublas::zero_vector<double>(scaled_sample.size()));
			//cerr << r << endl;
			//cerr << scaled_sample << endl;
		}
	}
}
