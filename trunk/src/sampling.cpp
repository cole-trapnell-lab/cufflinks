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
                                 std::vector<Eigen::VectorXd>& samples, 
                                 int num_samples,
                                 bool no_zeros)
{
	for (int i = 0; i < num_samples; ++i)
	{
        // TODO: we should switch the multinormal generator over to Eigen for
        // consistency as part of the push to drop uBLAS.
		boost::numeric::ublas::vector<double> r = generator.next_rand();
        
        Eigen::VectorXd scaled_sample = Eigen::VectorXd::Zero(r.size());
		
		for (int j = 0; j < scaled_sample.size(); ++j) 
        {
            scaled_sample(j) = r(j);
            
            if (scaled_sample(j) < 0)
            {
                scaled_sample(j) = 0.0;
                //scaled_sample(j) = -scaled_sample(j);
            }
		}
		
		double m = scaled_sample.sum();
		if (m && !isnan(m))
		{
			for (int j = 0; j < scaled_sample.size(); ++j) 
            {
				scaled_sample(j) = scaled_sample(j) / m;
			}
			if (no_zeros)
            {
                bool has_zero = false;
                for (int j = 0; j < scaled_sample.size(); ++j)
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
            samples.push_back(Eigen::VectorXd::Zero(scaled_sample.size()));
			//cerr << r << endl;
			//cerr << scaled_sample << endl;
		}
	}
}
