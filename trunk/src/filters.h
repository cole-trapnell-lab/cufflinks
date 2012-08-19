#ifndef FILTERS_H
#define FILTERS_H

/*
 *  filters.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 10/27/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "scaffolds.h"
#include "genes.h"

void filter_junk_isoforms(vector<shared_ptr<Abundance> >& transcripts,
						  vector<double>& abundances,
                          const vector<shared_ptr<Abundance> >& mapped_transcripts,
                          double locus_mass);


void filter_introns(int bundle_length,
					int bundle_left,
					vector<Scaffold>& hits, 
					double fraction,
					bool filter_on_intron_overlap,
					bool filter_with_intron_doc);

// Designed to strip out remaining pre-mrna genes, assembled repeats, and 
// fragments from isoforms too short to be reliably quantitated.
void filter_junk_genes(vector<Gene>& genes);


void filter_hits(int bundle_length, int bundle_left, vector<Scaffold>& hits);

void clip_by_3_prime_dropoff(vector<Scaffold>& scaff);

#endif
