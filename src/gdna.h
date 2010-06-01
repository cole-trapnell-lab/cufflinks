#ifndef GDNA_H
#define GDNA_H
#include "GBase.h"

char ntComplement(char c);

//in-place reverse complement of a nucleotide (sub)sequence
char* reverseComplement(char* seq, int slen=0);

bool ntCompTableInit();

#endif
