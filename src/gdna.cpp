#include "gdna.h"
#include <string.h>

#define IUPAC_DEFS "AaCcTtGgUuMmRrWwSsYyKkVvHhDdBbNnXx-*"
#define IUPAC_COMP "TtGgAaCcAaKkYyWwSsRrMmBbDdHhVvNnXx-*"

unsigned char ntCompTable[256];

static bool gdna_ntCompTableReady=ntCompTableInit();

char ntComplement(char c) {
 return ntCompTable[(int)c];
 }

//in place reverse complement of nucleotide (sub)sequence
char* reverseComplement(char* seq, int slen) {
   if (slen==0) slen=strlen(seq);
   //reverseChars(seq,len);
   int l=0;
   int r=slen-1;
   register char c;
   while (l<r) {
      c=seq[l];seq[l]=seq[r];
      seq[r]=c;   //this was: swap(str[l],str[r]);
      l++;r--;
      }
   for (int i=0;i<slen;i++) seq[i]=ntComplement(seq[i]);
   return seq;
 }

bool ntCompTableInit() {
       //if (gdna_ntCompTableReady) return true;
       char n[]=IUPAC_DEFS;
       char c[]=IUPAC_COMP;
       int l=strlen(IUPAC_DEFS);
       ntCompTable[0]=0;
       for (int ch=1;ch<256;ch++) {
          ntCompTable[ch]=0;
          for (int i=0;i<l;i++)
                if (ch==n[i]) {
                  ntCompTable[ch]=c[i];
                  break;
                  }
          if (ntCompTable[ch]==0)
              ntCompTable[ch]='N';
          }
      //gdna_ntCompTableReady=true;
      return true;
     }




