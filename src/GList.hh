//---------------------------------------------------------------------------
/*
Sortable collection of pointers to objects
*/

#ifndef GListHH
#define GListHH

#include "GBase.h"
//#include "assert.h"

#ifdef __LINE__
#define SLISTINDEX_ERR "GList error (%s:%d):Invalid list index: %d\n"
#define TEST_INDEX(x) \
 if (x<0 || x>=fCount) GError(SLISTINDEX_ERR, __FILE__,__LINE__, x)
#else
#define SLISTINDEX_ERR "GList error:Invalid list index: %d\n"
#define TEST_INDEX(x) \
 if (x<0 || x>=fCount) GError(SLISTINDEX_ERR, x, __FILE__,__LINE__)
#endif

#define SLISTCAPACITY_ERR "GList error: invalid capacity: %d\n"
#define SLISTCOUNT_ERR "GList error: invalid count: %d\n"
#define SLISTSORTED_ERR "Operation not allowed on a sorted list!\n"
#define SLISTUNSORTED_ERR "Operation not allowed on an unsorted list!\n"

// ------ macros:
#define BE_UNSORTED if (fCompareProc!=NULL) { GError(SLISTSORTED_ERR); return; }
#define BE_SORTED if (fCompareProc==NULL) { GError(SLISTUNSORTED_ERR); return; }

#define MAXLISTSIZE INT_MAX-1

#define SORTED (fCompareProc!=NULL)
#define UNSORTED (fCompareProc==NULL)
#define FREEDATA (fFreeProc!=NULL)
/* #define TEST_INDEX(x) assert(x>=0 && x<fCount); \
     if (x<0 || x>=fCount) GError(SLISTINDEX_ERR, x) */


//template for array of objects; GVec is not sortable, 
// so it doesn't require comparison operators defined
template <class OBJ> class GVec {
  protected:
    OBJ* fArray;
    int fCount;
    int fCapacity;
  public:
    GVec(int init_capacity=20);
    GVec(GVec<OBJ>& array); //copy constructor
    const GVec<OBJ>& operator=(GVec& array); //copy operator
    virtual ~GVec();
    void idxInsert(int idx, OBJ& item);
    void Grow();
    void Grow(int idx, OBJ& item);
    void Reverse(); //WARNING: will break the sort order if SORTED!
    int Add(OBJ* item); // simply append to the end of fArray, reallocating as needed
    int Add(OBJ& item) { return Add(&item); } //both will CREATE a new OBJ and COPY to it
                       // using OBJ new operator=
    void Add(GVec<OBJ>& list); //append copies of all items from another list
    OBJ& Get(int idx) {
          TEST_INDEX(idx);
          return fArray[idx];
          }
    OBJ& operator[](int i) {
          TEST_INDEX(i);
          return fArray[i];
          }
    void Clear();
    void Insert(int idx, OBJ* item);
    void Delete(int index);
    void Replace(int idx, OBJ& item); //Put, use operator= to copy
    void Exchange(int idx1, int idx2);
    void Swap(int idx1, int idx2)  { Exchange(idx1, idx2); } 
    int  Capacity() { return fCapacity; }
    //this will reject identical items in sorted lists only!
    void setCapacity(int NewCapacity);
    int  Count() { return fCount; }
    void setCount(int NewCount); //?! - will trim or expand the array as needed
    void Move(int curidx, int newidx);
};

// GArray is the sortable collection, but requires the comparison operators to be defined
template <class OBJ> class GArray:public GVec<OBJ> {
  protected:
    bool fUnique;
    static int DefaultCompareProc(OBJ& item1, OBJ& item2) {
      //the comparison operators MUST be defined for OBJ class!
      if ( item1 > item2) return 1;
        else return (item2 > item1) ? -1 : 0 ;
      }
  public:
    typedef int CompareProc(OBJ& item1, OBJ& item2);
  protected:
    CompareProc* fCompareProc;
    void qSort(int L, int R);
  public:
    GArray(CompareProc* cmpFunc=NULL);
    GArray(bool sorted, bool unique=false);
    GArray(int init_capacity, bool sorted, bool unique=false);
    GArray(GArray<OBJ>& array); //copy constructor
    const GArray<OBJ>& operator=(GArray& array);
    //~GArray();
    //assignment operator
    void setSorted(CompareProc* cmpFunc);
    //sort the array if cmpFunc not NULL or changes
    int Add(OBJ* item); // specific implementation if sorted
    int Add(OBJ& item) { return Add(&item); } //both will CREATE a new OBJ and COPY to it
                       // using OBJ new operator=
    void Add(GArray<OBJ>& list); //add copies of all items from another list
    //this will reject identical items in sorted lists only!
    void setUnique(bool beUnique) { fUnique = beUnique; };
    void Sort(); //explicit sort may be requested
    bool Sorted() { return fCompareProc!=NULL; }
    void Replace(int idx, OBJ& item); //Put, use operator= to copy
    int  Unique() { return fUnique; }
    int IndexOf(OBJ& item);
         //this needs the == operator to have been defined for OBJ
    bool Found(OBJ& item, int& idx); // for sorted arrays only;
         //search by content; if found, returns true and idx will be the index
         //of the first item found matching for which CompareProc returns 0
    bool Exists(OBJ& item); //same as above without existing index info
    //unsorted only, place item at position idx:
    void Move(int curidx, int newidx);
    void Insert(int idx, OBJ* item);
    void Insert(int idx, OBJ& item) { Insert(idx,&item); }
};

//------- template for array of pointers to objects ---------
template <class OBJ> class GPVec {
  protected:
    OBJ** fList; //pointer to an array of pointers to objects
    int fCount; //total number of entries in list
    int fCapacity; //current allocated size
    GFreeProc* fFreeProc; //useful for deleting objects
    //---
    void Expand();
    void Grow();
    void Grow(int idx, OBJ* newitem);
    void setCount(int NewCount); //will trim/expand the array as needed
  public:  
    static void DefaultFreeProc(pointer item) {
      delete (OBJ*)item;
      }
    virtual ~GPVec();
    GPVec(int init_capacity=10, bool free_elements=true); //also the default constructor
    GPVec(GPVec<OBJ>& list); //copy constructor?
    GPVec(GPVec<OBJ>* list); //kind of a copy constructor
    const GPVec<OBJ>& operator=(GPVec& list);
    OBJ* Get(int i);
    OBJ* operator[](int i) { return this->Get(i); }
    void Reverse(); //reverse pointer array; WARNING: will break the sort order if sorted!
    void freeItem(int idx); //calls fFreeProc (or DefaultFreeProc) on fList[idx] and sets NULL there, doesn't pack!
                      //it will free even if fFreeProc is NULL!
    void setFreeItem(GFreeProc *freeProc) { fFreeProc=freeProc; }
    void setFreeItem(bool doFree) {
       if (doFree) fFreeProc=DefaultFreeProc;
             else  fFreeProc=NULL;
       }
    // -- stack usage:
    int Push(OBJ* item) { return Add(item); }
    OBJ* Pop();// Stack use; removes and returns last item,but does NOT FREE it
    OBJ* Shift(); //Queue use: removes and returns first item, but does NOT FREE it
    void deallocate_item(OBJ* item); //forcefully call fFreeProc or delete on item
    void Clear();
    void Exchange(int idx1, int idx2);
    void Swap(int idx1, int idx2)  { Exchange(idx1, idx2); }
    OBJ* First() { return (fCount>0)?fList[0]:NULL; }
    OBJ* Last()  { return (fCount>0)?fList[fCount-1]:NULL;}
    bool isEmpty() { return fCount==0; }
    bool notEmpty() { return fCount>0; }
    int Capacity() { return fCapacity; }
    int Count()   { return fCount; }
    void setCapacity(int NewCapacity);
    int Add(OBJ* item); //simply append the pointer copy
    void Add(GPVec<OBJ>& list); //add all pointers from another list
    void Insert(int idx, OBJ* item);
    void Move(int curidx, int newidx);
    void Put(int idx, OBJ* item);
    void Pack();
    void Delete(int index); //also frees the item if fFreeProc!=NULL, and shifts the successor items
    void Forget(int idx); //simply places a NULL at fList[idx], nothing else
    int RemovePtr(pointer item); //always use linear search to find the pointer! calls Delete() if found
    int IndexOf(pointer item); //a linear search for pointer address only!
 };
 
template <class OBJ> class GList:public GPVec<OBJ> {
  protected:
    bool fUnique;
    GCompareProc* fCompareProc; //a pointer to a Compare function
    static int DefaultCompareProc(const pointer item1, const pointer item2) {
      //the comparison operators MUST be defined for OBJ class!
      if (*((OBJ*)item1) > *((OBJ*)item2)) return 1;
        else if (*((OBJ*)item2) > *((OBJ*)item1)) return -1;
                                             else return  0;
      }
    void QuickSort(int L, int R);
  public:
    void sortInsert(int idx, OBJ* item);
    GList(GCompareProc* compareProc=NULL); //free by default
    GList(GCompareProc* compareProc, //unsorted by default
        GFreeProc *freeProc,
        bool beUnique=false);
    GList(bool sorted, bool free_elements=true, bool beUnique=false);
    GList(int init_capacity, bool sorted, bool free_elements=true, bool beUnique=false);
    GList(GList<OBJ>& list); //copy constructor?
    GList(GList<OBJ>* list); //kind of a copy constructor
    const GList<OBJ>& operator=(GList& list);
    //void Clear();
    //~GList();
    void setSorted(GCompareProc* compareProc);
       //sorted if compareProc not NULL; sort the list if compareProc changes !
    bool Sorted() { return fCompareProc!=NULL; }
    void setSorted(bool sorted) {
     if (sorted) {
         if (fCompareProc!=&DefaultCompareProc) {
             fCompareProc=&DefaultCompareProc;
             Sort();
             }
          }
      else fCompareProc=NULL;
      }
    int Add(OBJ* item); //-- specific implementation if sorted
    void Add(GList<OBJ>& list); //add all pointers from another list

    OBJ* AddIfNew(OBJ* item, bool deleteIfFound=true, int* fidx=NULL);
    // default: delete item if Found() (and pointers are not equal)!
    //returns the equal (==) object if it's in the list already
    //or the item itself if it is unique and actually added

    int AddedIfNew(OBJ* item);
    // if Found(item) (and pointers are not equal) delete item and returns -1
    // if added, returns the new item index


    int Unique() { return fUnique; }
    //this will reject identical items in sorted lists only!
    void setUnique(bool beUnique) { fUnique = beUnique; };

    GCompareProc* GetCompareProc() {return fCompareProc;}
    int IndexOf(OBJ* item); //this has a specific implementation for sorted lists
               //if list is sorted, item data is located by binary search
               //based on the Compare function
               //if not, a linear search is performed, but
               //this needs the == operator to have been defined for OBJ
    
    void Put(int idx, OBJ* item, bool re_sort=false);
    bool Found(OBJ* item, int & idx); // sorted only;
               //search by content; if found, returns true and idx will be the index
               //of the first item found matching for which GTCompareProc returns 0
    bool Exists(OBJ* item); //same as above without existing index info
    bool Exists(OBJ& item); //same as above without existing index info
    void Sort(); //explicit sort may be requested using this function
    int Remove(OBJ* item); //search for pointer, using binary search if sorted
    void Insert(int idx, OBJ* item); //unsorted only, place item at position idx
    void Move(int curidx, int newidx);
}; //GList 


//basic template for a Stack of pointers (implemented as a linked list)
template <class OBJ> class GStack {
 protected:
   struct StackOBJ {
      OBJ* obj;
      StackOBJ* prev;
      };
   int fCount; //total number of elements in stack
   StackOBJ* base;
   StackOBJ* top;
 public:
   GStack(OBJ* po=NULL) {
      base=NULL;
      top=NULL;
      fCount=0;
      if (po!=NULL) Push(po);
      }
   ~GStack() {
      while (fCount>0) Pop();
      }
   bool isEmpty() { return fCount==0; }
   int Size() { return fCount; }
   int Count() { return fCount; }
   OBJ* Pop() {
      if (top==NULL) return NULL;
      fCount--;
      StackOBJ* ctop=top;
      if (top==base) base=NULL;
      OBJ* r=top->obj;
      top=top->prev;
      GFREE(ctop);
      return r;
      }
   OBJ* Push(OBJ* o) {
      fCount++;
      StackOBJ* ctop=top; //could be NULL
      GMALLOC(top, sizeof(StackOBJ));
      top->obj=o;
      top->prev=ctop;
      if (base==NULL) base=top;
      return o;
      }
   OBJ* Top() { return ((top==NULL)? NULL : top->obj); }
   OBJ* Base() { return ((base==NULL)? NULL : base->obj); }
};

//-------------------- TEMPLATE IMPLEMENTATION-------------------------------

template <class OBJ> GVec<OBJ>::GVec(int init_capacity) {
  fCount=0;
  fCapacity=0;
  fArray=NULL;
  setCapacity(init_capacity);
}

template <class OBJ> GVec<OBJ>::GVec(GVec<OBJ>& array) { //copy constructor
 this->fCount=array.fCount;
 this->fCapacity=array.fCapacity;
 this->fArray=NULL;
 if (this->fCapacity>0) {
    //GMALLOC(fArray, fCapacity*sizeof(OBJ));
    fArray=new OBJ[this->fCapacity];
    }
 this->fCount=array.fCount;
 // uses OBJ operator=
 for (int i=0;i<this->fCount;i++) fArray[i]=array[i];
 }

template <class OBJ> GArray<OBJ>::GArray(GArray<OBJ>& array):GVec<OBJ>(0) { //copy constructor
 this->fCount=array.fCount;
 this->fCapacity=array.fCapacity;
 this->fArray=NULL;
 if (this->fCapacity>0) {
    //GMALLOC(this->fArray, this->fCapacity*sizeof(OBJ));
    this->fArray=new OBJ[this->fCapacity];
    }
 this->fCount=array.fCount;
 fUnique=array.fUnique;
 fCompareProc=array.fCompareProc;
 // uses OBJ operator=
 for (int i=0;i<this->fCount;i++) this->fArray[i]=array[i];
 }

template <class OBJ> const GVec<OBJ>& GVec<OBJ>::operator=(GVec<OBJ>& array) {
 if (&array==this) return *this;
 Clear();
 fCount=array.fCount;
 fCapacity=array.fCapacity;
 if (fCapacity>0) {
    //GMALLOC(fArray, fCapacity*sizeof(OBJ));
    fArray=new OBJ[this->fCapacity];
    }
 fCount=array.fCount;
 // uses OBJ operator=
 for (int i=0;i<fCount;i++) {
   fArray[i]=array[i];
   }
 return *this;
}

template <class OBJ> const GArray<OBJ>& GArray<OBJ>::operator=(GArray<OBJ>& array) {
 if (&array==this) return *this;
 GVec<OBJ>::Clear();
 this->fCount=array.fCount;
 this->fUnique=array.fUnique;
 this->fCapacity=array.fCapacity;
 if (this->fCapacity>0) {
    //GMALLOC(this->fArray, this->fCapacity*sizeof(OBJ));
    this->fArray=new OBJ[this->fCapacity];
    }
 this->fCompareProc=array.fCompareProc;
 this->fCount=array.fCount;
 // uses OBJ operator=
 for (int i=0;i<this->fCount;i++) {
   this->fArray[i]=array[i];
   }
 return *this;
}

template <class OBJ> GArray<OBJ>::GArray(CompareProc* cmpFunc):GVec<OBJ>(0) {
  fCompareProc = cmpFunc;
  fUnique = false; //only affects sorted lists
}

template <class OBJ> GArray<OBJ>::GArray(bool sorted, bool unique):GVec<OBJ>(0) {
  fUnique=unique;
  fCompareProc=sorted? &DefaultCompareProc : NULL;
}

template <class OBJ> GArray<OBJ>::GArray(int init_capacity,
                        bool sorted, bool unique):GVec<OBJ>(init_capacity) {
  fUnique=unique;
  fCompareProc=sorted? &DefaultCompareProc : NULL;
}

template <class OBJ> GVec<OBJ>::~GVec() {
 this->Clear();
}


template <class OBJ> void GVec<OBJ>::setCapacity(int NewCapacity) {
  if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
    GError(SLISTCAPACITY_ERR, NewCapacity);
    //error: capacity not within range
  if (NewCapacity!=fCapacity) {
   if (NewCapacity==0) {
      delete[] fArray;
      }
    else {
      //GREALLOC(fArray, NewCapacity*sizeof(OBJ));
      OBJ* oldArray=fArray;
      fArray=new OBJ[NewCapacity];
      for (int i=0;i<this->fCount;i++) {
        fArray[i] = oldArray[i];
        }
      delete[] oldArray;
      }
   fCapacity=NewCapacity;
   }
}

template <class OBJ> void GVec<OBJ>::Clear() {
  setCount(0);
  setCapacity(0); //so the array itself is deallocated too!
}

/*
template <class OBJ> void GArray<OBJ>::Clear() {
  CompareProc* fcmp=fCompareProc;
  fCompareProc=NULL;
  GVec<OBJ>::setCount(0);
  GVec<OBJ>::setCapacity(0); //so the array itself is deallocated too!
  fCompareProc=fcmp;
}
*/
template <class OBJ> void GArray<OBJ>::setSorted(CompareProc* cmpFunc) {
  CompareProc* old_proc=fCompareProc;
  fCompareProc=cmpFunc;
  if (fCompareProc!=old_proc && fCompareProc!=NULL)
       Sort(); //new compare method
}

template <class OBJ> void GVec<OBJ>::Grow() {
 int delta;
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
  setCapacity(fCapacity + delta);
}

template <class OBJ> void GVec<OBJ>::Reverse() {
  int l=0;
  int r=fCount-1;
  OBJ c;
  while (l<r) {
     c=fArray[l];fArray[l]=fArray[r];
     fArray[r]=c;
     l++;r--;
     }
}

template <class OBJ> void GVec<OBJ>::Grow(int idx, OBJ& item) {
 int delta;
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
 int NewCapacity=fCapacity+delta;
  if (NewCapacity <= fCount || NewCapacity >= MAXLISTSIZE)
    GError(SLISTCAPACITY_ERR, NewCapacity);
    //error: capacity not within range

  if (NewCapacity!=fCapacity) {
    if (NewCapacity==0) {
      //GFREE(fArray);
      delete[] fArray;
      fArray=NULL;
      }
    else { //add the new item
      if (idx==fCount) { //append item
         //GREALLOC(fArray, NewCapacity*sizeof(OBJ));
         setCapacity(NewCapacity);
         fArray[idx]=item;
         }
       else { //insert item at idx
        OBJ* newList;
        //GMALLOC(newList, NewCapacity*sizeof(OBJ));
        newList=new OBJ[NewCapacity];
        //copy data before idx
        //memmove(&newList[0],&fArray[0], idx*sizeof(OBJ));
        // operator= required!
        for (int i=0;i<idx;i++) {
          newList[i]=fArray[i];
          }
        newList[idx]=item;
        //copy data after idx
        //memmove(&newList[idx+1],&fArray[idx], (fCount-idx)*sizeof(OBJ));
        for (int i=idx+1;i<=fCount;i++) {
          newList[i]=fArray[i-1];
          }
        //memset(&newList[fCount+1], 0, (NewCapacity-fCount-1)*sizeof(OBJ));
        //data copied:
        //GFREE(fArray);
        delete[] fArray;
        fArray=newList;
        }
      fCount++;
      }
   fCapacity=NewCapacity;
   }
}

template <class OBJ> int GArray<OBJ>::IndexOf(OBJ& item) {
 int result=0;
 if (Found(item, result)) return result;
                     else return -1;
 }

template <class OBJ> bool GArray<OBJ>::Exists(OBJ& item) {
 int result=0;
 if (Found(item, result)) return true;
                     else return false;
 }


template <class OBJ> int GVec<OBJ>::Add(OBJ* item) {
 if (item==NULL) return -1;
 int result=fCount;
 if (result==fCapacity) Grow();
 fArray[result] = *item; //OBJ::operator= must copy OBJ properly!
 fCount++;
 return result;
}

template <class OBJ> int GArray<OBJ>::Add(OBJ* item) {
 if (item==NULL) return -1;
 int result;
 if (SORTED) {
   if (Found(*item, result))
      if (fUnique) return -1; //cannot add a duplicate!
   //Found sets result to the position where the item should be!
   this->idxInsert(result, *item);
   }
  else {
   if (fUnique && Found(*item,result)) return -1; //set behaviour
   result = this->fCount;
   if (result==this->fCapacity) GVec<OBJ>::Grow();
   this->fArray[result] = *item; //operator=, copies the item
   this->fCount++;
   }
 return result;
}


template <class OBJ> void GVec<OBJ>::Add(GVec<OBJ>& list) {
  if (list.Count()==0) return;
  //simply copy
  setCapacity(fCapacity+list.fCount);
  int s=fCount;
  for (int i=0;i<list.fCount;i++)
           fArray[s+i]=list.fArray[i];
  fCount+=list.fCount;
}


template <class OBJ> void GArray<OBJ>::Add(GArray<OBJ>& list) {
  if (list.Count()==0) return;
  if (SORTED) {
    for (int i=0;i<list.fCount;i++) Add(&list[i]);
    }
  else { //simply copy
    this->setCapacity(this->fCapacity+list.fCount);
    int s=this->fCount;
    for (int i=0;i<list.fCount;i++)
           this->fArray[s+i]=list.fArray[i];
    this->fCount+=list.fCount;
    }
}

template <class OBJ> bool GArray<OBJ>::Found(OBJ& item, int& idx) {
 //search the list by using CompareProc (if defined)
 //or == operator for a non-sortable list
 //for sorted lists, even when the result is false, the idx is
 //set to the closest matching object!
 int i;
 idx=-1;
 if (this->fCount==0) { idx=0;return false;}
 if (SORTED) { //binary search based on CompareProc
   //do the simplest tests first:
   if ((*fCompareProc)(this->fArray[0],item)>0) {
                       idx=0;
                       return false;
                       }
   if ((*fCompareProc)(item, this->fArray[this->fCount-1])>0) {
                       idx=this->fCount;
                       return false;
                       }

   int l=0;
   int h = this->fCount - 1;
   int c;
   while (l <= h) {
       i = (l + h) >> 1;
       c = (*fCompareProc)(this->fArray[i], item);
       if (c < 0)  l = i + 1;
         else {
            h = i - 1;
            if (c == 0) { //found!
                 idx=i;
                 return true;
                }
            }
       } //while
   idx = l;
   return false;
   }
 else {//not sorted: use linear search
   // needs == operator to compare user defined objects !
   i=0;
   while (i<this->fCount) {
      if (this->fArray[i]==item) { //requires operator==
         idx=i;
         return true;
         }
      i++;
      }
   return false;
   }
}

template <class OBJ> void GVec<OBJ>::idxInsert(int idx, OBJ& item) {
 //idx must be the new position this new item must have
 //so the allowed range is [0..fCount]
 //the old idx item all the above will be shifted to idx+1
 if (idx<0 || idx>fCount) GError(SLISTINDEX_ERR, idx);
 if (fCount==fCapacity) { //need to resize the array
    Grow(idx, item); //expand and also copy/move data and insert the new item
    return;
    }
 //move data around to make room for the new item
 if (idx<fCount) {
      //copy after-idx items (shift up) 
      //memmove(&newList[idx+1],&fArray[idx], (fCount-idx)*sizeof(OBJ));
      for (int i=fCount; i>idx; i--) {
          fArray[i]=fArray[i-1];
          }
      }
 fArray[idx]=item;
 fCount++;
}

template <class OBJ> void GVec<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 idxInsert(idx, item);
}

template <class OBJ> void GArray<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 BE_UNSORTED; //forbid this operation on sorted data
 idxInsert(idx, item);
}


template <class OBJ> void GVec<OBJ>::Move(int curidx, int newidx) {
 if (curidx!=newidx || newidx>=fCount)
     GError(SLISTINDEX_ERR, newidx);
 OBJ tmp=fArray[curidx]; //copy constructor here
 fArray[curidx]=fArray[newidx];
 fArray[newidx]=tmp;
}


template <class OBJ> void GArray<OBJ>::Move(int curidx, int newidx) {
 BE_UNSORTED; //cannot do this in a sorted list!
 if (curidx!=newidx || newidx>=this->fCount)
     GError(SLISTINDEX_ERR, newidx);

 OBJ tmp=this->fArray[curidx]; //copy constructor here
 this->fArray[curidx]=this->fArray[newidx];
 this->fArray[newidx]=tmp;
}

template <class OBJ> void GVec<OBJ>::Replace(int idx, OBJ& item) {
 TEST_INDEX(idx);
 fArray[idx]=item;
}

template <class OBJ> void GVec<OBJ>::Exchange(int idx1, int idx2) {
 TEST_INDEX(idx1);
 TEST_INDEX(idx2);
 OBJ item=fArray[idx1];
 fArray[idx1]=fArray[idx2];
 fArray[idx2]=item;
}


template <class OBJ> void GArray<OBJ>::Replace(int idx, OBJ& item) {
 //TEST_INDEX(idx);
 if (idx<0 || idx>=this->fCount) GError(SLISTINDEX_ERR, __FILE__,__LINE__, idx);
 this->fArray[idx]=item;
 if ( SORTED ) Sort(); //re-sort ! this could be very expensive, don't do it
}


template <class OBJ> void GVec<OBJ>::Delete(int index) {
 TEST_INDEX(index);
 fCount--;
 while (index<fCount) {
    //move higher elements if any (shift down)
    //memmove(&fArray[index], &fArray[index+1], (fCount-index)*sizeof(OBJ));
    fArray[index]=fArray[index+1];
    index++;
    }
}

template <class OBJ> void GVec<OBJ>::setCount(int NewCount) {
  if (NewCount<0 || NewCount > MAXLISTSIZE)
     GError(SLISTCOUNT_ERR, NewCount);
  if (NewCount > fCapacity) setCapacity(NewCount);
  //if (NewCount > fCount)
  //  memset(&fArray[fCount], 0, (NewCount - fCount) * sizeof(OBJ));
  fCount = NewCount;
}

template <class OBJ> void GArray<OBJ>::qSort(int l, int r) {
 int i, j;
 OBJ p,t;
 do {
    i = l; j = r;
    p = this->fArray[(l + r) >> 1];
    do {
      while (fCompareProc(this->fArray[i], p) < 0) i++;
      while (fCompareProc(this->fArray[j], p) > 0) j--;
      if (i <= j) {
        t = this->fArray[i];
        this->fArray[i] = this->fArray[j];
        this->fArray[j] = t;
        i++; j--;
        }
      } while (i <= j);
    if (l < j) qSort(l, j);
    l = i;
    } while (i < r);
}

template <class OBJ> void GArray<OBJ>::Sort() {
 if (this->fArray!=NULL && this->fCount>0 && fCompareProc!=NULL)
     qSort(0, this->fCount-1);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*=> GPVec and GList implementation -- sortable array of pointers to OBJ

template <class OBJ> GPVec<OBJ>::GPVec(GPVec& list) { //copy constructor
 fCount=list.fCount;
 fCapacity=list.fCapacity;
 if (fCapacity>0) {
      GMALLOC(fList, fCapacity*sizeof(OBJ*));
      }
 fFreeProc=list.fFreeProc;
 fCount=list.fCount;
 memcpy(fList, list.fList, fCount*sizeof(OBJ*));
 //for (int i=0;i<list.Count();i++) Add(list[i]);
}

template <class OBJ> GPVec<OBJ>::GPVec(GPVec* plist) { //another copy constructor
 fCount=0;
 fCapacity=plist->fCapacity;
 fList=NULL;
 if (fCapacity>0) {
     GMALLOC(fList, fCapacity*sizeof(OBJ*));
     }
 fFreeProc=plist->fFreeProc;
 fCount=plist->fCount;
 memcpy(fList, plist->fList, fCount*sizeof(OBJ*));
 //for (int i=0;i<list->fCount;i++) Add(plist->Get(i));
}

template <class OBJ> const GPVec<OBJ>& GPVec<OBJ>::operator=(GPVec& list) {
 if (&list!=this) {
     Clear();
     fFreeProc=list.fFreeProc;
     //Attention: the object *POINTERS* are copied,
     // but the actual object content is NOT duplicated
     for (int i=0;i<list.Count();i++) Add(list[i]);
     }
 return *this;
}


template <class OBJ> void GPVec<OBJ>::Add(GPVec<OBJ>& list) {
  if (list.Count()==0) return;
  //simply copy the pointers! -- the objects will be shared
  setCapacity(fCapacity+list.fCount);
  memcpy( & (fList[fCount]), list.fList, list.fCount*sizeof(OBJ*));
  fCount+=list.fCount;
}


template <class OBJ> GList<OBJ>::GList(GList<OBJ>& list):GPVec<OBJ>(list) { //copy constructor
 fUnique=list.fUnique;
 fCompareProc=list.fCompareProc;
}

template <class OBJ> GList<OBJ>::GList(GList<OBJ>* plist):GPVec<OBJ>(0) { //another copy constructor
 this->fCapacity=plist->fCapacity;
 this->fList=NULL;
 if (this->fCapacity>0) {
     GMALLOC(this->fList, this->fCapacity*sizeof(OBJ*));
     }
 fUnique=plist->fUnique;
 fCompareProc=plist->fCompareProc;
 this->fFreeProc=plist->fFreeProc;
 this->fCount=plist->fCount;
 memcpy(this->fList, plist->fList, this->fCount*sizeof(OBJ*));
 //for (int i=0;i<list->fCount;i++) Add(plist->Get(i));
}

template <class OBJ> void GList<OBJ>::Add(GList<OBJ>& list) {
  if (list.Count()==0) return;
  if (SORTED) {
    for (int i=0;i<list.Count();i++) Add(list[i]);
    }
  else { //simply copy
    this->setCapacity(this->fCapacity+list.fCount);
    memcpy( & (this->fList[this->fCount]), list.fList, list.fCount*sizeof(OBJ*));
    this->fCount+=list.fCount;
    }
}


template <class OBJ> GList<OBJ>::GList(GCompareProc* compareProc,
       GFreeProc* freeProc, bool beUnique) {
  fCompareProc = compareProc;
  this->fFreeProc    = freeProc;
  fUnique = beUnique; //only affects sorted lists
}

template <class OBJ> GList<OBJ>::GList(GCompareProc* compareProc) {
  fCompareProc = compareProc;
  this->fFreeProc = GPVec<OBJ>::DefaultFreeProc;
  fUnique = false; //only affects sorted lists
}

template <class OBJ> void GPVec<OBJ>::Reverse() {
  int l=0;
  int r=fCount-1;
  OBJ* c;
  while (l<r) {
     c=fList[l];fList[l]=fList[r];
     fList[r]=c;
     l++;r--;
     }
}

template <class OBJ> GList<OBJ>::GList(bool sorted,
    bool free_elements, bool beUnique) {
  if (sorted) {
     if (free_elements) {
        fCompareProc=&DefaultCompareProc;
        this->fFreeProc = GPVec<OBJ>::DefaultFreeProc;
        fUnique=beUnique;
        }
       else {
        fCompareProc=&DefaultCompareProc;
        this->fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
   else {
     if (free_elements) {
        fCompareProc=NULL;
        this->fFreeProc=GPVec<OBJ>::DefaultFreeProc;
        fUnique=beUnique;
        }
      else {
        fCompareProc=NULL;
        this->fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
}

template <class OBJ> GPVec<OBJ>::GPVec(int init_capacity, bool free_elements) {
  fCount=0;
  fCapacity=0;
  fList=NULL;
  fFreeProc=(free_elements) ? DefaultFreeProc : NULL;
  setCapacity(init_capacity);
}


template <class OBJ> GList<OBJ>::GList(int init_capacity, bool sorted,
    bool free_elements, bool beUnique):GPVec<OBJ>(init_capacity, free_elements) {
  if (sorted) {
      fCompareProc=&DefaultCompareProc;
      fUnique=beUnique;
      }
   else {
      fCompareProc=NULL;
      fUnique=beUnique;
      }
}

template <class OBJ> GPVec<OBJ>::~GPVec() {
 this->Clear();//this will free the items if fFreeProc is defined
}

/*
template <class OBJ> GList<OBJ>::~GList() {
 this->Clear();//this will free the items if fFreeProc is defined
}
*/

template <class OBJ> void GPVec<OBJ>::setCapacity(int NewCapacity) {
  if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
    GError(SLISTCAPACITY_ERR, NewCapacity);
    //error: capacity not within range
  if (NewCapacity!=fCapacity) {
   if (NewCapacity==0) {
      GFREE(fList);
      }
    else {
      GREALLOC(fList, NewCapacity*sizeof(OBJ*));
      }
   fCapacity=NewCapacity;
   }
}

template <class OBJ> void GPVec<OBJ>::deallocate_item(OBJ* item) {
 if (item==NULL) return;
 if (FREEDATA) {
   (*fFreeProc)(item);
   }
 else {
  delete item;
  }
}

template <class OBJ> void GPVec<OBJ>::Clear() {
 if (FREEDATA) {
   for (int i=0; i<fCount; i++) {
     (*fFreeProc)(fList[i]);
     }
   }
 setCount(0);
 setCapacity(0); //so the array itself is deallocated too!
}

template <class OBJ> void GPVec<OBJ>::Exchange(int idx1, int idx2) {
//Warning: this will BREAK sort order for sorted GList
 TEST_INDEX(idx1);
 TEST_INDEX(idx2);
 OBJ* item=fList[idx1];
 fList[idx1]=fList[idx2];
 fList[idx2]=item;
}

template <class OBJ> void GPVec<OBJ>::Expand() {
 if (fCount==fCapacity) Grow();
 //return this;
}

template <class OBJ> OBJ* GPVec<OBJ>::Get(int idx) {
 TEST_INDEX(idx);
 return fList[idx];
}

template <class OBJ> const GList<OBJ>& GList<OBJ>::operator=(GList& list) {
 if (&list!=this) {
     GPVec<OBJ>::Clear();
     fCompareProc=list.fCompareProc;
     this->fFreeProc=list.fFreeProc;
     //Attention: the object pointers are copied directly,
     //but the actual objects are NOT duplicated
     for (int i=0;i<list.Count();i++) Add(list[i]);
     }
 return *this;
}

template <class OBJ> void GList<OBJ>::setSorted(GCompareProc* compareProc) {
 GCompareProc* old_proc=fCompareProc;
 fCompareProc=compareProc;
 if (fCompareProc!=old_proc && fCompareProc!=NULL)
       Sort(); //new compare method
}

template <class OBJ> void GPVec<OBJ>::Grow() {
 int delta;
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
  setCapacity(fCapacity + delta);
}

template <class OBJ> void GPVec<OBJ>::Grow(int idx, OBJ* newitem) {
 int delta;
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
 // setCapacity(fCapacity + delta);
 int NewCapacity=fCapacity+delta;
  if (NewCapacity <= fCount || NewCapacity > MAXLISTSIZE)
    GError(SLISTCAPACITY_ERR, NewCapacity);
    //error: capacity not within range
  if (NewCapacity!=fCapacity) {
    if (NewCapacity==0) {
      GFREE(fList);
      }
    else  {//add the new item
      if (idx==fCount) {
         GREALLOC(fList, NewCapacity*sizeof(OBJ*));
         fList[idx]=newitem;
         }
       else {
        OBJ** newList;
        GMALLOC(newList, NewCapacity*sizeof(OBJ*));
        //copy data before idx
        memmove(&newList[0],&fList[0], idx*sizeof(OBJ*));
        newList[idx]=newitem;
        //copy data after idx
        memmove(&newList[idx+1],&fList[idx], (fCount-idx)*sizeof(OBJ*));
        memset(&newList[fCount+1], 0, (NewCapacity-fCount-1)*sizeof(OBJ*));
        //data copied:
        GFREE(fList);
        fList=newList;
        }
      fCount++;
      }
   fCapacity=NewCapacity;
   }
}

template <class OBJ> int GPVec<OBJ>::IndexOf(pointer item) {
 int result=-1;
 for (int i=0;i<fCount;i++) {
     if (item==(pointer)fList[i]) return i;
     }
 return -1;
 }

template <class OBJ> int GList<OBJ>::IndexOf(OBJ* item) {
 int result=0;
 if (Found(item, result)) return result;
                     else return -1;
 }

template <class OBJ> bool GList<OBJ>::Exists(OBJ& item) {
 int result=0;
 if (Found(&item, result)) return true;
                      else return false;
 }

template <class OBJ> bool GList<OBJ>::Exists(OBJ* item) {
 int result=0;
 if (Found(item, result)) return true;
                      else return false;
 }

template <class OBJ> int GPVec<OBJ>::Add(OBJ* item) {
 int result;
 if (item==NULL) return -1;
 result = fCount;
 if (result==fCapacity) this->Grow();
 fList[result]=item;
 fCount++;
 return fCount-1;
}

template <class OBJ> int GList<OBJ>::Add(OBJ* item) {
 int result;
 if (item==NULL) return -1;
 if (SORTED) {
   if (Found(item, result))
      if (fUnique) return -1; //duplicates forbidden
   //Found sets result to the position where the item should be!
   sortInsert(result, item);
   }
  else {
   if (fUnique && Found(item,result)) return -1; //set behaviour
   result = this->fCount;
   if (result==this->fCapacity) GPVec<OBJ>::Grow();
   this->fList[result]=item;
   this->fCount++;
   }
 return result;
}

//by default, it deletes the item if it has an equal in the list!
//returns the existing equal (==) object if it's in the list already
//or returns the item itself if it's unique (and adds it)
template <class OBJ> OBJ* GList<OBJ>::AddIfNew(OBJ* item,
                                     bool deleteIfFound, int* fidx) {
 int r;
 if (Found(item, r)) {
    if (deleteIfFound && (pointer)item != (pointer)(this->fList[r])) {
       this->deallocate_item(item);
       }
    if (fidx!=NULL) *fidx=r;
    return this->fList[r]; //found
    }
 //not found:
 if (SORTED) {
   //Found() set result to the position where the item should be inserted:
   sortInsert(r, item);
   }
  else {
   r = this->fCount;
   if (r==this->fCapacity) GPVec<OBJ>::Grow();
   this->fList[r]=item;
   this->fCount++;
   }
 if (fidx!=NULL) *fidx=r;
 return item;
}

//if item is found already in the list DELETE it and return -1
//otherwise the item is added and its index is returned
template <class OBJ> int GList<OBJ>::AddedIfNew(OBJ* item) {
 int r;
 if (Found(item, r)) {
    if ((pointer)item != (pointer)(this->fList[r])) {
        this->deallocate_item(item);
        }
    return -1;
    }
 //not found:
 if (SORTED) {
   //Found() set r to the position where the item should be inserted:
   sortInsert(r, item);
   }
  else {
   r = this->fCount;
   if (r==this->fCapacity) GPVec<OBJ>::Grow();
   this->fList[r]=item;
   this->fCount++;
   }
 return r;
}


template <class OBJ> bool GList<OBJ>::Found(OBJ* item, int& idx) {
 //search the list by using CompareProc (if defined)
 //or == operator for a non-sortable list
 //for sorted lists, even when the result is false, the idx is
 //set to the closest matching object!
 int i;
 idx=-1;
 if (this->fCount==0) { idx=0;return false;}
 if (SORTED) { //binary search based on CompareProc
   //do the simple test first:

   if ((*fCompareProc)(this->fList[0],item)>0) {
                       idx=0;
                       return false;
                       }
   if ((*fCompareProc)(item, this->fList[this->fCount-1])>0) {
                       idx=this->fCount;
                       return false;
                       }

   int l, h, c;
   l = 0;
   h = this->fCount - 1;
   while (l <= h) {
       i = (l + h) >> 1;
       c = (*fCompareProc)(this->fList[i], item);
       if (c < 0)  l = i + 1;
         else {
            h = i - 1;
            if (c == 0) {
                 idx=i;
                 return true;
                }
            }
       } //while
   idx = l;
   return false;
   }
 else {//not sorted: use linear search
   // needs == operator to compare user defined objects !
   i=0;
   while (i<this->fCount) {
      if (*this->fList[i]==*item) {
         idx=i;
         return true;
         }
      i++;
      }
   return false;
   }
}

template <class OBJ> void GList<OBJ>::sortInsert(int idx, OBJ* item) {
 //idx must be the new position this new item must have
 //so the allowed range is [0..fCount]
 //the old idx item all the above will be shifted to idx+1
 if (idx<0 || idx>this->fCount) GError(SLISTINDEX_ERR, idx);
 if (this->fCount==this->fCapacity) {
    GPVec<OBJ>::Grow(idx, item);
    //expand and also copy/move data and insert the new item
    return;
    }
 //room still left, just move data around and insert the new one
 if (idx<this->fCount) //copy/move pointers only!
      memmove(&(this->fList[idx+1]), &(this->fList[idx]), (this->fCount-idx)*sizeof(OBJ*));
 this->fList[idx]=item;
 this->fCount++;
}

template <class OBJ> void GPVec<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 if (idx<0 || idx>fCount) GError(SLISTINDEX_ERR, idx);
 if (fCount==fCapacity) {
   Grow(idx, item);
   return;
   }
 if (idx<fCount)
      memmove(&fList[idx+1], &fList[idx], (fCount-idx)*sizeof(OBJ*));
 fList[idx]=item;
 fCount++;
}

template <class OBJ> void GList<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 BE_UNSORTED; //cannot do that with a sorted list!
 GPVec<OBJ>::Insert(idx,item);
}

template <class OBJ> void GPVec<OBJ>::Move(int curidx, int newidx) {
 //BE_UNSORTED; //cannot do that in a sorted list!
 if (curidx!=newidx || newidx>=fCount)
     GError(SLISTINDEX_ERR, newidx);
 OBJ* p;
 p=Get(curidx);
 //this is a delete:
 fCount--;
 if (curidx<fCount)
    memmove(&fList[curidx], &fList[curidx+1], (fCount-curidx)*sizeof(OBJ*));
 //-this was instead of delete
 Insert(newidx, p);
}

template <class OBJ> void GList<OBJ>::Move(int curidx, int newidx) {
 BE_UNSORTED; //cannot do this in a sorted list!
 GPVec<OBJ>::Move(curidx,newidx);
}

template <class OBJ> void GPVec<OBJ>::Put(int idx, OBJ* item) {
 //WARNING: this will never free the replaced item!
 TEST_INDEX(idx);
 fList[idx]=item;
}

template <class OBJ> void GList<OBJ>::Put(int idx, OBJ* item, bool re_sort) {
 //WARNING: this will never free the replaced item!
 // this may BREAK the sort order unless the "re_sort" parameter is given
 if (idx<0 || idx>this->fCount) GError(SLISTINDEX_ERR, idx);
 this->fList[idx]=item;
 if (SORTED && item!=NULL && re_sort) Sort(); //re-sort
}


template <class OBJ> void GPVec<OBJ>::Forget(int idx) {
 TEST_INDEX(idx);
 fList[idx]=NULL; //user should free that somewhere else
}

template <class OBJ> void GPVec<OBJ>::freeItem(int idx) {
  TEST_INDEX(idx);
  if (fFreeProc!=NULL) {
      (*fFreeProc)(fList[idx]);
      }
    else this->DefaultFreeProc(fList[idx]);
  fList[idx]=NULL;
}

template <class OBJ> void GPVec<OBJ>::Delete(int index) {
 TEST_INDEX(index);
 if (fFreeProc!=NULL && fList[index]!=NULL) {
   (*fFreeProc)(fList[index]); //freeItem
   }
 fList[index]=NULL;
 fCount--;
 if (index<fCount) //move higher elements if any
   memmove(&fList[index], &fList[index+1], (fCount-index)*sizeof(OBJ*));
}

//Stack usage:
template <class OBJ> OBJ* GPVec<OBJ>::Pop() {
 if (fCount<=0) return NULL;
 fCount--;
 OBJ* o=fList[fCount];
 fList[fCount]=NULL;
 return o;
}

//Queue usage:
template <class OBJ> OBJ* GPVec<OBJ>::Shift() {
 if (fCount<=0) return NULL;
 fCount--;
 OBJ* o=fList[0];
 if (fCount>0)
   memmove(&fList[0], &fList[1], (fCount)*sizeof(OBJ*));
 fList[fCount]=NULL; //not that it matters..
 return o;
}

template <class OBJ> int GList<OBJ>::Remove(OBJ* item) {
//removes an item if it's in our list
 int result=IndexOf(item);
 if (result>=0) GPVec<OBJ>::Delete(result);
 return result;
}

//linear search for the pointer address
template <class OBJ> int GPVec<OBJ>::RemovePtr(pointer item) {
if (item==NULL) return -1;
for (int i=0;i<fCount;i++)
   if ((pointer)fList[i] == item) {
       Delete(i);
       return i;
       }
return -1; //not found
}

template <class OBJ> void GPVec<OBJ>::Pack()  {//also frees items!
 for (int i=fCount-1; i>=0; i--)
    if (fList[i]==NULL) Delete(i); //shift rest of fList content accordingly
}

template <class OBJ> void GPVec<OBJ>::setCount(int NewCount) {
  if (NewCount<0 || NewCount > MAXLISTSIZE)
     GError(SLISTCOUNT_ERR, NewCount);
  if (NewCount > fCapacity) setCapacity(NewCount);
  if (NewCount > fCount)
    memset(fList[fCount], 0, (NewCount - fCount) * sizeof(OBJ*));
  fCount = NewCount;
}

template <class OBJ> void GList<OBJ>::QuickSort(int L, int R) {
 int I, J;
 OBJ* P;
 OBJ* T;
 do {
    I = L;
    J = R;
    P = this->fList[(L + R) >> 1];
    do {
      while (fCompareProc(this->fList[I], P) < 0) I++;
      while (fCompareProc(this->fList[J], P) > 0) J--;
      if (I <= J) {
        T = this->fList[I];
        this->fList[I] = this->fList[J];
        this->fList[J] = T;
        I++;
        J--;
        }
      }
    while (I <= J);
    if (L < J) QuickSort(L, J);
    L = I;
    }
 while (I < R);

}

template <class OBJ> void GList<OBJ>::Sort() {
 if (this->fList!=NULL && this->fCount>0 && fCompareProc!=NULL)
     QuickSort(0, this->fCount-1);
}

//---------------------------------------------------------------------------
#endif
