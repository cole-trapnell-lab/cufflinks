//---------------------------------------------------------------------------
/*
Sortable collection of pointers to objects
*/

#ifndef _GVec_HH
#define _GVec_HH

#include "GBase.h"
//#include "assert.h"

#ifdef __LINE__
#define GVEC_INDEX_ERR "GVec error (%s:%d):Invalid index: %d\n"
#define TEST_INDEX(x) \
 if (x<0 || x>=fCount) GError(GVEC_INDEX_ERR, __FILE__,__LINE__, x)
#else
#define GVEC_INDEX_ERR "GVec error: invalid index: %d\n"
#define TEST_INDEX(x) \
 if (x<0 || x>=fCount) GError(GVEC_INDEX_ERR, x, __FILE__,__LINE__)
#endif

#define GVEC_CAPACITY_ERR "GVec error: invalid capacity: %d\n"
#define GVEC_COUNT_ERR "GVec error: invalid count: %d\n"

#define MAXLISTSIZE INT_MAX-1

#define FREEDATA (fFreeProc!=NULL)

//basic template for array of objects;
//so it doesn't require comparison operators to be defined
template <class OBJ> class GVec {
  protected:
    OBJ* fArray;
    int fCount;
    int fCapacity;
    void qSort(int L, int R, GCompareProc* cmpFunc);
  public:
    GVec(int init_capacity=2);
    GVec(GVec<OBJ>& array); //copy constructor
    const GVec<OBJ>& operator=(GVec<OBJ>& array); //copy operator
    virtual ~GVec();
    void idxInsert(int idx, OBJ& item);
    void Grow();
    void Grow(int idx, OBJ& item);
    void Reverse(); //WARNING: will break the sort order if SORTED!
    int Add(OBJ* item); // simply append to the end of fArray, reallocating as needed
    int Add(OBJ item) { return Add(&item); } //both will CREATE a new OBJ and COPY to it
                       // using OBJ copy operator=
    // -- stack/queue usage:
    //int Push(OBJ& item) { return Add(&item); }
    int Push(OBJ item) { return Add(&item); }
    OBJ Pop();// Stack use; removes and returns a copy of the last item
    OBJ Shift(); //Queue use: removes and returns a copy of the first item

    void Add(GVec<OBJ>& list); //append copies of all items from another list

    OBJ& Get(int idx) {
          TEST_INDEX(idx);
          return fArray[idx];
          }
    OBJ& operator[](int i) {
          TEST_INDEX(i);
          return fArray[i];
          }
    OBJ& Last() {
         TEST_INDEX(fCount-1);
         return fArray[fCount-1];
         }
    OBJ& First() {
         TEST_INDEX(0);
         return fArray[0];
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

    void setCount(int NewCount);         // will trim or expand the array as needed
    void setCount(int NewCount, OBJ& v); //same as setCount(), but if expanding, new objects are set to v
    void Resize(int NewCount) { setCount(NewCount); }
    void Resize(int NewCount, OBJ& v) { setCount(NewCount, v); }

    void Move(int curidx, int newidx);
    bool isEmpty() { return fCount==0; }
    bool notEmpty() { return fCount>0; }

    void Sort(GCompareProc* cmpFunc);
};

//---- template for dynamic array of object pointers
//---- it's faster than GVec<OBJ*> and has item deallocation awareness
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
    void qSort(int L, int R, GCompareProc* cmpFunc);
  public:  
    static void DefaultFreeProc(pointer item) {
      delete (OBJ*)item;
      }
    virtual ~GPVec();
    GPVec(int init_capacity=2, bool free_elements=true); //also the default constructor
    GPVec(GPVec<OBJ>& list); //copy constructor?
    GPVec(GPVec<OBJ>* list); //kind of a copy constructor
    const GPVec<OBJ>& operator=(GPVec<OBJ>& list);
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
    int IndexOf(pointer item); //a linear search for pointer address!
    void Sort(GCompareProc* cmpFunc);
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

template <class OBJ> GVec<OBJ>::~GVec() {
 this->Clear();
}


template <class OBJ> void GVec<OBJ>::setCapacity(int NewCapacity) {
  if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
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
        fArray[i] = oldArray[i]; // operator=
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

template <class OBJ> void GVec<OBJ>::Grow() {
 int delta;
 if (fCapacity > 64 ) {
   delta = (fCapacity > 0xFFF) ? 0x100 : (fCapacity>>4);
 }
 else {
   delta = (fCapacity>8) ? (fCapacity>>2) : 1 ;
 }
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
 /*
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
 */
 if (fCapacity > 64 ) {
   delta = (fCapacity > 0xFFF) ? 0x100 : (fCapacity>>4);
 }
 else {
   delta = (fCapacity>8) ? (fCapacity>>2) : 1 ;
 }
 int NewCapacity=fCapacity+delta;
  if (NewCapacity <= fCount || NewCapacity >= MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
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


template <class OBJ> int GVec<OBJ>::Add(OBJ* item) {
 if (item==NULL) return -1;
 if (fCount==fCapacity) Grow();
 fArray[fCount] = *item; //OBJ::operator= must copy OBJ properly!
 fCount++;
 return fCount-1;
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


//Stack usage:
template <class OBJ> OBJ GVec<OBJ>::Pop() {
 if (fCount<=0) GError("Error: invalid GVec::Pop() operation!\n");
 fCount--;
 //OBJ o(fArray[fCount]); //copy constructor
 //o=fList[fCount];
 //fArray[fCount]=NULL;
 return fArray[fCount]; //copy of the last element
}

//Queue usage:
template <class OBJ> OBJ GVec<OBJ>::Shift() {
 if (fCount<=0) GError("Error: invalid GVec::Shift() operation!\n");
 fCount--;
 OBJ o(fArray[0]); //copy constructor
 if (fCount>0)
   memmove(&fArray[0], &fArray[1], (fCount)*sizeof(OBJ));
 //fList[fCount]=NULL; //not that it matters..
 return o;
}

template <class OBJ> void GVec<OBJ>::idxInsert(int idx, OBJ& item) {
 //idx must be the new position this new item must have
 //so the allowed range is [0..fCount]
 //the old idx item all the above will be shifted to idx+1
 if (idx<0 || idx>fCount) GError(GVEC_INDEX_ERR, idx);
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


template <class OBJ> void GVec<OBJ>::Move(int curidx, int newidx) {
 if (curidx!=newidx || newidx>=fCount)
     GError(GVEC_INDEX_ERR, newidx);
 OBJ tmp=fArray[curidx]; //copy constructor here
 fArray[curidx]=fArray[newidx];
 fArray[newidx]=tmp;
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
     GError(GVEC_COUNT_ERR, NewCount);
  if (NewCount > fCapacity) setCapacity(NewCount);
  fCount = NewCount;
}

template <class OBJ> void GVec<OBJ>::setCount(int NewCount, OBJ& v) {
  if (NewCount<0 || NewCount > MAXLISTSIZE)
     GError(GVEC_COUNT_ERR, NewCount);
  if (NewCount > fCapacity)	setCapacity(NewCount);
  if (NewCount>fCount) {
	for (int i=fCount;i<NewCount;i++)
	   fArray[i]=v;
  }
  fCount = NewCount;
}


template <class OBJ> void GVec<OBJ>::qSort(int l, int r, GCompareProc* cmpFunc) {
 int i, j;
 OBJ p,t;
 do {
    i = l; j = r;
    p = this->fArray[(l + r) >> 1];
    do {
      while (cmpFunc(&(this->fArray[i]), &p) < 0) i++;
      while (cmpFunc(&(this->fArray[j]), &p) > 0) j--;
      if (i <= j) {
        t = this->fArray[i];
        this->fArray[i] = this->fArray[j];
        this->fArray[j] = t;
        i++; j--;
        }
      } while (i <= j);
    if (l < j) qSort(l, j, cmpFunc);
    l = i;
    } while (i < r);
}

template <class OBJ> void GVec<OBJ>::Sort(GCompareProc* cmpFunc) {
 if (this->fArray!=NULL && this->fCount>0 && cmpFunc!=NULL)
     qSort(0, this->fCount-1, cmpFunc);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//*=> GPVec implementation

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

template <class OBJ> GPVec<OBJ>::GPVec(int init_capacity, bool free_elements) {
  fCount=0;
  fCapacity=0;
  fList=NULL;
  fFreeProc=(free_elements) ? DefaultFreeProc : NULL;
  setCapacity(init_capacity);
}

template <class OBJ> GPVec<OBJ>::~GPVec() {
 this->Clear();//this will free the items if fFreeProc is defined
}

template <class OBJ> void GPVec<OBJ>::setCapacity(int NewCapacity) {
  if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
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

template <class OBJ> void GPVec<OBJ>::Grow() {
 int delta;
 if (fCapacity > 64 ) {
   delta = (fCapacity > 0xFFF) ? 0x100 : (fCapacity>>4);
 }
 else {
   delta = (fCapacity>8) ? (fCapacity>>2) : 1 ;
 }
  setCapacity(fCapacity + delta);
}

template <class OBJ> void GPVec<OBJ>::Grow(int idx, OBJ* newitem) {
 int delta;
 if (fCapacity > 64 ) {
   delta = (fCapacity > 0xFFF) ? 0x100 : (fCapacity>>4);
 }
 else {
   delta = (fCapacity>8) ? (fCapacity>>2) : 1 ;
 }
 // setCapacity(fCapacity + delta);
 int NewCapacity=fCapacity+delta;
  if (NewCapacity <= fCount || NewCapacity > MAXLISTSIZE)
    GError(GVEC_CAPACITY_ERR, NewCapacity);
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

template <class OBJ> int GPVec<OBJ>::Add(OBJ* item) {
 int result;
 if (item==NULL) return -1;
 result = fCount;
 if (result==fCapacity) this->Grow();
 fList[result]=item;
 fCount++;
 return fCount-1;
}

template <class OBJ> void GPVec<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 if (idx<0 || idx>fCount) GError(GVEC_INDEX_ERR, idx);
 if (fCount==fCapacity) {
   Grow(idx, item);
   return;
   }
 if (idx<fCount)
      memmove(&fList[idx+1], &fList[idx], (fCount-idx)*sizeof(OBJ*));
 fList[idx]=item;
 fCount++;
}

template <class OBJ> void GPVec<OBJ>::Move(int curidx, int newidx) {
 //BE_UNSORTED; //cannot do that in a sorted list!
 if (curidx!=newidx || newidx>=fCount)
     GError(GVEC_INDEX_ERR, newidx);
 OBJ* p;
 p=Get(curidx);
 //this is a delete:
 fCount--;
 if (curidx<fCount)
    memmove(&fList[curidx], &fList[curidx+1], (fCount-curidx)*sizeof(OBJ*));
 //-this was instead of delete
 Insert(newidx, p);
}

template <class OBJ> void GPVec<OBJ>::Put(int idx, OBJ* item) {
 //WARNING: this will never free the replaced item!
 TEST_INDEX(idx);
 fList[idx]=item;
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
     GError(GVEC_COUNT_ERR, NewCount);
  if (NewCount > fCapacity) setCapacity(NewCount);
  if (NewCount > fCount)
    memset(fList[fCount], 0, (NewCount - fCount) * sizeof(OBJ*));
  fCount = NewCount;
}

template <class OBJ> void GPVec<OBJ>::qSort(int L, int R, GCompareProc* cmpFunc) {
 int I, J;
 OBJ* P;
 OBJ* T;
 do {
    I = L;
    J = R;
    P = this->fList[(L + R) >> 1];
    do {
      while (cmpFunc(this->fList[I], P) < 0) I++;
      while (cmpFunc(this->fList[J], P) > 0) J--;
      if (I <= J) {
        T = this->fList[I];
        this->fList[I] = this->fList[J];
        this->fList[J] = T;
        I++;
        J--;
        }
      }
    while (I <= J);
    if (L < J) qSort(L, J, cmpFunc);
    L = I;
    }
 while (I < R);
}

template <class OBJ> void GPVec<OBJ>::Sort(GCompareProc* cmpFunc) {
 if (this->fList!=NULL && this->fCount>0 && cmpFunc!=NULL)
     qSort(0, this->fCount-1, cmpFunc);
}

//---------------------------------------------------------------------------
#endif
