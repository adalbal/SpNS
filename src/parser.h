#ifndef PARSER_HEADER
#define PARSER_HEADER

#include <vector>

using namespace std;

//=============================================================================
// Memory interface
//=============================================================================
template <typename ValueType>
inline ValueType* GimmeMem(int N, const char* label = NULL){
    if(N==0) return NULL;
    if(N<0) crash("GimmeMem: wrong input N = %i, label = %s\n", N, label==NULL ? "" : label);
    ValueType* p=NULL;
    try{p = new ValueType[N];}
    catch (bad_alloc &ba){
        if(p!=NULL) delete[] p;
        crash("GimmeMem: memory fuckup! int %d, name %s, text: %s\n", N, label?label:"noname",ba.what());
    }
    return p;
}
 
//=============================================================================
// Default data type for block 2D arrays NxM on a base of flat array
//=============================================================================
template <typename ValueType> class BlockArray {
public:
    int N; // number of blocks
    int M; // size of block
    int alias; // является ли ссылкой на другой BlockArray - нужно для OpenMP, чтобы задействовать исходный экземпляр CE
    ValueType * V; // data array

    BlockArray(){M=0;N=0;alias=0;V=NULL;}
    BlockArray(int n, int m){
        if((n<=0)||(m<=0))crash("BlockArray constructor: wrong input %d %d \n", n,m);
        V = GimmeMem<ValueType>(n*m, "BlockArrayConstructor"); M=m; N=n; alias=0;
    }
    ~BlockArray(){ if(M>0 && !alias && V) FreeMem(V); }

    inline bool Allocated(){return (M>0); }
    
    inline int Alloc(int n, int m, const char* label = NULL){
        if((n<0)||(m<0))crash("BlockArray Alloc: wrong input %d %d \n", n,m);
        if(M>0)crash("BlockArray Alloc: already allocated!\n");
        V = GimmeMem<ValueType>(n*m, label); M=m; N=n;
        return 0;
    }
    inline void Dealloc(){ if((M>0)&&(alias==0)&& V) FreeMem(V); N=0; M=0; V=NULL; alias=0; }

    inline int Alias(BlockArray<ValueType>& object){ // Ссылка на другой BlockArray
        if(M>0)crash("BlockArray Alias: already allocated!\n");
        if(object.M==0)crash("BlockArray Alias: not yet allocated!\n");
        M=object.M;
        N=object.N;
        V=object.V;
        alias=1;
        return 0;
    }
    inline BlockArray <ValueType>& operator=(BlockArray<ValueType>& object){
        if(this != &object){
            if((M==0)&&(object.M>0)) Alloc(object.N,object.M); //if not yet allocated
            if(object.M!=M) crash("BlockArray: assignment of arrays with different size of block %d %d\n", object.M,M);
            if(object.N!=N) crash("BlockArray: assignment of arrays with different size\n", object.N, N);
            for(int i=0; i<N*M; i++) V[i]=object.V[i];
        }
        return *this;
    }
    inline BlockArray <ValueType>& operator+=(BlockArray<ValueType>& object){
        if(this != &object){
//            if((M==0)&&(object.M>0)) Alloc(object.N,object.M); //if not yet allocated
            if(object.M!=M) crash("BlockArray: assignment of arrays with different size of block %d %d\n", object.M,M);
            if(object.N!=N) crash("BlockArray: assignment of arrays with different size\n", object.N, N);
            for(int i=0; i<N*M; i++) V[i]+=object.V[i];
        }
        return *this;
    }
    inline BlockArray <ValueType>& operator=(ValueType x){
        if(M==0) crash("BlockArray =: not allocated!\n");
        for(int i=0; i<N*M; i++) V[i]=x;
        return *this;
    }

    inline ValueType *operator[](int i){return V+i*M;}

    // prototypes for  additional debug functionality for BlockArray (in control.cpp)
    void PrintToFile(const char* fname);
    void ReadFromFile(const char* fname);
    double CompareWithFile(const char* fname, int *imax = NULL, int write = 0);
};


#define MaxStrLen 256

#define GT 1    // ">"
#define GE 2    // ">="
#define LT 3    // "<"
#define LE 4    // "<="
#define UNLIM 5 // unlimited

#define TYPE_INT 'i'
#define TYPE_DOUBLE 'd'
#define TYPE_WORD 'w'
#define TYPE_STRING 's'
#define TYPE_RATIO 'r'
#define TYPE_UNDEFTYPE 'u'
#define TYPE_BOOL 'f'

#define IO_CRASH 1      // crash if cant open file
#define IO_DONTCRASH 0  // dont crash if cant open file

#define tiny    1e-16
#define tinyflt 1e-8
#define huge    1e+16
#define hugeflt 1e+8

// Parser interfaces
void GetNextLine(char* line, char*&pf);
void GetNextLine(char* line, FILE* pf);
bool GetNextWord(char *word,char * &pline, int lowerCase=1);

void toLowerCase(char* str);
bool CompareWords(const char* word1, const char* word2);

int GetIntFromWord(const char* word);
double GetDoubleFromWord(const char* word);

bool GetIntParameter(int &par,char *ParName,char * &line);
bool GetDoubleParameter(double &par,char *ParName,char * &line);
bool GetRatioParameter(double &par,char *ParName,char * &line);
bool GetStringParameter(char* par,char *ParName,char * &line);
bool GetIntFromLine(int &par,char * &line);
bool GetIntIntFromLine(int &par1,int &par2,char * &line);
bool GetDoubleArrayFromLine(double *par,int size,char * &line);

inline int round_noisette(double value) { return int(value + 0.5*(value>0.0 ? 1.0 : -1.0)); }

//=============================================================================
// Classes for initialization of parameters
//=============================================================================

class CLimitation{ // Class which describes the limitation of the parameter
private:
    double lowLimit;
    double upLimit;
    char parType;
    int symLowLim;
    int symUpLim;

public:
    CLimitation():
        lowLimit(0.0),upLimit(0.0),parType(TYPE_UNDEFTYPE),
        symLowLim(UNLIM),symUpLim(UNLIM){}

    CLimitation(CLimitation &aLimitation){
        lowLimit  = aLimitation.lowLimit;
        upLimit   = aLimitation.upLimit;
        symLowLim = aLimitation.symLowLim;
        symUpLim  = aLimitation.symUpLim;
    }

    void operator=(const CLimitation &aLimitation);
    void SetLimits(const int _symLowLim,const double _lowLim,const int _symUpLim,const double _upLim);
    void SetLimits(const int _symLowLim,const int _lowLim,const int _symUpLim,const int _upLim);
    bool CheckValue(const int value);
    bool CheckValue(const double value);
};


class CParameter{ // Class which describes the parameter
private:
    CLimitation limit;
    int *ptrInt;        // Pointer to the int parameter
    double *ptrDouble;  // Pointer to the double parameter
    char *ptrChar;      // Pointer to the char parameter
    int crashIt;
    char parName[MaxStrLen];
    char parType;
    bool isSet;

    inline void _init(){
        ptrInt=NULL;
        ptrDouble=NULL;
        ptrChar=NULL;
        isSet=false;
        crashIt=0;
        parType=TYPE_UNDEFTYPE;
    }

public:
    CParameter(){_init();};
    CParameter(const CParameter &aPar);

    void Set(int &par,const char *_parName,const int _parType,const int _crashIt,const int _symLowLim,
             const int _lowLim,const int _symUpLim,const int _upLim);
    void Set(double &par,const char *_parName,const int _parType,const int _crashIt,const int _symLowLim,
             const double _lowLim,const int _symUpLim,const double _upLim);
    void Set(char *par,const char *_parName,const int _parType,const int _crashIt);

    CParameter(int &par,const char *_parName,const int _parType,const int _crashIt,const int _symLowLim,
             const int _lowLim,const int _symUpLim,const int _upLim){
      _init();
      Set(par, _parName, _parType, _crashIt, _symLowLim, _lowLim, _symUpLim, _upLim);
    }
    CParameter(double &par,const char *_parName,const int _parType,const int _crashIt,const int _symLowLim,
             const double _lowLim,const int _symUpLim,const double _upLim){
      _init();
      Set(par, _parName, _parType, _crashIt, _symLowLim, _lowLim, _symUpLim, _upLim);
    }
    CParameter(char *par,const char *_parName,const int _parType,const int _crashIt){
      _init();
      Set(par,_parName,_parType,_crashIt);
    }

    inline bool CheckValue(const int value){return limit.CheckValue(value);}
    inline bool CheckValue(const double value){return limit.CheckValue(value);}
    bool GetParameter(char *line);
    void PrintParameter();
    inline bool GetIsSet(void){ return isSet; }
    inline bool GetCrashIt(void){ return crashIt==0 ? false : true ; }
    inline char* GetParName(void){ return parName; }
    inline char GetParType(void) { return parType; }
};


// Class for FileBuffer
struct FileBuffer {
    char *V; char* Zone; char* NextZone;
    int size;

    FileBuffer() { size=0; V=NULL; }
    FileBuffer(const char *fName, int CrashIt = 1) { size=0; V=NULL; LoadFile(fName, CrashIt); }
    void Close() { if(size>0) delete[]V; size = 0; V=NULL; }
    ~FileBuffer() { if(size>0) delete[]V; size = 0; V=NULL; }

    int LoadFile(const char *fName, int CrashIt = 1); // Returns 1 if file cannot be opened
    int MakeZone(const char *KeyWord); // Returns 1 if EOF
//    void Rewind() { Zone = V; NextZone = NULL; } // Не работает! Разбивка на зоны меняет содержимое
    void GotoNextZone() { Zone = NextZone; NextZone = NULL; }
    int PrintKeywords(); // Resurns number of unknown Keywords
};

// Class which manages initialization of list of parameters
class CParamInitManager:public std::vector<CParameter>
{
private:
public:
    typedef std::vector<CParameter> CParamsVectorBase;

    CParamInitManager(){};
    void RequestParameter(int &par,const char *_parName,const int _parType=TYPE_INT,const int _crashIt=IO_CRASH,
        const int _symLowLim=UNLIM,const int _lowLim=0,const int _symUpLim=UNLIM,const int _upLim=0){
        CParamsVectorBase::push_back(CParameter(par,_parName,_parType,_crashIt,_symLowLim,_lowLim,_symUpLim,_upLim));
    }
    void RequestParameter(double &par,const char *_parName,const int _parType=TYPE_DOUBLE,const int _crashIt=IO_CRASH,
        const int _symLowLim=UNLIM,const double _lowLim=0.0,const int _symUpLim=UNLIM,const double _upLim=0.0){
        CParamsVectorBase::push_back(CParameter(par,_parName,_parType,_crashIt,_symLowLim,_lowLim,_symUpLim,_upLim));
    }
    void RequestParameter(char *par,const char *_parName,const int _parType=TYPE_WORD,const int _crashIt=IO_CRASH){
        CParamsVectorBase::push_back(CParameter(par,_parName,_parType,_crashIt));
    }
    int ReadParamsFromFile(const char *fileName, int CrashIt = IO_CRASH);
    void ReadParamsFromBuffer(FileBuffer &FB, int CommentIt=1);
    CParameter& operator[] (const char* parname);
};
 
#endif
