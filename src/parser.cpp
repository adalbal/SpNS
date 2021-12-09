#include "config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "parser.h"

//=============================================================================
// Getting next line from file <pf>
//=============================================================================
void GetNextLine(char* pline,FILE* pf){
    do {
        if(fgets(pline,1000,pf)==NULL) pline[0]='\0';
    } while(pline[0]=='#');
    char* ppline = pline;
    while(*ppline>=0x20 || *ppline==9) ppline++;
    *ppline = 0;
}


//=============================================================================
// Getting next line from buffer <pf>
//=============================================================================
void GetNextLine(char* pline,char* &pf){
    do {
        char* ppline = pline;
        if(!*pf) { *pline = 0; return; }
        do *(ppline++) = *(pf++); while((pf[-1]>=0x20 || pf[-1]==9 || pf[-1]<0) && (pf[-1]!='&'));
        if(!pf[-1]) pf--;
        *(ppline-1)=0;
    } while((pline[0]==0 || pline[0]=='#') && *pf);
}


//=============================================================================
// Converting strint <str> to lowercase
//=============================================================================
void toLowerCase(char* str){
    char differ = 'A'-'a';
    char ch;
    int ii = strlen(str);
    for (int i=0; i <ii; i++){
        strncpy(&ch,str+i,1);
        if(ch>='A' && ch<='Z'){
            ch = ch-differ;
            str[i] = ch;
        }
    }
}

//=============================================================================
// Getting next word (array of synbols) from string and writing it to <word>
//=============================================================================
bool GetNextWord(char *word,char * &pline, int lowerCase){
    if( pline==NULL ) return false;
    int j=0, m=0, p=0, wrdLen=0, numSymbols=0;
    bool endOfWord = false;
    for(j=0; !endOfWord; j++){
        if( ((pline[j]==' ')||(pline[j]=='\t'))&&(numSymbols==0) ) continue;
        if( (pline[j]==' ')||(pline[j]=='\0')||(pline[j]=='\n')||(pline[j]=='\r')||(pline[j]=='\t') ){
            for(p=m; p<j; p++){
                if( (pline[p]==' ')||(pline[p]=='\t') ) continue;
                word[wrdLen] = pline[p];
                wrdLen++;
            }
            word[wrdLen] = '\0';
            endOfWord = true;
        }
        else{
            numSymbols++;
        }
    }

    if(lowerCase) if( wrdLen>0 ) toLowerCase(word);
    if(wrdLen==0 || word[0]=='#'){ pline = NULL; return false; }
    for(int i=j-1; (pline[i]==' ')||(pline[i]=='\0')||(pline[i]=='\t'); i++){
        j++;
        if( pline[i]=='\0' ){ pline = NULL; return true; }
    }
    pline = pline+(j-1);
    return true;
}

 
//=============================================================================
// Comparing two words
//=============================================================================
bool CompareWords(const char* word1, const char* word2){
    char *buf1 = GimmeMem<char>(strlen(word1)+1);
    strcpy(buf1,word1);
    toLowerCase(buf1);
    char *buf2 = GimmeMem<char>(strlen(word2)+1);
    strcpy(buf2,word2);
    toLowerCase(buf2);
    int retval = strcmp(buf1,buf2);
    delete[] buf1; delete []buf2;
    return retval==0;
}
//=============================================================================
// Converting word (array of symbols) to integer
//=============================================================================
int GetIntFromWord(const char* word){
    int Result; // number which will contain the result
    int Succeeded = sscanf( word, "%d", &Result );
    if( !Succeeded || Succeeded == EOF ){ // check if something went wrong during the conversion
        crash("GetIntFromWord: Can not convert %s to integer!\n",word);
        Result = 0;
    }
    return Result;
}

//=============================================================================
// Converting word (array of symbols) to double
//=============================================================================
double GetDoubleFromWord(const char* word){
    double Result; // number which will contain the result
    int Succeeded = sscanf( word, "%lf", &Result );
    if( !Succeeded || Succeeded == EOF ){ // check if something went wrong during the conversion
        crash("GetDoubleFromWord: Can not convert %s to double!\n",word);
        return 0;
    }
    return Result;
}
//=============================================================================
// Getting integer parameter <par> with from <line>
//=============================================================================
bool GetIntFromLine(int &par,char * &line){
    char word[1000];
    char *pline = line;
    if( GetNextWord(word,pline) ){
        par = GetIntFromWord(word);
        line = NULL;
        return true;
    }
    else{
        line = pline;
        return false;
    }
}

//=============================================================================
// Getting two integer parameters <par1>,<par2> from <line>
//=============================================================================
bool GetIntIntFromLine(int &par1,int &par2,char * &line){
    char word[1000];
    char *pline = line;
    if( GetNextWord(word,pline) ){
        par1 = GetIntFromWord(word);
        if( GetNextWord(word,pline) ){
            par2 = GetIntFromWord(word);
        }
        else crash("Can not read two integer parameters !\n");
        line = NULL;
        return true;
    }
    else{
        line = pline;
        return false;
    }
}

//=============================================================================
// Getting double array <par> with size <size> from <line>
//=============================================================================
bool GetDoubleArrayFromLine(double *par,int size,char * &line){
    char word[1000];
    char *pline = line;
    for(int i=0;i<size;i++){
        if( GetNextWord(word,pline) )  par[i] = GetDoubleFromWord(word);
        else{
            line = pline;
            return false;
        }
    }
    line = NULL;
    return true;
}

//=============================================================================
// Realization methods of classes for initialization of parameters
//=============================================================================
void CLimitation::operator =(const CLimitation &aLimitation){
    lowLimit = aLimitation.lowLimit;
    symLowLim = aLimitation.symLowLim;
    symUpLim = aLimitation.symUpLim;
    upLimit = aLimitation.upLimit;
    parType = aLimitation.parType;
}

//=============================================================================
// Setting the limitations of the parameter
//=============================================================================
void CLimitation::SetLimits(const int _symLowLim,const double _lowLim,const int _symUpLim,const double _upLim){
    if( _symLowLim!=UNLIM && _symUpLim!=UNLIM && _lowLim>_upLim ) crash(
            "SetLimits: lower limit is greater upper limit!\n");
    if( _symLowLim==LT || _symLowLim==LE ) crash("SetLimits: wrong symbol for lower limit!\n");
    if( _symUpLim==GT || _symUpLim==GE ) crash("SetLimits: wrong symbol for upper limit!\n");

    parType = TYPE_DOUBLE;

    switch( _symLowLim ){
    case UNLIM: symLowLim = UNLIM; break;
    case GT:    symLowLim = GT; lowLimit=_lowLim; break;
    case GE:    symLowLim = GE; lowLimit=_lowLim; break;
    }
    switch( _symUpLim ){
    case UNLIM: symUpLim = UNLIM; break;
    case LT:    symUpLim = LT; upLimit=_upLim; break;
    case LE:    symUpLim = LE; upLimit=_upLim; break;
    }
}


//=============================================================================
// Setting the limitations of the parameter
//=============================================================================
void CLimitation::SetLimits(const int _symLowLim,const int _lowLim,const int _symUpLim,const int _upLim){
    SetLimits(_symLowLim,static_cast<double>(_lowLim),_symUpLim,static_cast<double>(_upLim));
    parType = TYPE_INT;
}

//=============================================================================
// Checking the value of parameter
//=============================================================================
bool CLimitation::CheckValue(const double value){
    bool lowGood = false;
    bool upGood = false;
    if( symLowLim==UNLIM ) lowGood = true;
    if( symLowLim==GT && (value>lowLimit) ) lowGood = true;
    if( symLowLim==GE && (value>lowLimit || fabs(value-lowLimit)<tiny) ) lowGood = true;
    if( symUpLim==UNLIM ) upGood = true;
    if( symUpLim==LT && (value<upLimit)) upGood = true;
    if( symUpLim==LE && (value<upLimit || fabs(value-upLimit)<tiny) ) upGood = true;
    return (lowGood && upGood);
}

//=============================================================================
// Checking the value of parameter
//=============================================================================
bool CLimitation::CheckValue(const int value){
    bool lowGood = false;
    bool upGood = false;
    if( symLowLim==UNLIM ) lowGood = true;
    if( symLowLim==GT && value>round_noisette(lowLimit) ) lowGood = true;
    if( symLowLim==GE && value>=round_noisette(lowLimit) ) lowGood = true;
    if( symUpLim==UNLIM ) upGood = true;
    if( symUpLim==LT && value<round_noisette(upLimit) ) upGood = true;
    if( symUpLim==LE && value<=round_noisette(upLimit) ) upGood = true;
    return (lowGood && upGood);
}

//=============================================================================
// Copy constructor
//=============================================================================
CParameter::CParameter(const CParameter &aPar){
    ptrDouble = aPar.ptrDouble;
    ptrInt = aPar.ptrInt;
    ptrChar = aPar.ptrChar;
    limit = aPar.limit;
    for(int i=0;i<256;i++)parName[i]=aPar.parName[i]; // hardcode constant is not nice!
    isSet = aPar.isSet;
    crashIt = aPar.crashIt;
    parType = aPar.parType;
}

//=============================================================================
// Setting parameters attributes
//=============================================================================
void CParameter::Set(double &par, const char *_parName,const int _parType,const int _crashIt, const int _symLowLim,
                     const double _lowLim, const int _symUpLim,
                     const double _upLim){
    ptrDouble = &par;
    if(strlen(_parName)>MaxStrLen){crash("CParameter::Set: line too long %d \n",strlen(_parName));return;}
    strcpy(parName,_parName);
    limit.SetLimits(_symLowLim,_lowLim,_symUpLim,_upLim);
    crashIt = _crashIt;
    parType = (char) _parType;
}

//=============================================================================
// Setting parameters attributes
//=============================================================================
void CParameter::Set(int &par, const char *_parName,const int _parType,const int _crashIt, const int _symLowLim,
                     const int _lowLim, const int _symUpLim,
                     const int _upLim){
    ptrInt = &par;
    if(strlen(_parName)>MaxStrLen){crash("CParameter::Set: line too long %d \n",strlen(_parName));return;}
    strcpy(parName,_parName);
    limit.SetLimits(_symLowLim,_lowLim,_symUpLim,_upLim);
    crashIt = _crashIt;
    parType = (char) _parType;
}

//=============================================================================
// Setting parameters attributes
//=============================================================================
void CParameter::Set(char *par, const char *_parName,const int _parType,const int _crashIt){
    ptrChar = par;
    if(strlen(_parName)>MaxStrLen){crash("CParameter::Set: line too long %d \n",strlen(_parName));return;}
    strcpy(parName,_parName);
    limit.SetLimits(UNLIM,0,UNLIM,0);
    crashIt = _crashIt;
    parType = (char) _parType;
}

//=============================================================================
// Getting parameter from line from txt file
//=============================================================================
bool CParameter::GetParameter(char *line){
    double numerator = 1.0, denominator = 1.0;
    char word[1000];
    char* pline = line;
    if(!GetNextWord(word,pline)) return false;
    if(!CompareWords(word+(parType==TYPE_BOOL),parName)) return false;

    bool Result = true;
    switch( parType ){
    case TYPE_INT:
        if(!GetNextWord(word,pline)) crash("Can not read parameter \"%s\"!\n",parName);
        *ptrInt = GetIntFromWord(word);
        if( !limit.CheckValue(*ptrInt) )
            crash("GetParameter: parameter %s=%d is out of limitations!\n", parName, *ptrInt);
        pprintf("%30s: %d\n",parName,*ptrInt);
        break;
    case TYPE_DOUBLE:
        if(!GetNextWord(word,pline)) crash("Can not read parameter \"%s\"!\n",parName);
        *ptrDouble = GetDoubleFromWord(word);
        if( !limit.CheckValue(*ptrDouble) )
            crash("GetParameter: parameter %s=%.2f is out of limitations!\n",parName,*ptrDouble);
        if(*ptrDouble < 0.1) pprintf("%30s: %e\n",parName,*ptrDouble);
        else pprintf("%30s: %.2lf\n",parName,*ptrDouble);
        break;
    case TYPE_STRING:
        sprintf(ptrChar, "%s", pline);
        pprintf("%30s: %s\n",parName, pline);
        break;
    case TYPE_WORD:
        if(!GetNextWord(word,pline,0)) crash("Can not read parameter \"%s\"!\n",parName);
        strcpy(ptrChar, word);
        pprintf("%30s: %s\n",parName, word);
        break;
    case TYPE_RATIO:
        if(!GetNextWord(word,pline)) crash("Can not read parameter \"%s\"!\n",parName);
        numerator = GetDoubleFromWord(word);
        pprintf("%30s: %2.1lf/",parName,numerator);
        if(!GetNextWord(word,pline)) crash("Can not read parameter \"%s\"!\n",parName);
        denominator = GetDoubleFromWord(word);
        pprintf("%2.1lf\n",denominator);
        if( fabs(denominator) < tiny ) crash("Can not read parameter \"%s\": denominator=0!\n",parName);
        *ptrDouble = numerator / denominator;
        if( !limit.CheckValue(*ptrDouble) )
            crash("GetParameter: parameter %s=%.2f is out of limitations!\n",parName,*ptrDouble);
        break;
    case TYPE_BOOL:
        if(word[0]=='+') { *ptrInt=1; pprintf("%30s: (enabled)\n",parName); }
        else if(word[0]=='-') { *ptrInt=0; pprintf("%30s: (disabled)\n",parName); }
        else Result = false;
        break;
    default:
        crash("GetParameter error: unknown parType %i!\n", parType);
    }
    if( Result && isSet ) crash("GetParameter: double definition of parameter %s!\n",parName);
    if( Result && !isSet ) isSet = true;
    return Result;
}

//=============================================================================
// Printing parameter to log file
//=============================================================================
void CParameter::PrintParameter(){
    switch( parType ){
    case TYPE_INT:
        pprintf("%30s: %d\n",parName,*ptrInt); break;
    case TYPE_DOUBLE:
        pprintf("%30s: %.2f\n",parName,*ptrDouble); break;
    case TYPE_STRING:
    case TYPE_WORD:
        pprintf("%30s: %s\n",parName,ptrChar); break;
    case TYPE_RATIO:
        pprintf("%30s: %.2f\n",parName,*ptrDouble); break;
    case TYPE_BOOL:
        pprintf("%30s: %s\n",parName, *ptrInt ? "true" : "false"); break;
    }
}

//=============================================================================
// Reading params from text file
// <fileName> -- name of the file
// <CrashIt> -- flag to crash if file doesn't exist
//=============================================================================
int CParamInitManager::ReadParamsFromFile(const char *fileName, int CrashIt){
    FileBuffer FB(fileName, CrashIt);
    if(FB.V==NULL) return 1;
    ReadParamsFromBuffer(FB);
    FB.PrintKeywords();
    return 0;
}
//=============================================================================

//=============================================================================
// Reading all file to the buffer. Returns 1 if error
int FileBuffer::LoadFile(const char* fName, int CrashIt) {
//=============================================================================
    FILE* pFILE;
    char fileName[1000];
    if(size>0) crash("FileBuffer: can't read to buffer because it is in use\n");

    sprintf(fileName,"%s",fName);
    pFILE = fopen(fileName, "rb");
    if(!pFILE) { if(CrashIt)crash("param file %s not found\n",fName); return 0; }
    pprintf("FileBuffer: loading %s\n",fName);

// Gets size and allocates mem
    fseek (pFILE, 0, SEEK_END);
    size = ftell(pFILE) + 1;
    fseek (pFILE, 0, SEEK_SET);
    V = GimmeMem<char>(size, "parser");

// Reading to mem
    for(int i=0; i<size-1; i++) V[i] = (char)fgetc(pFILE);
    V[size-1] = 0; // eof
    fclose(pFILE);

    Zone = V; NextZone = NULL;
    return 0;
}
//=============================================================================

//=============================================================================
int FileBuffer::PrintKeywords() {
//=============================================================================
    char line[1000];
    char word[1000];
    char *pbuf = Zone;
    int printed = 0;
    while(*pbuf) {
        GetNextLine(line,pbuf);
        char* pline = line;
        if(GetNextWord(word, pline, 0)) {
            pprintf0("=UNKNOWN KEYWORD= %s\n", word);
            printed++;
        }
    }
    return printed;
}
//=============================================================================

//=============================================================================
void CParamInitManager::ReadParamsFromBuffer(FileBuffer& FB, int CommentIt){
    if( CParamsVectorBase::empty() ) crash("ReadParamsFromFile: params list is empty!\n");
    if(FB.Zone!=NULL) {
        char line[1000];
        char *pbuf = FB.Zone;
        while(*pbuf) {
            char* pbufold = pbuf;
            GetNextLine(line,pbuf);
            for(CParamsVectorBase::iterator i=CParamsVectorBase::begin(); i!=CParamsVectorBase::end(); ++i)
                if( i->GetParameter(line) ) {
                    if(CommentIt) {
                        while(pbufold!=pbuf) {
                            if(*pbufold > 0x20) *pbufold = '#';
                            pbufold++;
                        }
                    }
                    break;
                }
        }
    }
    int printed = 0;
    for(CParamsVectorBase::iterator i=CParamsVectorBase::begin(); i!=CParamsVectorBase::end(); ++i)
        if( !i->GetIsSet() ){
            if(i->GetParType()!=TYPE_UNDEFTYPE && !printed) {
                pprintf("Default used for parameters:\n");
                printed = 1;
            }
            if( !i->GetCrashIt() ) i->PrintParameter();
            else crash("Parameter %s is not set!\n",i->GetParName());
        }
}
//=============================================================================


//=============================================================================
CParameter& CParamInitManager::operator[] (const char* parname) {
//=============================================================================
    CParamsVectorBase::iterator i;
    for(i=CParamsVectorBase::begin(); i!=CParamsVectorBase::end(); ++i)
        if(CompareWords(i->GetParName(), parname)) return *i;
    crash("CParamInitManager: error, can't find parameter %s in the list!\n", parname);
    return *i;
}
//=============================================================================


 
