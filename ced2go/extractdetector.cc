#include "lcio.h"

#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCRunHeader.h" 

#include <cstdlib>

//static const char* FILEN = "recjob.slcio" ; // default file name 
static std::vector<std::string> FILEN ; 

using namespace std ;
using namespace lcio ;

int main(int argc, char** argv ){
    if(argc!=2){
        cout << " usage:  extractdetector <input-file>" << endl ;
        exit(1) ;
    }

    FILEN.push_back(argv[1]);
  
    LCReader* lcReader = LCFactory::getInstance()->createLCReader() ;
  
    lcReader->open( FILEN ) ;
  
    LCRunHeader *runHdr ;
  
    try{  
        while((runHdr = lcReader->readNextRunHeader()) != 0){
            cout << runHdr->getDetectorName() << endl;
        }
    }catch(IOException& e){
        cout << " io error when reading run data : " << e.what() << endl ;
    }
  
    lcReader->close() ;
    return 0 ;
}

  
