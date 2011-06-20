#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "GEAR.h"

/** Helper class that creates (is) a simple map of hardware channel to position of the pad center (x,y) 
 *  from a text file.
*/
class ChannelPosMap : public  std::map< std::pair<int,int> , gear::Vector2D> { 
  
public:
  
  
  ChannelPosMap(const std::string& file){
    
    std::ifstream f( file.c_str() , std::ios::in ) ;

    if( !  f.good() ){

      std::string mess( " \n cannot open file with channel to position mapping : \n   " ) ;
      mess +=  file ;
      throw Exception( mess ) ; 
    }

    while ( f.good() ){
      
      std::string s ;

      std::getline( f , s )  ;
      
      if( s[0] != '#' ){

	std::stringstream ss( s )  ;
      
	int rcu, channel ; 
	float x,y;

	ss >> rcu >> channel >> x  >> y ;
	
	//std::cout << "channel : " << channel << " rcu " << rcu << " x: " << x << " y: " << y << std::endl ;

	this->insert( std::make_pair( std::make_pair( channel, rcu )  , gear::Vector2D( x, y ) ) ) ;
      }


    }
    f.close() ;
    
  }
  
} ;


