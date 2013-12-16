// RWIS.h
// bridgemodel

#include<vector>
using std::vector;
typedef struct _RWISdata
{
  float surfaceTemp, airTemp, wind;
  int time;
}RWISdata;

class allRWIS
{
public:
  allRWIS(const char* RWISfile);
  ~allRWIS();
  RWISdata DataAtTime( int );
  RWISdata DataAtMinutesFromStart( int );
  int TotalMinutesInFile();

private :

 vector < RWISdata > dataFromFile;
};
