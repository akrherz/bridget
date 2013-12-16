//  RWIS.cpp
//  bridgemodel


#include<fstream>
#include<iostream>
#include <stdlib.h>
using namespace std;
#include "RWIS.h"


#define GetDaysInMonth( year, month, daysInMonth ) 		\
   switch ( month )					\
   {									\
      case 1:						\
      case 3:								\
      case 5:								\
      case 7:							\
      case 8:							\
      case 10:							\
      case 12:								\
         daysInMonth = 31;					\
         break;								\
      case 4:								\
      case 6:								\
      case 9:							\
      case 11:					\
         daysInMonth = 30;					\
         break;								\
      case 2:								\
         if (year % 4 == 0 &&				\
            (year % 100 != 0 || year % 400 == 0))		\
         {							\
            daysInMonth = 29;    \
         }								\
         else								\
         {								\
            daysInMonth = 28;			\
         }							\
         break;							\
   }

allRWIS::allRWIS(const char* RWISfile)
{  // this function takes in the RWIS file and converts its dates into minutes from the year 2000.
   ifstream rwisfile( RWISfile );
   RWISdata temporary;
   int year, month, day, hour, minute, minutes;
   long long date;
   int count = 0;

   if( rwisfile.good() == false )
   {
      cerr << "Can't open file " << RWISfile << endl;
      exit( 11 );
   }

   rwisfile >> date >> temporary.airTemp >> temporary.surfaceTemp >> temporary.wind;
   
   while ( rwisfile.good() )
   {
      //if bad rwis stuff comes up it skips it and goes to the next line of input
      while (temporary.airTemp == -99999 || temporary.surfaceTemp == -99999 || temporary.wind == -99999)
      {
         rwisfile >> date >> temporary.airTemp >> temporary.surfaceTemp >> temporary.wind;
         //cout<<"hi"<<temporary.airTemp<<"  "<<temporary.surfaceTemp<<"  "<<temporary.wind<<endl;
      }

      if ( temporary.surfaceTemp > 150 || temporary.surfaceTemp < -50 )
      {
         cerr << " bad RWIS surface temperature input  " <<temporary.surfaceTemp<< endl;
         exit(12);
      }
      if ( temporary.airTemp > 110 || temporary.airTemp < -50 )
      {
         cerr << " bad RWIS air temperature input  " <<temporary.airTemp<< endl;
         exit(13);
      }

      if ( temporary.wind > 110 || temporary.wind < 0 )
      {
         cerr << " bad RWIS wind input  " <<temporary.wind<< endl;
         exit(21);
      }

      minute = date % 100;
      date -= minute;
      date = date / 100;
      hour = date % 100;
      date -= hour;
      date = date / 100;
      day = date % 100;
      date -= day;
      date = date / 100;
      month = date % 100;
      date -= month;
      date = date / 100;
      year = date;

      minutes = 0;
      temporary.wind = temporary.wind * 0.515 ;
      temporary.surfaceTemp = (5/9.0) * (temporary.surfaceTemp - 32) + 273;
      temporary.airTemp = (5/9.0) * (temporary.airTemp - 32) + 273;
      // converts the date into the number of minutes from the year 2000
      month --;
      while ( year % 1999 > 0 )
      {
          while (month > 0 )
          {
             int daysInMonth;
             GetDaysInMonth(year, month, daysInMonth );
             minutes += daysInMonth * 1440;
             month --;
          }
         year --;
         month = 12;
      }
      minutes += day * 1440;
      minutes += hour * 60;
      minutes += minute;
      temporary.time = minutes;
      
      if ( dataFromFile.empty() == false && temporary.time < dataFromFile.back().time )
      {
         cerr << "RWIS data out of order  " << temporary.time<<"  "<<dataFromFile.back().time<<endl;
          exit (14);
      }
      
      
      dataFromFile.push_back( temporary );
      rwisfile >> date >> temporary.airTemp >> temporary.surfaceTemp >> temporary.wind;
      count++;
   }
   rwisfile.close();

   if( dataFromFile.empty() )
   {
      cerr << "Couldn't get RWIS data." << endl;
      exit( 15 );
   }
   
}

allRWIS::~allRWIS()//terminator
{ }

RWISdata allRWIS::DataAtTime( int requested)
{  // determines if a requested time is valid, needs interpolation, or can be retreived exactly. sends back to 
   // DataAtMinutesFromStart
   RWISdata entryBefore, entryAfter, returnValue; 
   unsigned int counter = 0;
   while (counter < dataFromFile.size() && dataFromFile[counter].time < requested)
   {
      counter++;
   }
   if (counter ==0 && requested < dataFromFile[0].time)
   {
      cerr << "requested time is before beginning of rwis file  " << requested<<"  " <<dataFromFile[0].time<<endl;
       exit (16);
   }
   
      if ( counter > 2 && requested > (dataFromFile[counter-1 ].time + 300))
      {
         cerr << "RWIS data too sparse -- 5 hour gap detected  "<< requested << "  " 
	      << dataFromFile[counter -1].time<< endl;
	 exit (17);
      }
   
   if (counter == dataFromFile.size())
   {
      cerr << " requested time beyond end of rwis file" << endl;
      exit (18);
   }
   if ( requested == dataFromFile[counter].time )
   {   // they are equal
      return dataFromFile[counter];
   }

   entryAfter = dataFromFile[counter];
   entryBefore = dataFromFile[counter-1];

   float airSlope, surfaceSlope, windSlope;
   int gapBeforeAfter, gapBeforeRequested;
   gapBeforeAfter = entryAfter.time - entryBefore.time;
   gapBeforeRequested = requested - entryBefore.time;

   airSlope = (entryAfter.airTemp - entryBefore.airTemp) / (float)gapBeforeAfter;
   surfaceSlope = (entryAfter.surfaceTemp - entryBefore.surfaceTemp) / (float)gapBeforeAfter;
   windSlope = (entryAfter.wind - entryBefore.wind) / (float)gapBeforeAfter;

   returnValue.airTemp = airSlope * gapBeforeRequested + entryBefore.airTemp;
   returnValue.surfaceTemp = surfaceSlope * gapBeforeRequested + entryBefore.surfaceTemp;
   returnValue.wind = windSlope * gapBeforeRequested + entryBefore.wind;
   returnValue.time = requested;

   return returnValue;
}


RWISdata allRWIS::DataAtMinutesFromStart( int minutes)
{  // finds the number of minutes from 2000 of the initialization timestep requested by main, calls a function to interpolate the meteorological values, and returns the next minute values to main
   if (dataFromFile.empty() )
   {
      cerr << "no info from RWIS file" << endl;
      exit (19);
   }

   int timeRequested = dataFromFile[0].time;
   timeRequested += minutes;

   return DataAtTime( timeRequested );
}

int allRWIS::TotalMinutesInFile()
// finds total minutes in file
{
   int beginTime;
   int endTime;
   beginTime = dataFromFile[0].time;
   endTime = dataFromFile.back().time;

   return ( endTime - beginTime );
}

