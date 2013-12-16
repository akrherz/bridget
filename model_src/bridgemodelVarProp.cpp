// Tina Greenfield
// This program is designed to calculate a concrete bridge's surface temperature.  It is initialized 
//with RWIS bridge and air temperatures, all else from MM5.  MM5 supplies atmospheric conditions 
//at regular intervals, which are used to force the bridge temperature.  MM5 values are not 
//altered in this program.

#include <iostream.h>
#include <fstream.h>
#include <string>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>

#include "RWIS.h"

float radiative( float, float, float, float, float,float );
float convection( float );
float initialize( float, float, float, float, float, double, int, int );
float precipitation( float, float, float, float, float, float, float, float &, float &, int &, int );
float bridgelevels( float, float, float, float, float, float, float, double, int, int );


int main (int argc, const char * argv[])
{
   if (argc != 4)
   {
      cerr << "usage: bridgemodel <rwis file> <mm5 file> <bridgeProperties>" << endl;
      return 1; 
   }

   allRWIS rwisfile( argv[1] );
   ifstream mm5file( argv[2] );
   ifstream bridgeProperties( argv[3] );
   // returned values from functions residual, convection, frost, precipitation, initialize,and bridgelevels
   float radiation, h, frostDepth, bridgeTemp, latentFlux, rwisAir;
   float initialTemp;
   int condition;
   
    
   //stuff from MM5
   float oSWradiationIn,oLWradiationIn, oprecipRate, otempA, owind, oqa;//"old" values for interpolation
   float nSWradiationIn,nLWradiationIn, nprecipRate, ntempA, nwind, nqa;//"new" values for interpolation
   float cSWradiationIn,cLWradiationIn, cprecipRate, ctempA, cwind, cqa;//"current" values from interpolation

   int minutes;
   int oldMinutes;
   float difference;

   int year, month, day, hour;
   int yearO, monthO, dayO, hourO;
   int junky;
   char junk;
   //char time[15];
   //char timestamp[15]; //for use in time output
   //char extraTime[9];
   

   //stuff from the property file
   float thickness, swAbs, lwAbs, kbridge;
   double thermalDiffusivity;
   

   
   if( mm5file.good() == false )
   {
      cerr << "Can't open file " << argv[2] << endl;
      exit( 1 );
   }
   if( bridgeProperties.good() == false )
   {
      cerr << "Can't open file " << argv[3] << endl;
      exit( 41 );
   }
   ofstream outfile( "pavetemp.out" , ios::out);
   if (!outfile)
   {
      cerr << "file could not be opened: "  << endl;
      exit (2);
   }
   
   bridgeProperties>>thickness>>kbridge>>thermalDiffusivity>>swAbs>>lwAbs;
   // initializing the bridge deck layer temperatures
   int initialcounter = 1;
   for ( int i = 0; i < rwisfile.TotalMinutesInFile(); i++ )
   {
      RWISdata currentData;
      currentData = rwisfile.DataAtMinutesFromStart(i);
      bridgeTemp = currentData.surfaceTemp;
      h = convection(currentData.wind);
      
      for ( int c = 0; c < 19; c++)
      {
         initialTemp = initialize( bridgeTemp, currentData.airTemp, h, kbridge, thickness, thermalDiffusivity,
	 c, initialcounter); 	 
      }

      initialcounter = initialcounter + 1;
      rwisAir = currentData.airTemp;
   }

   // All RWIS initialization is done.
   // Beginning of loop for calculations for forecasted bridge temperature
   int linecounter = 0;
   int check = 0;
   int z = 0;
   float zz = 0;

   while( mm5file>>year>>junk>>month>>junk>>day>>junk>>hour>>junk>>junky>>junk>>junky>>minutes>>junk>>
          ntempA>>nqa>>nwind>>nSWradiationIn>>nLWradiationIn>>nprecipRate )//setw(15)>>time>>setw(9)
   { 
      if ( linecounter == 0 && fabs(ntempA - rwisAir) > 15.0 )
      {
         cerr << "15+ gap MM5 & input air temperature   " << ntempA<<"  "<<rwisAir<<endl;
         return 23;
      }
      
      // mm5 data is hourly.  It needs to be interpolated for minute by minute.
      if ( linecounter != 0 )
      {
         // slope-point for the minute by minute data
         difference = minutes - oldMinutes; 
         while ( z < minutes )
         {  
	    zz = z - oldMinutes;
 
            cqa = ( nqa - oqa ) / ((float) difference) * (zz) + oqa;
            cwind = ( nwind - owind ) / ((float) difference) * (zz) + owind;
            cSWradiationIn = ( nSWradiationIn - oSWradiationIn ) / ((float) difference) * (zz) + oSWradiationIn;
            cLWradiationIn = ( nLWradiationIn - oLWradiationIn ) / ((float) difference) * (zz) + oLWradiationIn;
            cprecipRate = ( nprecipRate - oprecipRate ) / ((float) difference) * (zz) + oprecipRate;
            ctempA = ( ntempA - otempA ) / ((float) difference) * (zz) + otempA;

            // quality checks
            if ( ctempA > 410 || ctempA < 200 )
            {
              cerr << " bad MM5 air temperature input   " <<ctempA<< endl;
              return 4;
            }
            if ( cqa > 1 || cqa < 0 )
            {
               cerr << " bad MM5 humidity input   " << cqa<<endl;
              return 5;
            }
            if ( cwind > 70 || cwind < 0 )
            {
              cerr << " bad MM5 windspeed input   " <<cwind<< endl;
              return 6;
            }
            if ( cSWradiationIn> 1000 || cSWradiationIn < 0 )
            {
              cerr << " bad MM5 SWradiation input   "<<cSWradiationIn<< endl;
              return 7;
            }
            if ( cLWradiationIn> 2000 || cLWradiationIn < 50 )
            {
               cerr << " bad MM5 LWradiation input   " << cLWradiationIn << endl;
               return 8;
            }
	    
            if ( cprecipRate > 40 || cprecipRate < 0 )
            {
              cerr << " bad MM5 precipitation input   " <<cprecipRate<< endl;
              return 9;
            }

            //calling functions
	    
            latentFlux = precipitation( cprecipRate, h, ctempA, cqa, cwind, kbridge, 
	                               thickness,bridgeTemp, frostDepth, condition,  check );

            radiation = radiative(cSWradiationIn, cLWradiationIn, bridgeTemp, swAbs, lwAbs, ctempA);

            h = convection(cwind);
            bridgeTemp = bridgelevels(radiation, latentFlux, bridgeTemp,ctempA, h, thickness, kbridge, 
	                 thermalDiffusivity, check, z );
   
            check++;
	    // calculating dew point temperature for output purposes only.
	    float Td;
	    Td = 5.42*pow(10,3) / (log( 2.53 * pow(10,8) * 0.622 / ( cqa * 100.0)));

	    if(hourO == 24)
	    {
	       hourO = 0;
	       dayO++;
	    }
            //writing the calculated bridge temperature to file.  different outputs to make timestap correct
            outfile << yearO <<"-"<<monthO<<"-"<<dayO<<"_"<<hourO <<":"<<  ( (z%60) < 10 ? "0" : "" ) << (z%60)
                    << " " << ctempA<< " " << cwind << " "<<cSWradiationIn<<" "<<cLWradiationIn
		    <<" "<<h<<" "<<latentFlux<<" "<< bridgeTemp << " "<<frostDepth << " "<<Td << " ";

	    
            switch (condition)
            {
                case 1:
                   outfile<<"frosty";
                   break;
                case 2:
                   outfile<<"Icy/Snowy";
                   break;
                case 3:
                   outfile<< "Melting";
                   break;
                case 4:
                   outfile<<"Freezing";
                   break;
                case 5:
                   outfile<<"Wet";
                   break;
                case 0:
                   outfile<<"Dry";
                   break;
            }
            outfile<< endl;

            z++;// sends it through the minute loop again unless its at the next hour
         }
      }
      //turn the "new" data into the "old" data
      otempA = ntempA;
      oqa = nqa;
      owind = nwind;
      oSWradiationIn = nSWradiationIn;
      oLWradiationIn = nLWradiationIn;
      oprecipRate = nprecipRate;
      
      yearO = year;
      monthO = month;
      dayO = day;
      hourO = hour;
      
      linecounter = linecounter + 1;
      oldMinutes = minutes;
      z = oldMinutes;
      
      //strcpy(timestamp,time);
      
   } // end of bridge surface calculator loop. 
   outfile.close();        
   return 0;
}//end of main

float initialize( float bridgeTemp, float tempA, float h, float kbridge, float thickness, 
                  double thermalDiffusivity, int c, int i )
   { // this function is used to initialize the bridge layer temps before the forecasted temperatures are calculated.  
     //It uses the RWIS air and surface temperatures to force the lower level temperatures since only the surface 
     //temperature is actually measured.  It runs using the previous surface and air temperatures.
   
   static float nexttempB[19];// array to temporarily store the new values so they dont interfere with the following 
                              //calculations.
   static float tempB[19];  
   if ( i == 0 )// returns the calculated value to "bridgelevel" and skips any additional calculations
   {
      return tempB[c];
   }
   static float depth;
   depth = thickness / 19.0; // distance in meters between nodes
   static const float M = ((depth * depth) / (thermalDiffusivity * 60.0)); // unitless coefficient for heat transfer 
   //rates in concrete.    M = ((depth * depth) / (thermalDiffusivity * timestep));

   
   if ( i == 1)
   {
      tempB[c] = tempA ; // array to store the node temperatures.  Initially set to airTemp.
   }
   
   else
   {   
      nexttempB[0] = 1 / M * ( bridgeTemp + tempB[1] ) + ( 1 - ( 2 / M ) ) * tempB[0];
      if ( c > 0 && c < 18 )
      {
         nexttempB[c] = 1 / M * ( tempB[c-1] + tempB[c+1] ) + ( 1 - ( 2 / M ) ) * tempB[c];
      }
  
      nexttempB[18] =  (2 / M) * (h * depth / kbridge * tempA + tempB[17]) +
                 (1 - (2 / M) * (h * depth / kbridge + 1)) * tempB[18];
	 
      tempB[c] = nexttempB[c];
   }
   
   return tempB[c];
}

float precipitation( float precip, float h, float airTemp, float qa, float wind, float kbridge, float thick,
                     float &bridgeTemp, float &frostDepth, int &condition, int i )
{  // MM5 calculates precip rate in centimeters per hour, and the bridge will not hold all the water that falls on it.
   //This function divides the precip rate into kilograms per minute, calculates total precip. accumulation, and 
   //truncates the precip amount if more were to fall than the bridge will hold. The truncation will help keep 
   //unreasonable amounts of water from cooling the bridge through latent effects.  The evaporation rate and latent 
   //heat effects are calculated and the heat flux is returned to main.  The total frost depth is calculated.
   static const float waterDensity =  1000 ; //kg m-3
   static const float iceDensity = 917; // kg m-3
   static const float cp = 1007 ; // specific heat of air J/(kg k)
   static const float Le = 0.907; // Lewis # = thermal diffusivity of water / mass diffusivity  approximately constant.
                                  //Le = Lewis # raised to the 2/3 power.
   static const float R = 0.08314; // m3 bar/ kmol K  Universal gas constatant
   static const float M = 18.0; // Kg/Kmol  molecular weight of water.
   static const float epsilon = 0.622; // a constant.  Ratio of moleular weight of water/dry air
   static const float Lf = 3.34 * pow ( 10,5); // latent heat of freezing  J/kg  water
   static const float Lv = 2.5 * pow ( 10,6); // latent heat of condensation  J/kg  water
   static const float Cw = 4218 ;  // J/K kg  specific heat of water
   static const float concreteDensity = 2242.0;//kg/m3
   //static const float Kw = 0.613; // W/mK thermal conductivity of water
   static const float Cc = 880 ; // J/K kg specific heat of concrete
   static float thickness;// depth of the top concrete layer
   thickness = thick / 48.0;
   static float Mc;// mass of concrete in one node layer per unit area kg
   Mc = concreteDensity * thickness;  
   static float Tw; // temperature of the water already accumulated on the bridge
   static float precipDepth; // meters
   static float frozen; // depth of water that is frozen on the bridge
   static float totalDepth; // total depth of water accum.  Includes condensation/frost
   static float Mw; // per unit area  mass of existing accum.
   Mw = (precipDepth) * waterDensity ;
   float evapRate; // kg s-1 m-2
   float Psat;// pressure of water vapor at saturation
   float Patm;// pressure of water vapor in the surrounding air
   Patm = qa / epsilon;  // * 1000mb -- `in bars
   float airDensity; // density of air = p/(RT)  kg/m3  assuming pressure = 1000 mb
   airDensity = 100000/ ( 287 * airTemp);
   float hm; // mass transfer coefficient
   hm = h / (airDensity * cp * Le );
   float Mf; // per unit area   mass of water that just fell as precip
   Mf = precip * waterDensity ;
   float Te; // equilibrium temp of bridge and fallen precip
   float heatFlux; // latent heat flux lost or gained to the bridge.  Is returned to residual
   float freezeRate;  // kg/s of water freezing or melting.
   
   evapRate = 0; // makes sure evapRate is not carried over from last iteration
   precip = precip / 360000 ; // centimeters per hour to meters per minute

   if ( i == 0 ) // sets initial depth to zero
   {
      precipDepth = 0;
      frostDepth = 0;
   }
   
   precipDepth += precip;
   
   if ( precipDepth > 0.0025 ) // will need to refine this threshold later.  
   {
      precipDepth = 0.0025;
   }
   
   if ( precipDepth == 0 )
   {
      Tw = airTemp;
   }
   
   if (  precipDepth ==0 && bridgeTemp <= 273 ) // frost accumulation/evaporation possible
   {
      Psat = 6.1115 * exp( 17.67 * ( bridgeTemp - 273) / ( bridgeTemp - 29.5)) / 1000.0; // in bars
      //Psat = 6.1078 * exp( 21.875 * ( bridgeTemp - 273) / ( bridgeTemp - 7.66)) / 1000.0; // in bars
      //this is different from saturation over water! from Murray FW 1966 JAM 6
      evapRate = hm * M / R * ( Psat/bridgeTemp  - Patm / airTemp) ;  // kg/sec per unit area
      //frost will form if evapRate is negative -- water flux toward the bridge surface
      heatFlux =-evapRate * (Lv + Lf );
      frostDepth += ( -600 * evapRate / iceDensity); // depth after a minute of accumulation/evap
      if (frostDepth <= 0  )
      {
         frostDepth =0;
         heatFlux = 0;
      }
      
      totalDepth = frostDepth / 10.0;
      condition = 1;
      Tw = bridgeTemp;
   }
   
   if ( bridgeTemp > 273 && frostDepth>0 )
   {
      Tw = 273;
      condition = 3;
   }
   
   if ( precipDepth > 0 )
   {
      totalDepth = precipDepth + frostDepth / 10.0;
      frostDepth = 0;
   }

   if ( Tw > 273 &&  bridgeTemp > 273 ) // just evaporation or condensation
   {
      Te = ( Cw * (Mw * Tw  + Mf*airTemp) + Cc * Mc * bridgeTemp  ) / (Cw * (Mw + Mf) + Cc * Mc );
      Tw = bridgeTemp = Te;
      Psat = 6.112 * exp( 17.67 * ( bridgeTemp - 273) / ( bridgeTemp - 29.5)) / 1000.0; // in bars
      evapRate = hm * M / R * ( Psat / bridgeTemp - Patm / airTemp) ;  // kg/sec per unit area
      heatFlux = -evapRate * Lv;
   }
   
   if (totalDepth > 0)
   { 
      if ( Tw >= 273 && bridgeTemp <= 273 ) // freezing / condensation / evaporation
      {
         Te = ( Cw * (Mw * Tw  + Mf*airTemp) + Cc * Mc * bridgeTemp  ) / (Cw * (Mw + Mf) + Cc * Mc );
         Tw = bridgeTemp = Te;

         if ( Te >= 273 ) // no freezing -- evaporation or condensation
         {
            Psat = 6.112 * exp( 17.67 * ( bridgeTemp - 273) / ( bridgeTemp - 29.5)) / 1000.0; // in bars
            evapRate = hm * M / R * ( Psat / bridgeTemp - Patm / airTemp) ;  // kg/sec per unit area
            heatFlux = -evapRate * Lv;
            condition = 5;
         }
         else // freezing
         {
            heatFlux = kbridge * ( 273 - Te ) / (0.00762*60);
            freezeRate =  heatFlux / Lf;
            frozen = freezeRate/iceDensity * 60 + frozen;
            Tw = 273;
            condition = 4;
            if ( frozen > precipDepth )
            {
               condition = 2;
               if (frostDepth > 0.0 )
               {
                   condition = 1;
               }
               Tw = Te;
               frozen = totalDepth;
               heatFlux = 0;
            }
         }
      }
   
      else if ( Tw <= 273 && bridgeTemp <= 273 ) // snow on bridge  no evap, no melting
      {
         if ( precipDepth > 0.00063 ) // only 0.025 in water (~.25 in. of snow) allowed because of plowing
         {
            precipDepth = 0.00063;
         }
         if (frostDepth == 0)
	 {
	    heatFlux = 0 ;
	 }
         condition = 2;
         if (frostDepth > 0.0 )
         {
             condition = 1;
         }

      }
   
      else if ( Tw <= 273 && bridgeTemp > 273 ) // melting possible
      {
         Te = ( Cw * (Mw * Tw  + Mf*airTemp) + Cc * Mc * bridgeTemp  ) / (Cw * (Mw + Mf) + Cc * Mc );
         Tw = bridgeTemp = Te;

         if ( Te <= 273 ) // no melting
         {
            heatFlux = 0;
            condition = 2;
            if (frostDepth > 0.0 )
            {
                condition = 1;
            }

         }
         else // melting - heat taken from bridge
         {
            heatFlux = kbridge * ( 273 - Te ) / (0.00762*60);
            freezeRate =  heatFlux / Lf;
            frozen = freezeRate/waterDensity * 60 + frozen;
	    frostDepth = freezeRate/waterDensity * 60 + frostDepth;
            Tw = 273;
            condition = 3;
            if ( frozen < 0 ||  frostDepth <0)
            {
               Tw = Te;
               frozen = 0;
	       frostDepth =0;
               heatFlux = 0;
               condition = 5;      
            }
	    Tw = 273;
	    
         }
      }
   }

   // if statements make sure evaporation rate does not exceed totalDepth, the amount on the bridge

   if ( totalDepth > (evapRate * 60 / waterDensity) && frostDepth ==0)
   {
      totalDepth = totalDepth - (evapRate * 60 / waterDensity);
      precipDepth = precipDepth - (evapRate * 60 / waterDensity);
   }
   else if ( totalDepth <= 0 && precipDepth <= 0 )
   {
      heatFlux = 0;
      totalDepth = 0;
      precipDepth = 0;
      condition = 0;
   }
   else if (frostDepth == 0)
   {
      heatFlux = -1 * totalDepth * waterDensity * Lv / 60.0;
      totalDepth = 0;
      condition = 0;
   }
   return heatFlux;
}  

float radiative(float SW, float LW, float tempSfc, float absorptivitySW, float absorptivityLW, float airtemp)
{ 
    // This function is used to find the radiative energy flux available to influence the bridge temperature 
  
    static const float stephanBoltz = 5.67 * pow(10,-8); //(W m-2 K-4)
    float radiation; //leftover energy flux ( W m-2) received from atmosphere 
    
   // the following tries to improve the model LW  
   // Get rid of this if your input has good Longwave radiation values!!! 
LW = LW*1.11;
if(LW > (stephanBoltz * pow(airtemp,4.0)))
{
   LW =  (stephanBoltz * pow(airtemp,4.0));
}


    radiation = ( absorptivityLW  * LW ) + ( absorptivitySW  * SW ) -
               absorptivityLW * stephanBoltz * pow(tempSfc,4.0) ; // incoming radiation minus emitted radiation
	       
    return radiation;
}
float convection( float wind )
{  // this function finds the right convection heat transfer coefficient (in Wm-2) given windspeed.  It will be the 
   //average coefficient over a 10 meter span of bridge. h is computed for "tripped" flow, where turbulence begins 
   //at the very edge of the bridge. Immediately turbulent flow is a good assumption because of bridge 
   // rails and other bridge complexities.

   float h, reynoldsNumber;//, Cf, xCrit;
   static const float kair = 2.4 * pow(10,-2);// conductivity of air Wm-1k-1
   static const float length = 10 ; // meters.  length of bridge exposed to air flow
   static const float prandtlNumber = 0.71;// for air
   static const float kinematicViscosity = 1.4 * pow(10,-5) ;//m2 s-1  for air
   
   if (wind == 0)
   {
      wind =0.01; // keeps the reynolds number from going to zero.
   }  
   
   reynoldsNumber = wind * length / kinematicViscosity;

   h = kair / length *0.037* pow(reynoldsNumber,(0.8)) * pow(prandtlNumber,(1/3.0));

   return h;
}   

float bridgelevels (float radiation, float latentHeat, float tempSfc, float tempAir, float h, float thickness, 
                    float kbridge, double thermalDiffusivity,int l, int t)
{   /* n is used for node definition, it gives the vertical position of the node.  At n = 0,  is a node representing 
    the top surface of the bridge.  Since the bridge is much longer and wider than deep, a temperature gradient is 
    assumed to exist only in the verical. Max n = 19, total depth is determined by the property file.  
    The second array dimension is t, the time in minutes. The temperature of any node at any time is stored 
    in the two dimensional array, temperature. Nodes are initialized with  the  temperatures found in "initialize".
    
    static const float kbridge =1.4;//W/mK. varies 0.8129 to 1.401. heat conductivity of concrete
    thermalDiffusivity = 6.9677 * pow(10,-7) = k/heat cap. m2 s-1. concrete. varies 4.9032e-7 to 6.9677e-7
    */
    static float length = thickness / 20.0;// m  distance between nodes
    static float M = ((length * length) / (thermalDiffusivity * 60.0));//heat transfer rate coefficient, timestep 
                                                                       //= 60 seconds
    float dummy[6] = {0};// useless floats to pass to function "initialize"
    static float temperature[20][3000];

    if ( (1+(h*length/kbridge)) / M > 0.5 )
    {
       cerr <<"thickness, h, and thermal diffusivity combination creates numerical instability. Must be less than 0.5.  " 
            << (1+(h*length/kbridge)) / M << endl;
       exit( 32 );
    }
    if ( t > 3000 )
    {
       cerr << " file will be shortened to 3000 minutes " << endl;
       exit( 31 );
    }

    // sets initial bridge node temperatures

    if ( l == 0  )
    {
       int i = 0;
       temperature[0][t] = tempSfc;
       
       for ( int c = 1; c < 20; c++)
       {
          temperature[c][t] = initialize( dummy[0], dummy[1], dummy[2], dummy[3], dummy[4], dummy[5], (c-1), i );

	  if ( temperature[c][t] == 0 )
	  {
	     cerr<< " intialization = 0. insufficient rwis initialization" << endl;
	     exit(20);
	  }
       }
    }
    // for any other time, t
    // for top level, there is conduction, convection, and residual (radiative, and latent heat) processes at work
    else
    {  
       temperature[0][t-1] = tempSfc;//makes sure that the surface temperature from the latentHeat function is used
       temperature[0][t] = (2 / M) * ( h  * length / kbridge * tempAir + ( (radiation + latentHeat) *
                           length / kbridge ) + temperature[1][t-1]) + ( 1 - (2 / M) *
                           ( h * length / kbridge + 1 )) * temperature[0][t-1];
       
       // for lower levels: just conduction for middle layers, and convection for the lowest level.
       for (int n = 1; n < 19; n++ )
       { 
          temperature[n][t] = ( 1 / M * (temperature[n-1][t-1] + temperature[n+1][t-1] )
                                + (1 - 2 /M ) * temperature[n][t-1]);
      
       }
       
       temperature[19][t] =  (2 / M) * (h * length / kbridge * tempAir +
                              temperature[18][t-1]) + (1 - (2 / M) * 
                              (h * length / kbridge + 1)) * temperature[19][t-1];	      
	      
       if ( t > 0 && fabs(temperature[0][t-1] - temperature[0][t]) > 4 )
       {
          cerr << "  heating/cooling rate too high.  check inputs  "<< 
               temperature[0][t-1] <<"  "<<temperature[0][t] << "  "<<t<<endl;
          exit(25);
       }
    }
    
    if ( temperature[0][t] < 200 )
    {
       cerr<< "bad calculated surface temperature: check inputs   "<<temperature[0][t]<<"  time: "<<t<<endl;
       exit(10);
    }

    return temperature[0][t];
}

