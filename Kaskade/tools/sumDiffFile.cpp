/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <string>
#include <algorithm>
#include <cmath>
#include <time.h>

std::string numString(std::string const buffer, int* afterPos, int *ntyp)
/*
  ntyp: 1=int, 2=fixed (float), 3=scientific(double), 4=time, 5=percent
*/
{
  std::string const numerical="0123456789+-.eE";
  std::string result="";
  int i=0, strSize=buffer.size(), resultSize=0; 
  *afterPos=0; *ntyp=0;
  bool inNumber=false, hasNumber=false;
  while ( i<strSize )
  {
    if ( inNumber && ( ( buffer.substr(i,2) == "s " ) || ( buffer.substr(i,2) == "s\0" ) ) )
    {
      *ntyp=4;
      *afterPos=i+2;
      return result;
    };
    if ( inNumber && ( buffer[i]=='%' ) )
    {
      *ntyp=5;
      *afterPos=i+1;
      return result;
    };
    if ( inNumber && ( buffer[i] == 'e' || buffer[i] == 'E' ) )
    {
      if ( *ntyp!=3 )
        *ntyp=3;
      else
      {
        *ntyp=0;
        *afterPos=i+1;
        return "";
      }
    };
    if ( numerical.find(buffer[i]) != std::string::npos )
    {
      if ( ( buffer[i] == 'e' || buffer[i] == 'E' || buffer[i] == 's'  || buffer[i] == '.' ) && ! inNumber  )
      {
        i++; continue;
      }
      else
      {
        if ( numerical.find(buffer[i])<10 )
        {
          inNumber=true;
          hasNumber=true;
          if ( *ntyp==0 && result.size()<2 ) // multiple leading +- signs are no valid number
            *ntyp=1;
        };
        if ( buffer[i]=='.' )
        {
          if ( *ntyp!=2 )
          {
            *ntyp=2;
          }
          else  // multiple dots in number-string - no valid number
          {
            *ntyp=0;
            *afterPos=i+1;
            return "";
          }
        };
        result += buffer[i];
      }
    }
    else
    {
      inNumber=false;
      if ( hasNumber )
      {
        *afterPos=i;
        resultSize=result.size();
        if ( resultSize==1 && numerical.find(result)>9 )  // a single + - e E . is no number
        {
          *ntyp=0;
          return "";
        };
        if ( result[resultSize-1]=='.' )  // no digits following . -> number is treated as integer
        {
          *ntyp=1;
          return result.substr(0,resultSize-1);
        };
        if ( result[resultSize-1]=='e' || result[resultSize-1]=='E' ) // no exponent following - ignore e/E
        {
          *ntyp=2;
          return result.substr(0,resultSize-1);
        };
        return result;
      }
    };
    i++;
  };
  *afterPos=i;
  resultSize=result.size();
  if ( resultSize==1 && numerical.find(result)>9 )
  {
    *ntyp=0;
    return "";
  };
  if ( resultSize>1 )
    if ( result[resultSize-1]=='.' )
    {
      *ntyp=1;
      return result.substr(0,resultSize-1);
    };
  return result;
}

char typeOfLine(std::string const buffer)
{
  std::string const modes="cda", digit="0123456789";
  int charPos;
  if ( digit.find(buffer[0]) != std::string::npos )
  {
    int const bufLen=buffer.size();
    for (int i=0;i<bufLen;i++)
    {
      charPos=modes.find(buffer[i]);
      if ( charPos!=std::string::npos )
        return modes[charPos];
    }
  };
  return ' ';
}

int main(int const argc, char* const argv[])
{
  if ( argc != 2 && argc != 4 )
  {
    std::cout << "usage: sumDiffFile [ -tol floatvalue ] diffFile" << std::endl;
    exit(10);
  };
  
  int const maxInt=100, maxFloat=100, maxDouble=1000, maxTime=100;
  double const intThresh=1, floatThresh=1.0e-7, timeThresh=0.001, doubleThresh=1.0e-16;
  int line=0, pos=0, mode=0, dline, dpos, npos, numberType, bufferLength, line1c=0, line2c=0,
      int1c=0, int2c=0, float1c=0, float2c=0, double1c=0, double2c=0,
      time1c=0, time2c=0, intTotC=0, floatTotC=0, doubleTotC=0, timeTotC=0,
      intMaximum=0, intAbsDiff, excludedBlocks=0,
      intMaxLine=0, floatMaxLine=0, doubleMaxLine=0, timeMaxLine=0,
      intMaxPos=0, floatMaxPos=0, doubleMaxPos=0, timeMaxPos=0,
      intRelMaxLine=0, floatRelMaxLine=0, doubleRelMaxLine=0, timeRelMaxLine=0,
      intRelMaxPos=0, floatRelMaxPos=0, doubleRelMaxPos=0, timeRelMaxPos=0;
  double doubleMaximum=0.0, floatMaximum=0.0, timeMaximum=0.0,
         doubleRelMaximum=0.0, intRelMaximum=0.0, floatRelMaximum=0.0, timeRelMaximum=0.0,
         floatAbsDiff, timeAbsDiff, doubleAbsDiff, intRelDiff, floatRelDiff, timeRelDiff, doubleRelDiff,
         relTolerance=0.01;
  bool overflow, mismatch=false, eofin=false, exclude=false;
  char changeMode;
  std::array<int,maxInt> int1, int2, intLine, intPos;
  std::array<float,maxFloat> float1, float2, floatLine, floatPos;
  std::array<float,maxTime> time1, time2, timeLine, timePos;
  std::array<double,maxDouble> double1, double2, doubleLine, doublePos;
  std::string buffer, numberString;
  std::stringstream numberStream;
  std::fstream diffStream;
  std::string filename;
  time_t t;
  filename=argv[1];
  if ( filename == "-tol" )
  {
    filename=argv[3];
    numberStream.str(argv[2]);
    numberStream >> relTolerance;
    if ( numberStream.fail() )
    {
      std::cout << "*** Number format error on command line ***" << std::endl;
      exit (9);
    }
  };
  diffStream.open(filename,std::fstream::in);
  if ( ! diffStream.is_open() )
  {
    std::cout << "*** Fatal error: open of input-stream failed ***" << std::endl;
    exit(10);
  };

  time(&t);
  std::cout << "*** Results summary of " << ctime(&t);

  while ( ! diffStream.eof() )
  {
    std::getline(diffStream,buffer);
    eofin=diffStream.eof();
    bufferLength=buffer.size();
    line++;
    pos=0;
    dpos=1;
    changeMode=typeOfLine(buffer);
    while ( changeMode == 'd' || changeMode=='a' )
    {
      if ( changeMode=='d' )
      {
        std::cout << "*** deleted lines not summarized at line " << line << std::endl;
        while (true)
        {
          std::getline(diffStream,buffer);
          line++;
          if ( buffer[0]=='<' )
            std::cout << buffer << std::endl;
          else
            break;
        }
      };
      if ( changeMode=='a' )
      {
        std::cout << "*** added lines not summarized at line " << line << std::endl;
        while (true)
        {
          std::getline(diffStream,buffer);
          line++;
          if ( buffer[0]=='>' )
            std::cout << buffer << std::endl;
          else
            break;
        }
      };
      excludedBlocks++;
      eofin=diffStream.eof();
      changeMode=typeOfLine(buffer);
    }
    if ( changeMode=='c' || eofin )
    {
      if ( int1c==int2c )
      { 
        for (int j=0;j<int1c;j++)
        {
          intTotC++;
          intAbsDiff=std::abs(int1[j]-int2[j]);
          intRelDiff=((float) intAbsDiff)/std::max(0.5*(std::abs(int1[j])+std::abs(int2[j])),intThresh);
          if (intAbsDiff>intMaximum)
          {
            intMaximum=intAbsDiff;
            intMaxLine=intLine[j];
            intMaxPos=intPos[j];
          };
          if (intRelDiff>intRelMaximum)
          {
            intRelMaximum=intRelDiff;
            intRelMaxLine=intLine[j];
            intRelMaxPos=intPos[j];
          };
        }
      }
      else
      {
        std::cout << "*** Number of integer numbers to compare don't match, count1=" << int1c <<
          ",  count2=" << int2c << std::endl << "*** integer numbers excluded from compare in lines " <<
           line1c+1 << " through " << line2c << std::endl;
        exclude=true;
        
      };
      if ( time1c==time2c )
      {
        for (int j=0;j<time1c;j++)
        {
          timeTotC++;
          timeAbsDiff=std::fabs(time1[j]-time2[j]);
          timeRelDiff=timeAbsDiff/std::max(0.5*(std::fabs(time1[j])+std::fabs(time2[j])),timeThresh);
          if (timeAbsDiff>timeMaximum)
          {
            timeMaximum=timeAbsDiff;
            timeMaxLine=timeLine[j];
            timeMaxPos=timePos[j];
          };
          if (timeRelDiff>timeRelMaximum)
          {
            timeRelMaximum=timeRelDiff;
            timeRelMaxLine=timeLine[j];
            timeRelMaxPos=timePos[j];
          };
        }
      }
      else
      {
        std::cout << "*** Number of time numbers to compare don't match, count1=" << time1c <<
          ",  count2=" << time2c << std::endl << "*** time numbers excluded from compare in lines " <<
           line1c+1 << " through " << line2c << std::endl;
        exclude=true;
      };
      if ( float1c==float2c )
      {
        for (int j=0;j<float1c;j++)
        {
          floatTotC++;
          floatAbsDiff=std::fabs(float1[j]-float2[j]);
          floatRelDiff=floatAbsDiff/std::max(0.5*(std::fabs(float1[j])+std::fabs(float2[j])),floatThresh);
          if (floatAbsDiff>floatMaximum)
          {
            floatMaximum=floatAbsDiff;
            floatMaxLine=floatLine[j];
            floatMaxPos=floatPos[j];
          };
          if (floatRelDiff>floatRelMaximum)
          {
            floatRelMaximum=floatRelDiff;
            floatRelMaxLine=floatLine[j];
            floatRelMaxPos=floatPos[j];
          }
        };
      }
      else
      {
        std::cout << "*** Number of fixed point numbers to compare don't match, count1=" << float1c <<
          ",  count2=" << float2c << std::endl << "*** fixed point numbers excluded from compare in lines " <<
           line1c+1 << " through " << line2c << std::endl;
        exclude=true;
      };
      if ( double1c==double2c )
      {
        for (int j=0;j<double1c;j++)
        {
          doubleTotC++;
          doubleAbsDiff=std::fabs(double1[j]-double2[j]);
          doubleRelDiff=doubleAbsDiff/std::max(0.5*(std::fabs(double1[j])+std::fabs(double2[j])),doubleThresh);
          if (doubleAbsDiff>doubleMaximum)
          {
            doubleMaximum=doubleAbsDiff;
            doubleMaxLine=doubleLine[j];
            doubleMaxPos=doublePos[j];
          };
          if (doubleRelDiff>doubleRelMaximum)
          {
            doubleRelMaximum=doubleRelDiff;
            doubleRelMaxLine=doubleLine[j];
            doubleRelMaxPos=doublePos[j];
          };
        }
      }
      else
      {
        std::cout << "*** Number of floating point numbers to compare don't match, count1=" << double1c <<
          ",  count2=" << double2c << std::endl << "*** floating point numbers excluded from compare in lines " <<
           line1c+1 << " through " << line2c << std::endl;
        exclude=true;
      };
      if ( exclude )
        excludedBlocks++;
      if ( eofin )
        break;
      dline=0;  
      int1c=0; int2c=0;
      float1c=0; float2c=0; 
      double1c=0; double2c=0;
      time1c=0; time2c=0;
      overflow=false;
      mismatch=false;
      exclude=false;
      line1c=line;
    }
    else
    {
      if ( buffer[0]=='<' )
        mode=1; 
      else if ( buffer[0]=='>' )
        mode=2;
      else if ( buffer[0]=='-' )
        mode=0;  dline=0;
      if ( mode==1 || mode==2 )
      {
        line2c=line;
        while (dpos<bufferLength)
        {
          numberString=numString(buffer.substr(dpos),&npos,&numberType);
          pos++;
          numberStream.str(numberString);
          numberStream.seekg(0);
          dpos += npos;
          if ( numberType==0 )
            break;
          switch (numberType)
          {
            case 1:
            {
              if ( mode==1 )
              {
                overflow=(int1c==maxInt);
                if (overflow)
                  break;
                numberStream >> int1[int1c];
                intLine[int1c]=line;
                intPos[int1c]=pos;
                int1c++;  
              }
              else
              {
                overflow=(int2c==maxInt);
                if (overflow)
                  break;
                numberStream >> int2[int2c];
                int2c++;  
              };
            };
            break;
            case 2:
            {
            {
              if ( mode==1 )
              {
                overflow=(float1c==maxFloat);
                if (overflow)
                  break;
                numberStream >> float1[float1c];
                floatLine[float1c]=line;
                floatPos[float1c]=pos;
                float1c++;
              }
              else
              {
                overflow=(float2c==maxFloat);
                if (overflow)
                  break;
                numberStream >> float2[float2c];
               float2c++;
              };
            };
            break;
            case 3:
            {
              if ( mode==1 )
              {
                overflow=(double1c==maxDouble);
                if (overflow)
                  break;
                numberStream >> double1[double1c];
                doubleLine[double1c]=line;
                doublePos[double1c]=pos;
                double1c++;
              }
              else
              {
                overflow=(double2c==maxDouble);
                if (overflow)
                  break;
                numberStream >> double2[double2c];
                double2c++;
              };
            };
            break;
            case 4:
            {
              if ( mode==1 )
              {
                overflow=(time1c==maxTime);
                if (overflow)
                  break;
                numberStream >> time1[time1c];
                timeLine[time1c]=line;
                timePos[time1c]=pos;
                time1c++;
              }
              else
              {
                overflow=(time2c==maxTime);
                if (overflow)
                  break;
                numberStream >> time2[time2c];
                time2c++;
              };
            };
            break;
            case 5: /* evaluation of percent type numbers not yet implemented! */
            break;
            default:
            {
               std::cout << "*** Internal program error - invalid numberType=" << 
                            numberType << ", mode=" << mode << std::endl;
               exit (10);
            };
            break;
          }
          if (overflow)
          {
            std::cout << "*** Array overflow in line=" << line <<", numberType="
                      << numberType << std::endl;
            exit(9);
          }
        }
      }
     }
    }
  };
  if ( excludedBlocks!=0 )
    std::cout << std::endl;
  std::cout << "*** Summary of absolute deviations ***" << std::endl;
  if ( intMaxLine!=0 )
    std::cout << "Maximum integer        deviation is " << intMaximum << " at line " << intMaxLine  <<
                 " and position " << intMaxPos << std::endl;
  if ( floatMaxLine!=0 )
    std::cout << "Maximum fixed point    deviation is " << floatMaximum << " at line " << floatMaxLine  <<
                 " and position " << floatMaxPos << std::endl;
  if ( doubleMaxLine!=0 )
    std::cout << "Maximum floating point deviation is " << doubleMaximum << " at line " << doubleMaxLine  <<
                 " and position " << doubleMaxPos << std::endl;
  if ( timeMaxLine!=0 )
    std::cout << "Maximum time           deviation is " << timeMaximum << " at line " << timeMaxLine  <<
                 " and position " << timeMaxPos << std::endl << std::endl;
  std::cout << "*** Summary of relative deviations ***" << std::endl;
  if ( intRelMaxLine!=0 )
    std::cout << "Maximum integer        deviation is " << intRelMaximum << " at line " << intRelMaxLine  <<
                 " and position " << intRelMaxPos << std::endl;
  if ( floatRelMaxLine!=0 )
    std::cout << "Maximum fixed point    deviation is " << floatRelMaximum << " at line " << floatRelMaxLine  <<
                 " and position " << floatRelMaxPos << std::endl;
  if ( doubleRelMaxLine!=0 )
    std::cout << "Maximum floating point deviation is " << doubleRelMaximum << " at line " << doubleRelMaxLine  <<
                 " and position " << doubleRelMaxPos << std::endl;
  if ( timeRelMaxLine!=0 )
    std::cout << "Maximum time           deviation is " << timeRelMaximum << " at line " << timeRelMaxLine  <<
                 " and position " << timeRelMaxPos << std::endl;
  std::cout << "*** Summary of comparations ***" << std::endl;
  std::cout << "The number of integer numbers        compared is " << intTotC << std::endl <<
               "The number of fixed point numbers    compared is " << floatTotC << std::endl <<
               "The number of floating point numbers compared is " << doubleTotC << std::endl <<
               "The number of time numbers           compared is " << timeTotC << std::endl;
  if ( intMaxLine==0 && floatMaxLine==0 && doubleMaxLine==0 && timeMaxLine==0 )
    std::cout << "*** No deviations found          ***" << std::endl;
  if ( std::max(intRelMaximum,std::max(floatRelMaximum,doubleRelMaximum))>relTolerance )
  {
    std::cout << "*** Warning: Maximum relative error of integer, fixed point and floating point is larger than required precision "
              << relTolerance << "***" << std::endl;
    exit (1);
  };
  if ( excludedBlocks!=0 )
  {
    std::cout << std::endl << 
      "*** Warning: " << excludedBlocks << " compare blocks have been (partially) excluded from summary!" << std::endl;
    exit (2);
  }
}
