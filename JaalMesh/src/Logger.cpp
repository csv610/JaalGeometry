#include "Logger.hpp"
/*
using namespace log4cxx;
using namespace log4cxx::helpers;
LoggerPtr JLogger::logger = 0;
*/
using namespace std;

JLogger* JLogger::instance = nullptr;

///////////////////////////////////////////////////////////////////////////////

JLogger* JLogger ::getInstance()
{
    if( instance == nullptr)
        instance = new JLogger;
    return instance;

    /*
      if(logger == 0) {
         cout << "Creating basic log configuration" << endl;
         BasicConfigurator::configure();
         logger = Logger::getLogger("JMeshLogger");
     }
      infile.open(filename.c_str(), ios::out);
    */
}

void JLogger ::setFileName( const string &s)
{
    filename = s;
    if( infile.is_open() )
        infile.close();
    infile.open(filename.c_str());
}

void JLogger ::setInfo( const string &mesg)
{
    if( infile.is_open() )
        infile  << "Info: " << mesg << endl;;
}

void JLogger ::setWarn( const string &mesg)
{
    if( infile.is_open() )
        infile  << "Warn: " << mesg << endl;
}

void JLogger ::setError( const string &mesg)
{
    if( infile.is_open() )
        infile  << "Error: " << mesg << endl;
}
