#pragma once

#include <values.h>

#ifdef HAVE_LOG4CXX
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/helpers/exception.h>
using namespace log4cxx;
using namespace log4cxx::helpers;
#endif

#include <string>
#include <fstream>

#define JBreak() \
{       \
   cout << "Break in file " << __FILE__ << " at " << __LINE__ << endl; \
   getchar(); \
} \
 
#define JExit() \
{       \
   cout << "Exit in file " << __FILE__ << " at " << __LINE__ << endl; \
   exit(0); \
} \
 
#define JNoImpl() \
{  \
    cout << "Not yet implemented "  << __FILE__ << " at " << __LINE__ << endl; \
    exit(0);\
}\
 
using namespace std;

class JLogger {
//  static log4cxx::LoggerPtr logger;

    static JLogger* instance;
public:
    static JLogger* getInstance();

    void setFileName( const string &s);
    void setTrace(const string &);
    void setInfo(const string  &);
    void setWarn( const string &);
    void setError( const string &);
private:
    JLogger() {}
    JLogger(JLogger const&) {};     // copy constructor is private
    JLogger& operator=(JLogger const&) {}; // assignment operator is private
    string filename = "logfile.txt";
    ofstream infile;
};
