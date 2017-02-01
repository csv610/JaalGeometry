#pragma once

#include <sys/time.h>

class RegisterCounter
{
public:
    typedef unsigned long long int   register_counter_t;
    RegisterCounter() {
        rtotal = 0;
    }

    void  reset() {
        rtotal  = 0;
    }
    void  start() {
        rbegin = rdtsc();
    }

    void  stop()
    {
        rend  = rdtsc();
        rtotal +=  rend-rbegin;
    }
    register_counter_t  getCount() {
        return rtotal;
    }

private:
    register_counter_t  rtotal, rbegin, rend;
    register_counter_t rdtsc()
    {
        register_counter_t x;
        __asm__ volatile(".byte 0x0f,0x31" : "=A" (x));
        return x;
    }
};

class JStopWatch
{
public:
    JStopWatch() {
        rtotal = 0;
    }

    void  reset() {
        rtotal  = 0;
    }

    void  start() {
        gettimeofday(&tp, 0);
        rbegin = (double)tp.tv_sec + (1.E-06)*tp.tv_usec;
    }

    void  stop()
    {
        gettimeofday(&tp, 0);
        rend  = (double)tp.tv_sec + (1.E-06)*tp.tv_usec;
        rtotal +=  rend-rbegin;
    }

    double getSeconds() {
        return rtotal;
    }

private:
    struct  timeval tp;
    double  rtotal, rbegin, rend;
};

