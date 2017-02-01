#ifndef SIMP_ARR_H
#define SIMP_ARR_H

#include <stdlib.h>

#include <iostream>
using namespace std;

template <typename T>
class SimpleArray
{
private:
    T* arr;
    int arrSize;
    int arrAllocated;

public:
    SimpleArray() : arr(0) , arrSize(0), arrAllocated(0) {}

    SimpleArray( unsigned s ) :arrSize(s), arrAllocated(s)
    {
        resize(s);
    }

    ~SimpleArray() {
        clear();
    }

    T**  ptr()            {
        return &arr;
    }
    int& size()           {
        return arrSize;
    }
    int  size()     const {
        return arrSize;
    }
    int& capacity()       {
        return arrAllocated;
    }
    int  capacity() const {
        return arrAllocated;
    }

    int  resize( unsigned s )
    {
        if( s == 0) return 0;
        clear();
        arr = (T*)malloc(s*sizeof(T));
        for (unsigned i = 0; i < s; ++i)
            new (arr+i) T();
        arrSize   = s;
        arrAllocated = s;
        return 0;
    }

    void clear()
    {
        if( arr == 0) return;
        for (int i = 0; i < size(); ++i)
            arr[i].~T();
        free(arr);
        arr = 0;
        arrSize = 0;
        arrAllocated = 0;
    }

    typedef T* iterator;
    typedef const T* const_iterator;
    iterator       begin()       {
        return arr;
    }
    const_iterator begin() const {
        return arr;
    }
    iterator         end()       {
        return arr + arrSize;
    }
    const_iterator   end() const {
        return arr + arrSize;
    }


    T& operator[]( unsigned idx )       {
        return arr[idx];
    }
    T  operator[]( unsigned idx ) const {
        return arr[idx];
    }
};

#define ARRAY_INOUT( A ) A.ptr(), &A.capacity(), &A.size()
#define ARRAY_IN( A ) &A[0], A.size()

#endif
