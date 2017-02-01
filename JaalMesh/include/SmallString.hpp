#pragma once

class JSmallString {
public:
    JSmallString(const char* p);
    JSmallString(const String& s);

    JSmallString(String&& s);
    ~JSmallString();
    size_t size() const {
        return _size;
    }

    char* begin() {
        return _size < 16 ? str : _beg;
    }
    char* end() {
        return begin() + _size;
    }

private:
    unsigned _size;
    union {
        char* _beg;
        char str[16];
    };
};

inline JSmallString::JSmallString(const char* p) : _size(strlen(p))
{
    if (_size < 16) {
        // Small string optimization: store in static str
        memcpy(str,p,_size+1);
    } else {
        // Dynamically allocate string (like before)
        _beg(new char[_size+1]);
        memcpy(_beg,p,_size+1);
    }
}

inline JSmallString::JSmallString(const JSmallString& s) : _size(s.length())
{
    if (_size < 16) {
        // Small string optimization: store in static str
        memcpy(str,s.str,_size+1);
    } else {
        // Dynamically allocate string (like before)
        _beg(new char[_size+1]);
        memcpy(_beg,s._beg,_size+1);
    }
}


inline JSmallString::SString(JSmallString&& s) : _size(s.length())
{
    if (_size < 16) {
        memcpy(str,s.str,_size+1);
    } else {
        _beg = s._beg; // steal resource from s
        s._beg = nullptr; // prevent s from deleting on destruction
    }
}

inline JSmallString::~JSmallString()
{
    if(_size > 16) delete [] _beg;
}
}
