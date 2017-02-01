#pragma once

#include <string>
#include <vector>
#include <boost/any.hpp>

typedef std::string  AttribString;

class JAttribute {
public:
    static const int DENSE_ATTRIBUTE  = 0;
    static const int SPARSE_ATTRIBUTE = 1;

    static bool isChar(const boost::any & operand);
    static bool isFloat(const boost::any & operand);
    static bool isDouble(const boost::any & operand);
    static bool isInteger(const boost::any & operand);
    static bool isString(const boost::any & operand);

    JAttribute() {}

    JAttribute(const AttribString &s, const boost::any & p)
    {
        name = s;
        value = p;
    }
    ~JAttribute()
    {
    }

    template <typename T>
    JAttribute(T n)
    {
        value = n;
    }
    boost::any value;
    AttribString  name;
};

inline bool JAttribute :: isChar( const boost::any &operand)
{
    return operand.type() == typeid(char);
}

inline bool JAttribute :: isFloat( const boost::any &operand)
{
    return operand.type() == typeid(float);
}

inline bool JAttribute :: isDouble( const boost::any &operand)
{
    return operand.type() == typeid(double);
}

inline bool JAttribute :: isInteger( const boost::any &operand)
{
    return operand.type() == typeid(int);
}

inline bool JAttribute :: isString( const boost::any &operand)
{
    return operand.type() == typeid(std::string);
}

////////////////////////////////////////////////////////////////////////////////
class JAttributeManager 
{
public:
    typedef boost::shared_ptr<JAttribute> JAttributePtr;

    ~JAttributeManager()
    {
        attributes.clear();
    }

    void clearAll()
    {
        attributes.clear();
    }

    template<class T>
    int setAttribute(const AttribString &s, const T &val)
    {
        int nAttribs = attributes.size();
        for( int i = 0; i < nAttribs; i++) {
            if( attributes[i]->name == s ) {
                attributes[i]->value = val;
                return 0;
            }
        }
        JAttributePtr a(new JAttribute(s,val));
        attributes.push_back( a );
        return 0;
    }


    int getAttribute(const AttribString &s, boost::any &attrib) {
        int nAttribs = attributes.size();
        for( int i = 0; i < nAttribs; i++) {
            if( attributes[i]->name == s ) attrib = attributes[i];
            return 0;
        }
        return 1;
    }

    template<class T>
    int getAttribute(const AttribString &s, T &val) const
    {
        int nAttribs = attributes.size();
        for( int i = 0; i < nAttribs; i++) {
            if( attributes[i]->name == s ) {
                val = boost::any_cast<T>(attributes[i]->value);
                return 0;
            }
        }
        return 1;
    }

    bool hasAttribute(const AttribString &s) const
    {
        int nAttribs = attributes.size();
        for( int i = 0; i < nAttribs; i++)
            if( attributes[i]->name == s ) return 1;
        return 0;
    }

    void deleteAttribute( const AttribString &s)
    {
        int nSize = attributes.size();
        for( int i = 0; i < nSize; i++) {
            if( attributes[i]->name == s) {
                attributes.erase(attributes.begin() + i );
                return;
            }
        }
    }

    int getNumAttributes() const
    {
        return attributes.size();
    }

    int getAttributeNames( vector<AttribString> &names) const
    {
        names.clear();
        int nSize = attributes.size();
        for( int i = 0; i < nSize; i++)
            names.push_back( attributes[i]->name );
        return 0;
    }

    std::vector<JAttributePtr> attributes;
};

typedef boost::shared_ptr<JAttributeManager> JAttributeManagerPtr;

