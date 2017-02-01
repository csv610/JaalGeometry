#include  <boost/container/flat_map.hpp>

#include <vector>
#include <map>
#include <iostream>

using namespace std;

std::pair<int,int> max_frequency(const vector<int>  &items)
{
   map<int, unsigned>  freqCount;
   for( int v: items) {
        if( freqCount.find(v) == freqCount.end() ) freqCount[v] = 0;
        freqCount[v]++;
   }
   unsigned maxval = 0;
   for( auto it = freqCount.begin(); it != freqCount.end(); ++it)
        maxval = max( maxval, it->second );

   for( auto it = freqCount.begin(); it != freqCount.end(); ++it) {
        if( it->second == maxval) return std::make_pair(it->first, maxval);
   }
}

int main()
{
     vector<int> a(4);
     a[0] = 1;
     a[1] = 2;
     a[2] = 3;
     a[3] = 4;
     pair<int,int> maxfreq = max_frequency(a);
     cout << maxfreq.first << " " << maxfreq.second << endl;
    


}
