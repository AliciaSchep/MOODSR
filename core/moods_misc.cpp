#include "moods.h"
#include "moods_misc.h"
#include <iostream>


using std::vector;
using std::size_t;

namespace MOODS { namespace misc{

    // basically base-2 logarithm of a, rounded up
    unsigned int shift(unsigned int a)
    {
        unsigned int s = 0;
        unsigned int b = 1;
        
        while (b < a){
            s += 1;
            b = b << 1;
        }
        return s;
    }
        
    bits_t mask(unsigned int a){
        bits_t b = 1;
        
        while (b < a){
            b = b << 1;
        }
        return b-1;
    }

}}