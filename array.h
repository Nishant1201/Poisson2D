#ifndef ARRAY_H
#define ARRAY_H
#define NDEBUG 
#include<mm_malloc.h>
#define VECT_ALLIGN_SIZE 4096

using T = double;

template <typename T> class array2D
{
public:
    T * mat;
    std::size_t n1;
    std::size_t n2; 
    std::size_t size;
    
    array2D(const std::size_t m1,const std::size_t m2) 
    {
        mat=(T*) _mm_malloc(m2*m1*sizeof(T), VECT_ALLIGN_SIZE);
        n1 = m1;
        n2 = m2;
        size = m1*m2;
    }
    
//element
    T& operator()(const std::size_t& i, const std::size_t& j) { return this->mat[i*n2+j]; }
    const T& operator()(const std::size_t& i, const std::size_t& j) const { return this->mat[i*n2+j]; }
    T& operator[](const std::size_t& i) { return this->mat[i]; }
    const T& operator[](const std::size_t& i) const { return this->mat[i]; }    
};

#endif
