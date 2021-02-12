
#ifndef GRID_H
#define GRID_H

class node
{
public:
    node(const std::size_t sizeX1, const std::size_t sizeX2, const std::size_t PAD):
    phin(sizeX1,  sizeX2), phinp1(sizeX1, sizeX2),
    phi_exact(sizeX1, sizeX2), source(sizeX1, sizeX2)
    {
        nX[0] = sizeX1;                                                         // Total number of grid points along X
        nX[1] = sizeX2;                                                         // Total number of grid points along Y
        nB[0] = PAD;                                                            // Beginning of actual node along X
        nB[1] = PAD;                                                            // Beginning of actual node along Y
        nE[0] = sizeX1-PAD-1;                                                   // End of actual node along X
        nE[1] = sizeX2-PAD-1;                                                   // End of actual node along Y
    }
    std::size_t nX[2];
    std::size_t nB[2];
    std::size_t nE[2];
    array2D<T> phin;
    array2D<T> phinp1;
    array2D<T> phi_exact;
    array2D<T> source;
    T *coordX;                                                                   // Array for x-coordinates 
    T *coordY;                                                                   // Array for y-coordinates
};

#endif