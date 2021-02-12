#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <fstream>
#include <chrono>

#include "array.h"
#include "grid.h"
#include "initialize.cpp"
#include "boundaryconditions.cpp"
#include "solvers.cpp"
#include "errorCalculator.cpp"
#include "postprocess.cpp"
 
using namespace std;

int main()
{   
    // using T = double;

    ofstream out;
    
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> dif(0.0);
    
    std::size_t nX, nY;                                             // number of grid points along X,Y,Z
    std::size_t PAD;                                                // number of ghost nodes along each directions (X,Y,Z)
    nX = 50;                                                                        
    nY = 50;                    
    PAD = 0;
    std::size_t sizeX(nX+PAD*2), sizeY(nY+PAD*2);                   // Array size including ghost nodes and computational nodes (Array.h)

    node myGrid(sizeX,sizeY,PAD);                                   // node class (see file grid.h)

    // Set up the Grid
    T len = 1.0;                                                    // Domain length (along X,Y,Z)
    T dx = len/(nX-1);                                              // delta x - distance between grid points along X
    T dy = len/(nY-1);                                              // delta y - distance between grid points along Y
    
    // Set the parameters
    parameters param;                                           // class for various parameters used in problem (see file initialize.cpp) 
    setparameters(param, dx, dy);                              // set the parameters (function setparameters included in initialize.cpp)
    
    std::size_t iter = 0;
    std::size_t freq = 1;

    InitializeField(param, myGrid);                             // Domain initialization function, (initialize.cpp)

    start = std::chrono::system_clock::now();
    while(param.normL2 >= 1.0e-4 && iter <=10000)
    {   
        //solvePoissonJacobi(param, myGrid, iter);
        //solvePoissonGaussSeidel(param, myGrid, iter);
        solvePoissonSOR(param, myGrid, iter);
        calculatenorm(param, myGrid); 

        //std::cout.precision(12);
        //std::cout<<nX<<"    "<<iter<<"      "<<param.normL1<<"  "<<param.normL2<<"  "<<param.normLinfinity<<std::endl;
    }
    end = std::chrono::system_clock::now();
    dif = end - start;
    
    //writeSolutionfile(myGrid);
    std::cout.precision(12);
    std::cout<<nX<<"    "<<iter<<"      "<<dif.count()<<"  "<<param.normL1<<"  "<<param.normL2<<"  "<<param.normLinfinity<<std::endl;

   return 0;
}

