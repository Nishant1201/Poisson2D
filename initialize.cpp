
class parameters
{
public:
    T twoPI;
    T factor;
    T diffusivity;
    T dx,dy;
    T omega;
    T normL1;
    T normL2;
    T normLinfinity;
};

void setparameters(parameters &param, T dx, T dy)
{   
    param.twoPI = 2.0*M_PI;
    param.factor = 8.0*M_PI*M_PI;
    param.diffusivity = 1.0;
    param.dx = dx;
    param.dy = dy;
    param.omega = 1.9;
    param.normL1 = 1.0;
    param.normL2 = 1.0;
    param.normLinfinity = 1.0;
}

void InitializeField(parameters param, node &myGrid)
{
    myGrid.coordX = (T*) _mm_malloc(myGrid.nX[0]*sizeof(T), VECT_ALLIGN_SIZE);
    myGrid.coordY = (T*) _mm_malloc(myGrid.nX[1]*sizeof(T), VECT_ALLIGN_SIZE);
    
    for(std::size_t i=0; i<myGrid.nX[0]; i++)
    {
        myGrid.coordX[i] = -(myGrid.nB[0])*param.dx + i*param.dx;                       // [xmin-dx,xmax+dx], xmin = 0.0, xmax = 1.0;
    //    std::cout<<myGrid.coordX[i]<<"  ";
    }   // std::cout<<std::endl;
    for(std::size_t j=0; j<myGrid.nX[1]; j++){
        myGrid.coordY[j] = -(myGrid.nB[1])*param.dy + j*param.dy;                       // [ymin-dy,ymax+dy], ymin = 0.0, ymax = 1.0;
    //    std::cout<<myGrid.coordY[j]<<"  ";
    }   // std::cout<<std::endl;
    
    /*
    for(int i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {       
        myReal xval = sin(param.twoPI*myGrid.coordX[i]);
        for(int j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
        {   
            myGrid.phi_exact(i,j) = xval*sin(param.twoPI*myGrid.coordY[j]); 
            myGrid.source(i,j) = param.factor*myGrid.phi_exact(i,j); 
            myGrid.phin(i,j) = 0.0;
            myGrid.phinp1(i,j) = 0.0;
        }  
    }  
    //myGrid.source(int(myGrid.nE[1]/4),int(myGrid.nE[0]/4)) = 100.0;
    //myGrid.source(int(myGrid.nE[1]/4),int(myGrid.nE[0]/4)) = -100.0;
    */

    for(std::size_t i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {
        for(std::size_t j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
        {
            myGrid.phi_exact(i,j) = 1 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
            myGrid.source(i,j) = -6.0;
            myGrid.phin(i,j) = 0.0;
            myGrid.phinp1(i,j) = 0.0;
        }
    }  
}
