
void solvePoissonJacobi(parameters param, node &myGrid, std::size_t &iter)
{   
    //dirichletX(myGrid.phin, myGrid.nB, myGrid.nE);
    //dirichletY(myGrid.phin, myGrid.nB, myGrid.nE);
    //dirichletX(myGrid.phinp1, myGrid.nB, myGrid.nE);
    //dirichletY(myGrid.phinp1, myGrid.nB, myGrid.nE);

    for(std::size_t i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {
        myGrid.phin(i,myGrid.nB[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nB[1]]*myGrid.coordY[myGrid.nB[1]];
        myGrid.phinp1(i,myGrid.nB[1]) = myGrid.phin(i,myGrid.nB[1]);
        myGrid.phin(i,myGrid.nE[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nE[1]]*myGrid.coordY[myGrid.nE[1]];
        myGrid.phinp1(i,myGrid.nE[1]) = myGrid.phin(i,myGrid.nE[1]);
    }

    for(std::size_t j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
    {
        myGrid.phin(myGrid.nB[0],j) = 1.0 + myGrid.coordX[myGrid.nB[0]]*myGrid.coordX[myGrid.nB[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
        myGrid.phinp1(myGrid.nB[0],j) = myGrid.phin(myGrid.nB[0],j);
        myGrid.phin(myGrid.nE[0],j) = 1.0 + myGrid.coordX[myGrid.nE[0]]*myGrid.coordX[myGrid.nE[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
        myGrid.phinp1(myGrid.nE[0],j) = myGrid.phin(myGrid.nE[0],j);
    }

    for(std::size_t i=myGrid.nB[0]+1; i<=myGrid.nE[0]-1; i++)
    {
        for(std::size_t j=myGrid.nB[1]+1; j<=myGrid.nE[1]-1; j++)
        {
            myGrid.phinp1(i,j) = 0.25*(myGrid.source(i,j)*param.dx*param.dx + myGrid.phin(i+1,j) + myGrid.phin(i-1,j) 
                                    + myGrid.phin(i,j+1) + myGrid.phin(i,j-1) );    
        }
    }
    iter++;

    //dirichletX(myGrid.phinp1, myGrid.nB, myGrid.nE);
    //dirichletY(myGrid.phinp1, myGrid.nB, myGrid.nE);
    //dirichletX(myGrid.phin, myGrid.nB, myGrid.nE);
    //dirichletY(myGrid.phin, myGrid.nB, myGrid.nE);

    for(std::size_t i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {
        myGrid.phin(i,myGrid.nB[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nB[1]]*myGrid.coordY[myGrid.nB[1]];
        myGrid.phinp1(i,myGrid.nB[1]) = myGrid.phin(i,myGrid.nB[1]);
        myGrid.phin(i,myGrid.nE[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nE[1]]*myGrid.coordY[myGrid.nE[1]];
        myGrid.phinp1(i,myGrid.nE[1]) = myGrid.phin(i,myGrid.nE[1]);
    }

    for(std::size_t j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
    {
        myGrid.phin(myGrid.nB[0],j) = 1.0 + myGrid.coordX[myGrid.nB[0]]*myGrid.coordX[myGrid.nB[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
        myGrid.phinp1(myGrid.nB[0],j) = myGrid.phin(myGrid.nB[0],j);
        myGrid.phin(myGrid.nE[0],j) = 1.0 + myGrid.coordX[myGrid.nE[0]]*myGrid.coordX[myGrid.nE[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
        myGrid.phinp1(myGrid.nE[0],j) = myGrid.phin(myGrid.nE[0],j);
    }

    for(std::size_t i=myGrid.nB[0]+1; i<=myGrid.nE[0]-1; i++)
    {
        for(std::size_t j=myGrid.nB[1]+1; j<=myGrid.nE[1]-1; j++)
        {   
            myGrid.phin(i,j) = 0.25*(myGrid.source(i,j)*param.dx*param.dx + myGrid.phinp1(i+1,j) + myGrid.phinp1(i-1,j) 
                                    + myGrid.phinp1(i,j+1) + myGrid.phinp1(i,j-1) );
        }
    }
    iter++;
}

void solvePoissonGaussSeidel(parameters param, node &myGrid, std::size_t &iter)
{
    //dirichletX(myGrid.phin, myGrid.nB, myGrid.nE);
    //dirichletY(myGrid.phin, myGrid.nB, myGrid.nE);

    for(std::size_t i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {
        myGrid.phin(i,myGrid.nB[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nB[1]]*myGrid.coordY[myGrid.nB[1]];
        myGrid.phin(i,myGrid.nE[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nE[1]]*myGrid.coordY[myGrid.nE[1]];
    }

    for(std::size_t j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
    {
        myGrid.phin(myGrid.nB[0],j) = 1.0 + myGrid.coordX[myGrid.nB[0]]*myGrid.coordX[myGrid.nB[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
        myGrid.phin(myGrid.nE[0],j) = 1.0 + myGrid.coordX[myGrid.nE[0]]*myGrid.coordX[myGrid.nE[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
    }
    
    for(std::size_t i=myGrid.nB[0]+1; i<=myGrid.nE[0]-1; i++)
    {
        for(std::size_t j=myGrid.nB[1]+1; j<=myGrid.nE[1]-1; j++)
        {
            myGrid.phin(i,j) = 0.25*(myGrid.source(i,j)*param.dx*param.dx + myGrid.phin(i+1,j) + myGrid.phin(i-1,j) 
                                    + myGrid.phin(i,j+1) + myGrid.phin(i,j-1) );    
        }
    }
    iter++;
}

void solvePoissonSOR(parameters param, node &myGrid, std::size_t &iter)
{
    //dirichletX(myGrid.phin, myGrid.nB, myGrid.nE);
    //dirichletY(myGrid.phin, myGrid.nB, myGrid.nE);

    for(std::size_t i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {
        myGrid.phin(i,myGrid.nB[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nB[1]]*myGrid.coordY[myGrid.nB[1]];
        myGrid.phin(i,myGrid.nE[1]) = 1.0 + myGrid.coordX[i]*myGrid.coordX[i] + 2.0*myGrid.coordY[myGrid.nE[1]]*myGrid.coordY[myGrid.nE[1]];
    }

    for(std::size_t j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
    {
        myGrid.phin(myGrid.nB[0],j) = 1.0 + myGrid.coordX[myGrid.nB[0]]*myGrid.coordX[myGrid.nB[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
        myGrid.phin(myGrid.nE[0],j) = 1.0 + myGrid.coordX[myGrid.nE[0]]*myGrid.coordX[myGrid.nE[0]] + 2.0*myGrid.coordY[j]*myGrid.coordY[j];
    }
    
    for(std::size_t i=myGrid.nB[0]+1; i<=myGrid.nE[0]-1; i++)
    {
        for(std::size_t j=myGrid.nB[1]+1; j<=myGrid.nE[1]-1; j++)
        {
            myGrid.phin(i,j) = param.omega*0.25*(myGrid.source(i,j)*param.dx*param.dx + myGrid.phin(i+1,j) + myGrid.phin(i-1,j) 
                                    + myGrid.phin(i,j+1) + myGrid.phin(i,j-1) ) + (1.0-param.omega)*myGrid.phin(i,j);    
        }
    }
    iter++;
}
