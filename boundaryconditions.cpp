
// Dirichlet boundary conditions routine
template<typename T>
void dirichletX(array2D<T> A, std::size_t *nB, std::size_t *nE)
{
    for(std::size_t j=nB[1]; j<=nE[1]; j++)
    {
        A(nB[0],j) = 0.0;
        A(nE[0],j) = 0.0;
    }
}

template<typename T>
void dirichletY(array2D<T> A, std::size_t *nB, std::size_t *nE)
{
    for(std::size_t i=nB[0]; i<=nE[0]; i++)
    {
        A(i,nB[1]) = 0.0;
        A(i,nE[1]) = 0.0;
    }
}

//  Periodic boundary Conditions routine //
template<typename T>
void periodicX(array2D<T> A, std::size_t *nB, std::size_t *nE)
{
    std::size_t offset = nE[0] - nB[0] + 1;
    for(std::size_t i=0; i<nB[0]; i++)
        for(std::size_t j=0; j<=nB[1]+nE[1]; j++)
        {
            A(i,j) = A(i+offset,j);
            A(i+offset+nB[0],j) = A(i+nB[0],j);
        }            
}

template<typename T>
void periodicY(array2D<T> A, std::size_t *nB, std::size_t *nE)
{
    std::size_t offset = nE[1] - nB[1] + 1;
    for(std::size_t i=0; i<=nB[0]+nE[0]; i++)
        for(std::size_t j=0; j<nB[1]; j++)
            {
                A(i,j) = A(i,j+offset);
                A(i,j+offset+nB[1]) = A(i,j+nB[1]);
            }
}

//  Neumann boundary Conditions routine //
template<typename T>
void neumannX(array2D<T> A, std::size_t *nB, std::size_t *nE)
{
    for(std::size_t j=nB[1]; j<=nE[1]; j++)
    {
        A(nB[0],j) = A(nB[0]+1,j);
        A(nE[0],j) = A(nE[0]-1,j);
    }
}

template<typename T>
void neumannY(array2D<T> A, std::size_t *nB, std::size_t *nE)
{
    for(std::size_t i=nB[0]; i<=nE[0]; i++)
    {
        A(i,nB[1]) = A(i,nB[1]+1);
        A(i,nE[1]) = A(i,nE[1]-1);
    }
}