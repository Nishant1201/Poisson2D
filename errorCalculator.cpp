
void calculatenorm(parameters &param, node myGrid)
{
    param.normL2 =0.0;
    param.normL1 =0.0;
    param.normLinfinity = 0.0;
    
    for(std::size_t i=myGrid.nB[0]; i<=myGrid.nE[0]; i++)
    {
        for(std::size_t j=myGrid.nB[1]; j<=myGrid.nE[1]; j++)
        {
            T error = fabs(myGrid.phin(i,j)-myGrid.phi_exact(i,j));
            param.normL1 += error;
            param.normL2 += error*error;
            if(error>param.normLinfinity)
                param.normLinfinity = error;
        }
    }
    
    std::size_t normalization_factor = (myGrid.nX[0]-2*myGrid.nB[0])*(myGrid.nX[1]-2*myGrid.nB[1]);
    param.normL1 = param.normL1/normalization_factor;
    param.normL2 = sqrt(param.normL2/normalization_factor);
    param.normLinfinity = param.normLinfinity; 
}


