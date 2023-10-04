% 
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Code: YPEA126
% Project Title: Non-dominated Sorting Genetic Algorithm III (NSGA-III)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD (member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints,"
% in IEEE Transactions on Evolutionary Computation,
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Papaer URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

function y=Mutate(x,mu,sigma,VarMin,VarMax)

    nVar=numel(x);
    
    nMu=ceil(mu*nVar);

    j=randsample(nVar,nMu);
    
    if numel(sigma)>1
        sigma = sigma(j);
    end
  
    y=x;
    
    y(j)=x(j)+sigma.*randn(size(j)); % Mutando e criando o novo vetor de parametros 
    
    while (y(j)<VarMin(j)) | (VarMax(j)<y(j))
        y(j)=x(j)+sigma.*randn(size(j));
        %fprintf('%d ,erro, y= %4.2f \n',j,y(j));
    end
end