%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA110
% Project Title: Implementation of Shuffled Complex Evolution (SCE-UA)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, Shuffled Complex Evolution in MATLAB (URL: https://yarpiz.com/80/ypea110-shuffled-complex-evolution), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function [pop, SortOrder] = SortPopulation(pop)

    % Get Costs
    Costs = [pop.Cost];
    
    % Sort the Costs Vector
    [~, SortOrder] = sort(Costs);
    
    % Apply the Sort Order to Population
    pop = pop(SortOrder);

end