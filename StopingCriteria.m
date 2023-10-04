%
%
% function [nF ImpBest ImpAv MovPar MovObj MaxDist StdDev Diff FO_Best] = StopingCriteria(F)
 
%function [nF ImpBest ImpAv MovPar MovObj MaxDist StdDev Diff FO_Best] =
%StopingCriteria(F)
%     [nF ~]= size(F);  % numero de membros de F
%     FO=zeros(nF,2);
%     for k=1:nF        
%         FO(k,:) = [F(k,:).Cost];
%     end
%     FO = sortrows(FO,1);
%     FO = FO(isfinite(FO(:, 1)), :);
%     [nF ~]= size(FO);
%     if nF > 1
%         FO_med = [mean(FO(:,1)) mean(FO(:,2))];
% 
%         d = (FO(:,1).^2 + FO(:,2).^2).^0.5;
% 
%         % Critérios baseados em melhoria
%         [ImpBest Imin] = min(d);
%         FO_Best = FO(Imin,:);
% 
%         ImpAv = (FO_med(1)^2 + FO_med(2)^2)^0.5;
% 
%         % Critérios baseados no movimento
%         if(nF == 1)
%             MovPar = FO;
%         else
%             MovPar = max(FO);
%         end
% 
%         MovObj = FO_med;
% 
%         % Critérios baseado em distribuição
%         if(nF == 1)
%             MaxDist = 0;
%             StdDev = 0;
%             Diff = 0;
%         else
%             dist = sort(((FO(:,1)-FO(Imin,1)).^2 + (FO(:,2)-FO(Imin,2)).^2).^0.5);
%             dist = dist(dist>0);
%             dist_med = mean(dist);
%             [MaxDist,~] = max(dist);
% 
%             E = (dist - dist_med).^2;
%             StdDev = (sum(E)^0.5)/numel(E);
% 
%             Diff = (sum((FO(1,:)-FO(nF,:)).^2))^0.5;
%         end
%     elseif nF == 1
%         d = (FO(:,1).^2 + FO(:,2).^2).^0.5;        
%         ImpBest = d;        
%         ImpAv = d;
%         MovPar = FO;
%         MovObj = FO;
%         MaxDist = d;        
%         FO_Best = FO;        
%         StdDev = 0;
%         Diff = 0;
%     elseif nF == 0
%         FO = [nan,nan];
%         d = nan;
%         ImpBest = d;        
%         ImpAv = d;
%         MovPar = FO;
%         MovObj = FO;
%         MaxDist = d;        
%         FO_Best = FO;        
%         StdDev = 0;
%         Diff = 0;
%     end
% end


function [nF GD] = StopingCriteria(F)
    [nF ~]= size(F);  % numero de membros de F
    FO=zeros(nF,2);
    for k=1:nF        
        FO(k,:) = [F(k,:).Cost];
    end
    FO = sortrows(FO,1);
    FO = FO(isfinite(FO(:, 1)), :);
    [nF ~]= size(FO);
    d = (FO(:,1).^2 + FO(:,2).^2);
    % Critérios baseados em melhoria
    GD = (sum(d)^0.5)/nF;
end