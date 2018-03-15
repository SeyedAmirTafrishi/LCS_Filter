function [Edge] = Line(lambda,psi,Edge)%lambda(L,BS) , psi(L,R,Ty)

global ploti
set(0,'DefaultTextInterpreter','Latex');
c9 = Edge;
CB = 0;

%% Edge Remover by Psi
% Note: Psi here considered Nx4 not Nx5 maybe in future quadrangular would add that require 1 more
if psi == 0
else
    for i = 1:numel(psi(:,1)) %Counter for edges
        if psi(i,4) == 1
          j = 1;
          while j <= (numel(c9)/2) %Counter for edges
            if (((((psi(i,1)-c9(j,1))^2) + ((psi(i,2)-c9(j,2))^2))^(0.5)) <= psi(i,3))
              %The condition to check the boundery property
              c9(j,:) = [];
              j = j - 1;
            end
            j = j + 1;
          end
          %elseif psi(:,4)==2%WILL BE FILLED for SQUARE
          %WILL BE FILLED
        end
        
        th = 0:pi/50:2*pi;%for loop for creating circle
        xunit = (psi(i,3) + CB) * cos(th) + psi(i,2);%equation of circle :D
        yunit = (psi(i,3) + CB) * sin(th) + psi(i,1);
        hold on
        subplot(1,2,1)
        ploti = plot(xunit, yunit,'y');%Plot the boys :v
        xlim([1 640])
        ylim([1 480])
    end
end

hold on
subplot(1,2,1)
ploti = plot(c9(:,2),c9(:,1),'y.'); %edges

%% Edge grouper by BS of lambda
DV = 0;%Data Variable for categorizing the edges in tree form collections
f = 1; %the counter of new column edge groups
k = 0; %shifter triger for new column edge groups
i = 1;% simple counter
%-------in Lambda
for z = 1:(numel(lambda(:,1))) %Lambda Counter
    %------------- %Matrix Shifter
    if k == 1
      f = f + 2;%shift two for x and y
      k = 0;
    end
    hold on
    %--------------------- The Lambda edge comparator
    j = 1;
    g = 1;
    if lambda == 0
    else
      while j <= (numel(c9)/2)
        if (((((lambda(z,1)-c9(j,1))^2) + ((lambda(z,2)-c9(j,2))^2))^(0.5)) <= lambda(z,3)) %The condition to check the boundery property
          subplot(1,2,1)
          plot([lambda(z,2), c9(j,2)],[lambda(z,1), c9(j,1)], 'g'); %it plots the related lines when condition satisfied
          DV(g,f) = c9(j,1);%include the edges to the relevant group
          DV(g,f+1) = c9(j,2);
          c9(j,:) = [];
          k = 1;%Activation for new column couple
          g = g+1;
        else
          j = j+1;
        end
      end
    end
%--------------------------------------
end % foraaa\
%-----------------Boxing part :D
%%
%Out of Lambda Members
for j = 1:(numel(c9)/2) %Counter for edges
    f = f + 2;
    %plot(c9(j,2),c9(j,1),'r'); %it plots the related lines when condition satisfied
    DV(j,f) = c9(j,1);%include the edges to the relevant group
    DV(j,f+1) = c9(j,2);
end
hold on
Edge = DV;

end % end of function