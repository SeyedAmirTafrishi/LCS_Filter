function [Edge] = Line(lambda,psi,Edge)%lambda(L,BS) , psi(L,R,Ty)

global ploti
set(0,'DefaultTextInterpreter','Latex');
c9 = Edge;
CB = 0;

%% Edge Remover by Psi
% new C9 after step 1.
% primary analysis of edges
% Note: Psi here considered Nx4 not Nx5 maybe in future quadrangular would add that require 1 more
% Edge Remover by Psi Note: Psi here considered Nx4 Psi(X,Y,Radius,Type)
if psi == 0
	% skip as nothing to remove;
else
    for i = 1:numel(psi(:,1)) % Counter for edges
        if psi(i,4) == 1      % Check whether the Psi is Circle or Square
            j = 1;
            while j <= numel(c9(:,1)) % Counter for edges
                if (((((psi(i,1) - c9(j,1))^2) + ((psi(i,2)-c9(j,2))^2))^(0.5)) <= psi(i,3))
                    % The condition to check the boundery property
                    c9(j,:) = []; % Remove the Edge that has to be ignored
                    j = j - 1;
                end
                j = j + 1;
            end
        % elseif psi(:,4) == 2  %WILL BE FILLED for SQUARE
            % WILL BE FILLED
        end

        % can be delete
        th = 0:pi/50:2*pi; % for loop for creating circle
        xunit = (psi(i,3) + CB) * cos(th) + psi(i,2); % equation of circle :D
        yunit = (psi(i,3) + CB) * sin(th) + psi(i,1);
        hold on
        subplot(2,2,1)
        ploti = plot(xunit, yunit, 'y'); % Plot the boys :v
        xlim([1 640])
        ylim([1 480])
    end
end


subplot(2,2,1)
hold on
ploti = plot(c9(:,2), c9(:,1), 'y.'); %edges

%% Edge grouper by BS of lambda
% step 2. to collect information in a better way
% try to group the edges roughly
% location and boundry size
% will be used later in rebel and circle analysis
% Lambda group the edges in column by column

DV = 0; % Data variable for categorizing the edges in tree form collections
f = 1;  % the counter of new column edge groups
k = 0;  % shifter triger for new column edge groups;
        % After checking the Nth Lambda for all newly detected edges. It set k=1 to add new couple column
%-------in Lambda
for z = 1:numel(lambda(:,1)) % Lambda Counter
    %------------- % Matrix Shifter
    if k == 1
        f = f + 2;   % shift two for x and y
        k = 0;
    end
    hold on
    %--------------------- The Lambda edge comparator
    j = 1;
    g = 1;
    if lambda == 0
        % pass
    else
        while j <= numel(c9(:,1))
            % The condition to check the boundery property
            if (((((lambda(z,1) - c9(j,1))^2) + ((lambda(z,2) - c9(j,2))^2))^(0.5)) <= lambda(z,3))
                subplot(2,2,1)
                plot([lambda(z,2), c9(j,2)], [lambda(z,1), c9(j,1)], 'g'); % it plots the related lines when condition satisfied
                DV(g,f) = c9(j,1); % include the edges to the relevant group
                DV(g,f + 1) = c9(j,2);
                c9(j,:) = [];
                k = 1; % Activation for new column couple
                g = g + 1;
            else
                j = j + 1;
            end
        end
    end
end % end of for
%-----------------Boxing part :D

%%
% Out of Lambda members and psi
for j = 1:numel(c9(:,1)) % Counter for edges
    f = f + 2;
    %plot(c9(j,2),c9(j,1),'r'); % it plots the related lines when condition satisfied
    DV(j,f) = c9(j,1);   % include the edges to the relevant group
    DV(j,f+1) = c9(j,2);
end

hold on
Edge = DV;

end % end of function
