function F = FindTangenfordm(x, REF) %NOTEEEE: the function is obtained from FindTangentx1.m
X_d = REF(1,1); %the circle between A and B
Y_d = REF(1,2);
X_A = REF(1,3);
Y_A = REF(1,4); %Original coordinate
R = REF(1,5); % radius of circle between A and B

global xd_1 yd_2

% Solution 2
if (X_A>X_d && Y_A>Y_d) || (X_A<X_d && Y_A>Y_d) %+ + 
    xd_1 = X_d + ((R/R)*sqrt((R)^2-(x(1)-Y_d)^2));
    % save x_1
    F(1) = (R^2*(Y_A-x(1))*(x(1)-Y_d))+(R^2*(X_A-xd_1)*(xd_1-X_d));
    yd_2 = Y_d + ((R/R)*sqrt((R)^2-(x(2)-X_d)^2));  %change this guy
    % save y_2
    F(2) = (R^2*(Y_A-yd_2)*(yd_2-Y_d)) + (R^2*(X_A-x(2))*(x(2)-X_d));  
else %- +
    xd_1 = X_d + ((R/R)*sqrt((R)^2-(x(1)-Y_d)^2)); 
    % save x_1
    F(1) = (R^2*(Y_A-x(1))*(x(1)-Y_d)) + (R^2*(X_A-xd_1)*(xd_1-X_d));
    yd_2 = Y_d - ((R/R)*sqrt((R)^2-(x(2)-X_d)^2));  %change this guy
    % save y_2
    F(2) = (R^2*(Y_A-yd_2)*(yd_2-Y_d)) + (R^2*(X_A-x(2))*(x(2)-X_d));    
    
    %elseif if X_o=X_e    
end
 

end
