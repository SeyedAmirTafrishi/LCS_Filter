clc
clear
%% Our points

X_A = 1;%center of circle A
Y_A = 2;
X_d = 3.01; %The circle that is between A and B
Y_d = 4;%Original coordinate
R_d=1; %Radius of the circle between A and B
X_B=5; %Center of circle B
Y_B=9;

%% Find the Tangent point on the Circle of x_d,y_d
x0 = [0 0 0 0]; 
global xd_1 yd_2
% Solver may increase computation time, may solve it algebrically in future.
REF = [X_d,Y_d,X_A,Y_A,R_d];
f = @(x) FindTangenfordm(x,REF); % function of dummy variable y
%fsolve doesnt give multiple solutons
F = fsolve(f,x0);
Point1(1,1) = real(xd_1(1,1)); 
Point2(1,1) = real(F(1,1));
Point1(2,1) = real(F(1,2));
Point2(2,1) = real(yd_2(1,1));

%% Find the angles of circle respect to center of the Circle A
%[ angle ] = calculate_vector_angle( x1, y1, X_A, Y_A ) we cant use it
%since our initial zero angle line is not making true angle between Cricle
%x,d and y,d and C_b
AngleCtangL=asin(sqrt((Point2(1,1)-Y_d)^2+(Point1(1,1)-X_d)^2)/sqrt((Y_d-Y_A)^2+(X_d-X_A)^2))
AngleCtangU=asin(sqrt((Point2(2,1)-Y_d)^2+(Point1(2,1)-X_d)^2)/sqrt((Y_d-Y_A)^2+(X_d-X_A)^2))

% The Circle B can be anywhere and there is not tangent point! So we have
% to have at least 6 cases
%Our three sides of the traingle O_A,O_d,O_B
D_Ad=sqrt((X_A-X_d)^2+(Y_A-Y_d)^2);
D_AB=sqrt((X_A-X_B)^2+(Y_A-Y_B)^2);
D_Bd=sqrt((X_B-X_d)^2+(Y_B-Y_d)^2);
if D_Ad>=D_AB && D_AB>=D_Bd

elseif D_Ad>D_Bd && D_Bd>D_AB
    
elseif D_AB>=D_Ad && D_Ad>=D_Bd
    
elseif D_AB>D_Bd && D_Bd>D_Ad
    
elseif D_Bd>=D_Ad && D_Ad>=D_AB
    
elseif D_Bd>D_AB && D_AB>D_Ad
    
end
AngleCBtang=asin(sqrt((Point2(2,1)-Y_d)^2+(Point1(2,1)-X_d)^2)/sqrt((Y_d-Y_A)^2+(X_d-X_A)^2))







