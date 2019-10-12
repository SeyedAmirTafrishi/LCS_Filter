clc
clear
%% Our points

X_A = 6;%center of circle A
Y_A = 3.6
X_d = 3.2; %The circle that is between A and B
Y_d = 6;%Original coordinate
R_d=1; %Radius of the circle between A and B
X_B=5; %Center of circle B
Y_B=5;

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
AngleCtangL=atan2(sqrt((Point2(1,1)-Y_d)^2+(Point1(1,1)-X_d)^2),sqrt((Y_d-Y_A)^2+(X_d-X_A)^2))
AngleCtangU=asin(sqrt((Point2(2,1)-Y_d)^2+(Point1(2,1)-X_d)^2)/sqrt((Y_d-Y_A)^2+(X_d-X_A)^2))

% The Circle B can be anywhere and there is not tangent point! So we have
% to have at least 6 cases
%Our three sides of the traingle O_A,O_d,O_B
D_Ad=sqrt((X_A-X_d)^2+(Y_A-Y_d)^2);
D_AB=sqrt((X_A-X_B)^2+(Y_A-Y_B)^2);
D_Bd=sqrt((X_B-X_d)^2+(Y_B-Y_d)^2);
% if D_Ad>=D_AB && D_AB>=D_Bd
%     1
% AngleCBD=asin(D_Bd/D_Ad);
% elseif D_Ad>D_Bd && D_Bd>D_AB
%     2
% AngleCBD=asin(D_Bd/D_Ad);
% elseif D_AB>=D_Ad && D_Ad>=D_Bd
%     3
% Ahh=acos((D_AB^2+D_Ad^2-D_Bd^2)/(2*D_AB*D_Ad))   
% AngleCBD=asin(D_Bd/D_AB);    
% elseif D_AB>D_Bd && D_Bd>D_Ad
%     4
% AngleCBD=atan2(D_Bd,D_Ad);    
% elseif D_Bd>=D_Ad && D_Ad>=D_AB
%     5
% AngleCBD=acos(D_Ad/D_Bd);     
% elseif D_Bd>D_AB && D_AB>D_Ad
%     7
% AngleCBD=acos(D_Ad/D_Bd);   
% end
% AngleCBD
Ahh=acos((D_AB^2+D_Ad^2-D_Bd^2)/(2*D_AB*D_Ad))
th = 0:pi/50:2*pi;%for loop for creating circle
CB = 1;
xunit = (R_d ) * cos(th) + X_d;%equation of circle :D
yunit = (R_d ) * sin(th) + Y_d;
ploti = plot(xunit, yunit,'g');% Ellipse
hold on
plot(X_d,Y_d,'- *b','MarkerSize', 18,'LineWidth' , 2.5)
plot(Point1(1,1),Point2(1,1),'- xr','MarkerSize', 18,'LineWidth' , 2.5)
plot(Point1(2,1),Point2(2,1),'- xr','MarkerSize', 18,'LineWidth' , 2.5)
plot(X_A,Y_A,'- om','MarkerSize', 18,'LineWidth' , 2.5)
plot(X_B,Y_B,'- ob','MarkerSize', 18,'LineWidth' , 2.5)
plot([X_A Point1(1,1)],[Y_A Point2(1,1)],'y','MarkerSize', 18,'LineWidth' , 1.5)
plot([X_A Point1(2,1)],[Y_A Point2(2,1)],'y','MarkerSize', 18,'LineWidth' , 1.5)

