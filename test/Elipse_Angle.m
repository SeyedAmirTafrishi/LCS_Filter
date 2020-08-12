%%
clc
clear

X_e = 1;%Elipse center coordinate
Y_e = 2;
X_o = -3.01;
Y_o = 4;%Original coordinate
a_1 = 4; %elipse minor past
b_1 = .5; %elipse major past
%syms X Y
tic
x0 = [0 0 0 0]; 
global x_1 y_2
% Solver may increase computation time, may solve it algebrically in future.
REF = [X_e,Y_e,X_o,Y_o,a_1,b_1];
f = @(x) FindTangentx1(x,REF); % function of dummy variable y
%fsolve doesnt give multiple solutons
F = fsolve(f,x0);
Point1(1,1) = real(x_1(1,1));
Point2(1,1) = real(F(1,1));
Point1(2,1) = real(F(1,2));
Point2(2,1) = real(y_2(1,1));
%  
%%
%% Plot the results!
th = 0:pi/50:2*pi;%for loop for creating circle
CB = 1;
xunit = (a_1 ) * cos(th) + X_e;%equation of circle :D
yunit = (b_1 ) * sin(th) + Y_e;
ploti = plot(xunit, yunit,'g');% Ellipse
hold on
plot(X_o,Y_o,'- *b','MarkerSize', 18,'LineWidth' , 2.5)
plot(Point1(1,1),Point2(1,1),'- xr','MarkerSize', 18,'LineWidth' , 2.5)
plot(Point1(2,1),Point2(2,1),'- xr','MarkerSize', 18,'LineWidth' , 2.5)
plot(X_e,Y_e,'- om','MarkerSize', 18,'LineWidth' , 2.5)

%   R_e R_o
R_e = sqrt((Y_e-Y_o)^2 + (X_e-X_o)^2); %Orignal to Elipse
R_o1 = sqrt((Point2(1,1)-Y_o)^2 + (Point1(1,1)-X_o)^2); %Orignal to contact 1
R_o2 = sqrt((Point2(2,1)-Y_o)^2 + (Point1(2,1)-X_o)^2); %Orignal to contact 2
%for two different side!!! trignomatic
if R_e > R_o1 %Hypto
  alpha1 = acos(R_o1/R_e);
else
  alpha1 = acos(R_e/R_o1);   
end

if R_e > R_o2
  alpha2 = acos(R_o2/R_e);    
else
  alpha2 = acos(R_e/R_o2);     
end

%% Estimated Delta
Delta_r = 1.6; %Estimated shift on R_e
    if R_e > R_o1 %Hypto
      DRo1 = (R_e+Delta_r)*cos(alpha1) - R_o1;
    else
      DRo1 = (1/cos(alpha1))*((R_e+Delta_r) - R_o1*cos(alpha1));
    end

    if R_e > R_o2 %Hypto
      DRo2 = (R_e+Delta_r)*cos(alpha2) - R_o2;
    else
      DRo2 = (1/cos(alpha2))*((R_e+Delta_r) - R_o2*cos(alpha2));
    end
%% Finding New Tangent points
   m_ro1 = ((Point1(1,1)-X_o)/(Point2(1,1)-Y_o)); %the slope of tangent line Ro1
   m_ro2 = ((Point1(2,1)-X_o)/(Point2(2,1)-Y_o));%the slope of tangent line Ro2
   m_re = ((X_e-X_o)/(Y_e-Y_o)); %the slope of Re (elipse position)
   %% Solving for finding (X_n1, Y_n1) and (X_n2,Y_n2)
%syms X_n1 Y_n1 X_n2 Y_n2 % Slow solution
%yn0=0;
%xn0=0;
  % Solver may increase computation time, may solve it algebrically in future.
%REF1 = [Point1(1,1),Point2(1,1),DRo1,m_ro1];
%fn = @(x) NewTP(x,REF1); % function of dummy variable y
%fsolve doesnt give multiple solutons
%y_n1 = fsolve(fn,[xn0 yn0]) 
%FN1=(X_n1-Point1(1,1))^2+(Y_n1-Point2(1,1))^2==DRo1^2;
%FN2=(m_ro1*(Y_n1-Point2(1,1)))-(X_n1-Point1(1,1))==0;
%sold = solve([FN1, FN2], [X_n1, Y_n1]); % Solver may increase computation time, may solve it algebrically in future.
%Pointn1a = real(y_n1(1,2))
%Pointn1b = real(y_n1(1,1))
% FM1=(X_n2-Point1(2,1))^2+(Y_n2-Point2(2,1))^2==DRo2^2;
% FM2=(m_ro2*(Y_n2-Point2(2,1)))-(X_n2-Point1(2,1))==0;
% sold = solve([FM1, FM2], [X_n2, Y_n2]); % Solver may increase computation time, may solve it algebrically in future.
% Pointn2a = double(sold.X_n2)
% Pointn2b = double(sold.Y_n2)
% if sqrt((Pointn1b(1,1)-Y_o)^2+(Pointn1a(1,1)-X_o)^2) > R_o1+DRo1 %Length must be longer than existing one always to find right points
% X_nf1=Pointn1a(1,1);
% Y_nf1=Pointn1b(1,1);
% else
% X_nf1=Pointn1a(2,1);
% Y_nf1=Pointn1b(2,1);  
% end
% if sqrt((Pointn2b(1,1)-Y_o)^2+(Pointn2a(1,1)-X_o)^2) > R_o2+DRo2 %Length must be longer than existing one always to find right points
% X_nf2=Pointn2a(1,1);
% Y_nf2=Pointn2b(1,1);
% else
% X_nf2=Pointn2a(2,1);
% Y_nf2=Pointn2b(2,1);  
% end
if Y_o < Y_e % working Solution
  Y_nf1 = Point2(1,1)+sqrt(DRo1^2/(m_ro1^2+1));
  X_nf1 = Point1(1,1)+m_ro1*(Y_nf1-Point2(1,1));    
else
  Y_nf1 = Point2(1,1)-sqrt(DRo1^2/(m_ro1^2+1));
  X_nf1 = Point1(1,1)+m_ro1*(Y_nf1-Point2(1,1));    
end

if Y_o < Y_e % working Solution
  Y_nf2 = Point2(2,1)+sqrt(DRo2^2/(m_ro2^2+1));
  X_nf2 = Point1(2,1)+m_ro2*(Y_nf2-Point2(2,1));    
else
  Y_nf2 = Point2(2,1)-sqrt(DRo2^2/(m_ro2^2+1));
  X_nf2 = Point1(2,1)+m_ro2*(Y_nf2-Point2(2,1));    
end
%Algebric results
plot(X_nf1,Y_nf1,'- xy','MarkerSize', 18,'LineWidth' , 2.5) 
plot(X_nf2,Y_nf2,'- xy','MarkerSize', 18,'LineWidth' , 2.5)     
    
%% Solving Equations of Elipse for a_n and b_n
% % Estimated location of center of elipse
% syms X_ne Y_ne 
% EN1=(X_ne-X_e)^2+(Y_ne-Y_e)^2==(Delta_r)^2;
% EN2=(m_re*(Y_ne-Y_e))-(X_ne-X_e)==0;
% sold = solve([EN1, EN2], [X_ne, Y_ne]); % Solver may increase computation time, may solve it algebrically in future.
% XE = double(sold.X_ne);
% YE = double(sold.Y_ne);
% if sqrt((YE(1,1)-Y_o)^2+(XE(1,1)-X_o)^2) > R_e+Delta_r %Length must be longer than existing one always to find right points
% X_ef1=XE(1,1);
% Y_ef1=YE(1,1);
% else
% X_ef1=XE(2,1);
% Y_ef1=YE(2,1);  
% end
%% Solved Algebriclly
if Y_o < Y_e % working Solution
  Y_ef1 = Y_e+sqrt(Delta_r^2/(m_re^2+1));
  X_ef1 = X_e+m_re*(Y_ef1-Y_e);    
else
  Y_ef1 = Y_e-sqrt(Delta_r^2/(m_re^2+1));
  X_ef1 = X_e+m_re*(Y_ef1-Y_e);    
end

plot(X_ef1,Y_ef1,'- ok','MarkerSize', 18,'LineWidth' , 2.5) 

% syms a_n b_n
% Fa1=(a_n^2*(Y_o-Y_nf1)*(Y_nf1-Y_ef1))+(b_n^2*(X_o-X_nf1)*(X_nf1-X_ef1))==0;
% Fb2=((Y_nf1-Y_ef1)^2/(b_n^2))+((X_nf1-X_ef1)^2/(a_n^2))==1;
% 
% sol = solve([Fa1, Fb2], [a_n, b_n]); % Solver may increase computation time, may solve it algebrically in future.
% A_n = double(sol.a_n)
% B_n = double(sol.b_n)
if abs((X_o-X_nf1)) > .01
  B_n = sqrt(abs((((X_o-X_nf1)*(Y_nf1-Y_ef1)^2)-((X_nf1-X_ef1)*(Y_o-Y_nf1)*(Y_nf1-Y_ef1))) / ((X_o-X_nf1))))
  A_n = sqrt(abs((B_n^2*(X_o-X_nf1)*(X_nf1-X_ef1)) / ((Y_o-Y_nf1)*(Y_nf1-Y_ef1))))    
else %when it is zero Singular point
% deviate with 0.01
    
end

th = 0:pi/50:2*pi;%for loop for creating circle
CB = 1;
xunit = (abs(A_n) ) * cos(th) + X_ef1;%equation of circle :D
yunit = (abs(B_n) ) * sin(th) + Y_ef1;
ploti = plot(xunit, yunit,'b');% Ellipse

%% Approximate the point to largest, Output is a, b (minor and major radii) as well as the rounded X and Y (Think about property of variables Natural or DEcimal) 
 toc