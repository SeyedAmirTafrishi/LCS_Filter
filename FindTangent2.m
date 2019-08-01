function F = FindTangent(x)
X_e= 1;%Elipse center coordinate
Y_e= 3;
Y_o=2;%Original coordinate
X_o=3;
a_1=1; %elipse minor past
b_1=.5; %elipse major past
%global a_1 b_1 Y_e X_e X_o Y_o
%x_2=Y_e-((b_1/a_1)*sqrt((a_1)^2-(x(1)-X_e)^2))
% x_1=X_e-((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2))
% F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
%F(2) = ((x(2)-Y_e)^2/(b_1^2))+((x(1)-X_e)^2/(a_1^2))-1;
% %Solution 1 
% x_1=X_e+((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2))
% F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
%Solution 2
x_2=Y_e-((b_1/a_1)*sqrt((a_1)^2-(x(1)-X_e)^2));
F(1) = (a_1^2*(Y_o-x_2)*(x_2-Y_e))+(b_1^2*(X_o-x(1))*(x(1)-X_e));