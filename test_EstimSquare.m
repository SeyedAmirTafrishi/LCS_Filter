%%
clc
clear


%global x_1 y_2
x_1 = 0;
y_2 = 0;


X_e = 1;        % Elipse center coordinate
Y_e = 2;
X_o = -3.01;
Y_o = 4;        % Original coordinate
a_1 = 4;        % elipse minor past
b_1 = .5;       % elipse major past
Delta_r = 1.6;  % Estimated shift on R_e

tic
[Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r)
toc