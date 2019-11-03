%%
clc
clear


%global x_1 y_2
x_1 = 0;
y_2 = 0;
%Problem in MINUS?

X_e = 300;        % Elipse center coordinate
Y_e = 300;
X_o = 303+.1;
Y_o= 303+.1;        % Original coordinate
a_1 = 10.01;        % elipse minor past
b_1 = 5.01;       % elipse major past
Delta_r = 1.6;  % Estimated shift on R_e

tic
[Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r)
toc