% %%
% clc
% clear


%global x_1 y_2
x_1 = 0;
y_2 = 0;
%Problem in MINUS?

X_e = 389;        % Elipse center coordinate
Y_e = 290.0000;
X_o = 320.0000;
Y_o= 300.0000;        % Original coordinate
a_1 =  60.3278 ;        % elipse minor past
b_1 =   60.3278 ;       % elipse major past
Delta_r = 3;  % Estimated shift on R_e

% 487.0000  164.0000  320.0000  240.0000  160.3278  160.3278
tic
[Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r)
toc