% %%
% clc
% clear


%global x_1 y_2
x_1 = 0;
y_2 = 0;
%Problem in MINUS?

X_e = 400.2;        % Elipse center coordinate
Y_e = 364.1000;
X_o = 440.1000;
Y_o= 320.1000;        % Original coordinate
a_1 =  50.3278 ;        % elipse minor past
b_1 =   35.3278 ;       % elipse major past
Delta_r = 20;  % Estimated shift on R_e
X_o-X_e
Y_o-Y_e
% 487.0000  164.0000  320.0000  240.0000  160.3278  160.3278
tic
[Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r)
toc     


% % PRoblem
% X_e = 400.2;        % Elipse center coordinate
% Y_e = 324.1000;
% X_o = 500.1000;
% Y_o= 350.1000;        % Original coordinate
% a_1 =  50.3278 ;        % elipse minor past
% b_1 =   35.3278 ;       % elipse major past