function F = FindTangentx1(x,REF)
X_e=REF(1,1);%Elipse center coordinate
Y_e=REF(1,2);
X_o=REF(1,3);
Y_o=REF(1,4);%Original coordinate
a_1=REF(1,5); %elipse minor past
b_1=REF(1,6); %elipse major past

%x_2=Y_e-((b_1/a_1)*sqrt((a_1)^2-(x(1)-X_e)^2))
% x_1=X_e-((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2))
% F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
%F(2) = ((x(2)-Y_e)^2/(b_1^2))+((x(1)-X_e)^2/(a_1^2))-1;
% %Solution 1 
% x_1=X_e+((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2))
% F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
%Solution 2
if (X_o>X_e && Y_o>Y_e) || (X_o<X_e && Y_o>Y_e) %+ + 
x_1=X_e+((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2));
%save x_1
F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
y_2=Y_e+((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%save y_2
F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e))+(b_1^2*(X_o-x(2))*(x(2)-X_e));  

else %- +
x_1=X_e+((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%save x_1
F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
y_2=Y_e-((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%save y_2
F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e))+(b_1^2*(X_o-x(2))*(x(2)-X_e));    

%elseif if X_o=X_e    
end
 
end


