function F = FindTangentAlgebric( REF)
X_e = REF(1,1); %Elipse center coordinate
Y_e = REF(1,2);
X_o = REF(1,3);
Y_o = REF(1,4); %Original coordinate
a_1 = REF(1,5); %elipse minor past
b_1 = REF(1,6); %elipse major past

F(1) = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
C=(b^2/a^2)*((X_o-F(1))*(F(1)-X_e))+Y_o*Y_E;
A=roots([1 -(Y_o+Y_e) C])
F(2)=0
F(3)=0;
F(4)=0;


end
