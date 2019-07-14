%%
X_e= 1;%Elipse center coordinate
Y_e= 3;
Y_o=1;%Original coordinate
X_o=3;
a_1=1; %elipse minor past
b_1=.5; %elipse major past
syms X Y

F1=(a_1^2*(Y_o-Y)*(Y-Y_e))+(b_1^2*(X_o-X)*(X-X_e))==0;
F2=((Y-Y_e)^2/(b_1^2))+((X-X_e)^2/(a_1^2))==1;

sol = solve([F1, F2], [X, Y]);
Point1 = double(sol.X)
Point2 = double(sol.Y)
Px1=Point1(1,1);
Py1=Point2(1,1);
Px2=Point1(2,1);
Py2=Point2(2,1);

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
 %%   R_e R_o
 R_e=sqrt((Y_e-Y_o)^2+(X_e-X_o)^2) %Orignal to Elipse
 R_o1=sqrt((Point2(1,1)-Y_o)^2+(Point1(1,1)-X_o)^2) %Orignal to contact 1
 R_o2=sqrt((Point2(2,1)-Y_o)^2+(Point1(2,1)-X_o)^2) %Orignal to contact 2
 %for two different side!!! trignomatic
 if R_e>R_o1 %Hypto
 alpha1=acos(R_o1/R_e)
 else
 alpha1=acos(R_e/R_o1)   
 end
 if R_e>R_o2
 alpha2=acos(R_o2/R_e)    
 else
 alpha2=acos(R_e/R_o2)      
 end

 
 
 
 
 
 