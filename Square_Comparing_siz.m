clc
clear


%Square 1
a_1=5;
b_1=10;
Lx_1=50;
Ly_1=40;


%Square 2
a_2=11;
b_2=14;
Lx_2=50;
Ly_2=30;

%Find where Sq 1 is located around 2
if (Lx_1+a_1<Lx_2-a_2 || Lx_1-a_1>Lx_2+a_2 ) ||  (Ly_1-b_1>Ly_2+b_2 || Ly_1+b_1<Ly_2-b_2 )% Sq_1 is out of Sq_2
 disp("not in") % If out then good bye
    
else %Always in
 disp("is in") % if in who is bigger?
 % Find intersecting area!
 %need to know direction
 
end
hold on
plot(Lx_1,Ly_1,'- *r','MarkerSize', 18,'LineWidth' , 2.5)
hold on
plot([Lx_1+a_1 Lx_1+a_1],[Ly_1-b_1 Ly_1+b_1],'b')
plot([Lx_1-a_1 Lx_1+a_1],[Ly_1+b_1 Ly_1+b_1],'b')
plot([Lx_1-a_1 Lx_1-a_1],[Ly_1-b_1 Ly_1+b_1],'b')
plot([Lx_1-a_1 Lx_1+a_1],[Ly_1-b_1 Ly_1-b_1],'b')


hold on
plot(Lx_2,Ly_2,'- *g','MarkerSize', 18,'LineWidth' , 2.5)
hold on
plot([Lx_2+a_2 Lx_2+a_2],[Ly_2-b_2 Ly_2+b_2],'m')
plot([Lx_2-a_2 Lx_2+a_2],[Ly_2+b_2 Ly_2+b_2],'m')
plot([Lx_2-a_2 Lx_2-a_2],[Ly_2-b_2 Ly_2+b_2],'m')
plot([Lx_2-a_2 Lx_2+a_2],[Ly_2-b_2 Ly_2-b_2],'m')