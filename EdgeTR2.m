function dydt = EdgeTR2(t,y)
global t1 frame Vv
x1=y(1);
x2=y(2);
dx1=10; % it requires transformation (x1/frame)+ Vv worked O_O 
dx2=0;
dydt= [dx1;dx2];

