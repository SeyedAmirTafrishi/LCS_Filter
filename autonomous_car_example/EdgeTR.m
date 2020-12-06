function dydt = EdgeTR(t,y)
global t1 frame Vv Av
x1 = y(1);
x2 = y(2);
dx1 =1.2*Vv; % it requires transformation (x1/frame)+ Vv worked O_O  280
dx2 = Av;
dydt = [dx1;dx2];
