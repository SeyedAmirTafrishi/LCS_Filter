e=1;
% screen parameteres
SCREEN_X = 640;
SCREEN_Y = 480;
frame=1;
time_diff = 1/frame;%second devided by frame per sec in real activation
global ICX ICY
ICX = SCREEN_X / 2+eps;  %2
ICY = SCREEN_Y / 2+eps;  %1
En=[300 200 0 3 40 .02]
options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3]);
Vv=.4;
beta =  calculate_vector_angle( En(e,2), En(e,1), ICX, ICY)
R = (((En(e,1)-ICY)^2) + (En(e,2)-ICX)^2)^(0.5); %The R
x_0 = R;
x_1 = abs(En(e,6)+Vv)/2; % CHanged! :D
[T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options) %location of estimated E the 4 space is nutrilized to one since we want just vel
NEn(1,1) = -(ceil(Y1(end,1))-R)*sin((pi/180)*-beta) + (En(e,1)-ICY); %estimation of En x
NEn(1,2) = (ceil(Y1(end,1))-R)*cos((pi/180)*-beta) + (En(e,2)-ICX); %estimation of En y
hold on
plot(En(e,2),En(e,1),'rs')
hold on
plot(NEn(1,2) + ICX, NEn(1,1) + ICY, 'bs')
hold on
%plot( ICX, ICY, 'b*')