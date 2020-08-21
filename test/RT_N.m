% SCREEN_X = 640;
% SCREEN_Y = 480;
% fu=.1;
% fv=.1;
% s=3;
% OIz = SCREEN_X / 2+eps;  %2
% OIy = SCREEN_Y / 2+eps;  %1
% DAngy=0; %angular differences
% DAngz=0.4;
% 
% t=[0;0;0];
% Dv=0.03;
% Lpk1=[30;30]/Dv; %Inverse Depth
% 
% Rz=[cos(DAngz) -sin(DAngz) 0;sin(DAngz) cos(DAngz) 0;0 0 1];
% Ry=[cos(DAngy) 0 sin(DAngy);0 1 0;-sin(DAngy) 0 cos(DAngy)];
% R=Rz*Ry;
% L=[s*R s*t;zeros(1,3) 1];
% K=[Dv;Lpk1]/norm([Dv;Lpk1])
% y=inv(L)*[K;1] % Environment


SCREEN_X = 640;
SCREEN_Y = 480;
fu=.1;
fv=.1;
s=10; %Camera ratio
OIx = SCREEN_X / 2+eps;  %2
OIy = SCREEN_Y / 2+eps;  %1
DAngy=.3; %angular differences
DAngz=.5;

t=[0;0;0];
Dv=0.03;
Lpk1=[30;30]; %Inverse Depth
LPn=[Lpk1(1)-OIx;Lpk1(2)-OIy];
Rz=[cos(DAngz) -sin(DAngz) 0;sin(DAngz) cos(DAngz) 0;0 0 1];
Ry=[cos(DAngy) 0 sin(DAngy);0 1 0;-sin(DAngy) 0 cos(DAngy)];
R=Rz*Ry;
L=R*[s*LPn;s*.15] %we have to add image ratio! 
%adding velocity integral error! 
 
Lx=fu*(L(1)/L(3))+OIx
Ly=fv*(L(2)/L(3))+OIy


hold on
plot(L(1)+OIx,L(2)+OIy,'bs','LineWidth',3);
hold on 
plot(Lpk1,Lpk1,'ms');
hold on
plot(OIx,OIy,'b*');
hold on
plot(Lx,Ly,'g*');




