SCREEN_X = 1280;
SCREEN_Y = 960;
global ICX ICY b_config_plot_on
b_config_plot_on = true; %Ploting Graph
ICX = SCREEN_X / 2+eps;  %2
ICY = SCREEN_Y / 2+eps;  %1
MX=700;%E_n
MY=500;

% angleC = calculate_vector_angle( MX, MY, ICX, ICY );
% m = tan((pi/180)*(angleC)) % Slope of given angle of En respect to O frame
% m2= (MY-ICY)/(MX-ICX)
%             deltaT = sqrt(deltay^2 + deltax^2);
% EdgeX=710;
% EdgeY=510;
% AngularCom=((abs((-(EdgeX-ICX))+m2*(EdgeY-ICY)))/sqrt(1+m2^2))
% Dist=abs(((MX-ICX)*(ICY-EdgeY))-((ICX-EdgeX)*(MY-ICY)))/(sqrt((MX-ICX)^2+(MY-ICY)^2))
% %Acom=cos(AngularCom(1,1));
% %Bcom=sin(AngularCom(1,1));
% %AngularCom=abs((180/pi)*atan2(Bcom,Acom))  
% plot(ICX,ICY,'- *b','MarkerSize', 18,'LineWidth' , 1)
% hold on
% plot(MX,MY,'- om','MarkerSize', 18,'LineWidth' , 1)
% hold on
% plot(EdgeX,EdgeY,'- *g','MarkerSize', 18,'LineWidth' , 1)
%% Flow Optic (Field of Motion)


f=45; %Smaller means smaller distance 25 means 1 pixel
x=ICX-ICX;
y=ICY-ICY;
omx= 0.0087; %x-1
omy= -0.0627; %x-1
omz=-0.0017; 
V_x=-omy*f+omz*y+((omx*x*y)/(f))-((omy*x^2)/f)
V_y=+omx*f-omz*x-((omy*x*y)/(f))+((omx*y^2)/f)
xnew=x+V_x+ICX
ynew=y+V_y+ICY


plot(ICX,ICY,'- *b','MarkerSize', 14,'LineWidth' , 1)
hold on
plot(x+ICX,y+ICY,'- om','MarkerSize', 10,'LineWidth' , 1)
hold on
plot(xnew,ynew,'- *g','MarkerSize', 10,'LineWidth' , 1)

















