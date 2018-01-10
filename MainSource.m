clc
clear cam lambda psi En Er C alpha Trs Trcr Cr
cam = webcam('USB Camera'); %camera name  USB2.0 Camera USB Video Device
cam.Resolution = '640x480';
%preview(cam)
%Main intial conditions
%--------------------------- Algorithm Constants
lambda=0;
psi=0;%first is x next is y? O_
En=0;
Er=0;
C=0;
Cr=0;
delta=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
%----------------------------
alpha=[0 0 0 0 0 0 0];
global deltay deltaz Trs Trcr Trmax Av
global frame Vv
Trs=2;
Trcr=1;
Trmax=3;
delta=0;
Dv=.1;
%-------- Sensors 
m=mobiledev;
%---------
while(1)
    tic
img = snapshot(cam);%begin snap
i = img;
%Make image greyscale
if length(size(i)) == 3
	im =  double(i(:,:,2));
else
	im = double(i);
end

c9 = fast9(im, 25,1);
axis image
colormap(gray)
subplot(1,2,1)
hold on 
imshow(im / max(im(:)));
subplot(1,2,2)
hold on 
imshow(im / max(im(:)));
hold on
subplot(1,2,1)
plot(c9(:,1),c9(:,2),'r.'); %edges
hold  on
c9=[c9(:,2),c9(:,1)];
Edge=c9;
%%
%Transision V
hold  on
Edge = Line(lambda,psi,Edge);
frame=toc;
[a, t] = accellog(m);
[o, t] = orientlog(m);
Av=a(end,1)-a(end-10,1);
Vv=abs(cumtrapz(a(end-2:end,1)-a(end-10,1)));
Vv=Vv(end,1);
deltay=abs(180-abs(o(end,1)))+2;
deltaz=abs(o(end,3))+2;
[En,Er,C,Cr,psi,lambda,alpha,delta] = Circle(Edge,C,Cr,En,Er,psi,delta,Vv,Dv,lambda,alpha);
delta
hold on
subplot(1,2,1)
plot(En(:,2),En(:,1),'bs')
    xlim([1 640])
    ylim([1 480])
hold on 
th = 0:pi/50:2*pi;%for loop for creating circle
CB=1;
hold on
    if C==0
    else
    for i=1:1:(numel(C(:,1)))
    xunit = (C(i,3)+CB) * cos(th) + C(i,2);%equation of circle :D
    yunit = (C(i,3)+CB) * sin(th) + C(i,1);
    subplot(1,2,2)
    plot(xunit, yunit,'g');%Plot the boys :v 
    xlim([1 640])
    ylim([1 480])
    end
    end
hold on
CB=1;
    if Cr==0
    else
    for i=1:1:(numel(Cr(:,1)))
    xunit = (Cr(i,3)+CB) * cos(th) + Cr(i,2);%equation of circle :D
    yunit = (Cr(i,3)+CB) * sin(th) + Cr(i,1);
    subplot(1,2,2)
    plot(xunit, yunit,'r');%Plot the boys :v 
    xlim([1 640])
    ylim([1 480])
    end
    end
    size(En)
    size(Er)
clear figure
delta=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
end % dont touch while !!! 
%----------
clear(cam)