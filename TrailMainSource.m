clc
clear all
clear figure
%Main intial conditions
lambda=[240 140 200;500 90 80];
%psi=[150 180 120 1];%first is x next is y? O_O
psi=0;
En=[110 550 200 2 60.5241 1;158 225 100 4 130.7994 .3];
EOLD=En;
Er=0;
C=0;
alpha=0;
global deltay deltaz Trs Trcr
deltay=40;
deltaz=40;
Trs=5;
Trcr=2;
delta=0;
Vv=.01;
Dv=.1;
i = imread('g1.png');
%Make image greyscale
if length(size(i)) == 3
	im =  double(i(:,:,2));
else
	im = double(i);
end
c9 = fast9(im, 25,1);
axis image
colormap(gray)

imshow(im / max(im(:)));
hold on
plot(c9(:,1),c9(:,2),'r.'); %edges
Edge=c9;
%%

Edge = Line(lambda,psi,Edge);
[En,Er,C,psi,lambda,alpha] = Circle(Edge,C,En,Er,psi,delta,Vv,Dv,lambda,alpha);
plot(En(:,2),En(:,1),'bs')
hold on 
plot(EOLD(:,2),EOLD(:,1),'ys')%So Z is second and Y is first 
%----------
title('9 point FAST');
%clear(cam)