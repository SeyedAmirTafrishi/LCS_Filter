clc
clear lambda psi En Er C alpha Trs Trcr Cr ploti
clear delta
set(0,'DefaultTextInterpreter','Latex');
%preview(cam)
%Main intial conditions
%--------------------------- Algorithm Constants
lambda = 0;
psi = 0;%first is x next is y? O_
En = 0;
Er = 0;
C = 0;
Cr = 0;
delta = [0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0];
%----------------------------
alpha = [0 0 0 0 0 0 0];
global deltay deltaz Trs Trcr Trmax Av
global frame Vv
Trs = 3;
Trcr = 2;
Trmax = 5;
Dv = 0.05;
Av = 0.001;
Vv = 0.05;
deltay = 8;
deltaz = 8;
%----image series
drs = './Approaching_Boxi'; % in current directory
dr1 = dir([drs '/*.png']); % get only jpg
f1 = {dr1.name};% get only filenames to cell

for c = 1:length(f1) % for each image
  tic
i = imread([drs '/' f1{c}]);
%image series
%BL=100; %boundery layer Cylindrical-tin
%img = snapshot(cam);%begin snap
%
%i = img;
%Make image greyscale
figure(c)

%subplot(length(f1)/4,length(f1)/5,c)
%Make image greyscale
if length(size(i)) == 3
  im =  double(i(:,:,2));
else
  im = double(i);
end

c9 = fast9(im, 25, 1);
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
c9 = [c9(:,2),c9(:,1)];
if c == 1
  Size(c,1) = numel(c9(:,1));
else
  Size(c,1) = numel(c9(:,1))+Size(c-1,1);
end
Size(c,6) = numel(c9(:,1));
Edge = c9;
%--------Algo begins HERE ......!!!!!
Edge = Line(lambda,psi,Edge);
frame = 1; %Every Sec one frame! Works
[En,Er,C,Cr,psi,lambda,alpha,delta] = Circle(Edge,C,Cr,En,Er,psi,delta,Vv,Dv,lambda,alpha);
Size(c,2) = numel(En(:,1));
Size(c,3) = numel(Er(:,1));
Size(c,4) = numel(C(:,1));
Size(c,5) = numel(Cr(:,1));
%delta
%delta=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
hold on
subplot(1,2,1)
ploti = plot(En(:,2),En(:,1),'bs');
xlim([1 640])
ylim([1 480])
hold on
th = 0:pi/50:2*pi;%for loop for creating circle
CB=1;
hold on
  if C == 0
  else
    for i=1:1:(numel(C(:,1)))
      xunit = (C(i,3)+CB) * cos(th) + C(i,2);%equation of circle :D
      yunit = (C(i,3)+CB) * sin(th) + C(i,1);
      subplot(1,2,2)
      ploti = plot(xunit, yunit,'g');%Plot the boys :v
      xlim([1 640])
      ylim([1 480])
    end
  end
  CB = 1;
  if Cr == 0
  else
    for i=1:1:(numel(Cr(:,1)))
      xunit = (Cr(i,3)+CB) * cos(th) + Cr(i,2);%equation of circle :D
      yunit = (Cr(i,3)+CB) * sin(th) + Cr(i,1);
      subplot(1,2,2)
      ploti = plot(xunit, yunit,'r');%Plot the boys :v
      xlim([1 640])
      ylim([1 480])
    end
  end
  hold on
subplot(1,2,2)
txt = ['Frame ',num2str(c)];
title(txt,'FontSize',16)
hold on
xlabel('$C_n$, $E_r$ and $C_r$','FontSize',16)
hold on
subplot(1,2,1)
txt = ['Frame ',num2str(c)];
title(txt,'FontSize',16)
hold on
xlabel('$E_n$, $Edge$, $\lambda$ and $\psi$','FontSize',16)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%----------Plot
 temp = ['fig',num2str(c),'.fig'];
 saveas(gca,temp);
% %--------------------
end
clear figure
%clear(cam)
