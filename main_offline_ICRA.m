%% Line-Cirle-Square Filter
% This code is an offline code with a series of captured images when
% the camera on the robot was moving in constant velocity. 

% If you use this in published work, please cite:
%       Tafrishi, S.A., Xiaotian, D., Kandjani, V.E., 2020. Line-Circle-Square (LCS): A Multilayered Geometric Filter for Edge-Based Detection,
%arXiv:2008.09315.
%     The Bibtex entries are:
%
%@misc{tafrishi2020linecirclesquare,
%    title={Line-Circle-Square (LCS): A Multilayered Geometric Filter for Edge-Based Detection},
%    author={Seyed Amir Tafrishi and Xiaotian Dai and Vahid Esmaeilzadeh Kandjani},
%    year={2020},
%    eprint={2008.09315},
%    archivePrefix={arXiv},
%    primaryClass={cs.RO}
%    }
% @inproceedings{tafrishi2017line,
%  title={Line-Circle: A Geometric Filter for Single Camera Edge-Based Object Detection},
%   author={Tafrishi, Seyed Amir and Kandjani, Vahid E},
%   booktitle={2017 5th RSI International Conference on Robotics and Mechatronics (ICRoM)},
%   pages={588--594},
%   year={2017},
%   organization={IEEE}
% }


clc
clear lambda psi En Er C alpha Trs Trcr Cr ploti
clear delta
close all
%set(0,'DefaultTextInterpreter','Latex');

%% intial conditions
% screen parameteres
SCREEN_X = 1280;
SCREEN_Y = 960;
global ICX ICY b_config_plot_on
b_config_plot_on = true; %Ploting Graph
ICX = SCREEN_X / 2+eps;  %2
ICY = SCREEN_Y / 2+eps;  %1

% algorithm constants
lambda = 0;
psi = 0;

En = 0; % Normal edges
Er = 0; % Rebel edges
C  = 0; % Normal circles
Cr = 0; % Rebel circles
S = 0;  % squares

delta = zeros(5, 4); %[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];

%% algoritm parameters
alpha = [0 0 0 0 0 0 0];
global Trs Trcr Trmax TrsSq TrcrSq TrmaxSq
global frame Av Vv deltay deltax

frame = 1; %Every Sec one frame! Works

% trust parameters
Trs = 4;   % Circle Expert Trust
Trcr = 1;  % Circle Expert Critical Trust
Trmax = 7; % Circle Expert Maximum Trust 7 better

TrsSq = 5;   % Square Expert Trust
TrcrSq = 2;  % Square Expert Critical Trust
TrmaxSq = 8; % Square Expert Maximum Trust



%% main code begins
drs = './example_pictures'; % in current directory
dr1 = dir([drs '/*.png']);  % get all png files in the folder
f1 = {dr1.name};           % get filenames to cell

%mkdir('./results')          % dir for saving results
ShiftIM=1729+145; %145 
Kts=ShiftIM; %sensor counter
% loop for each image
fcr=1730; % Frame nnumber in Image
ffr=1730; % Image number folder
Size=zeros(fcr);
for c = 1:length(f1)
   tic
c=c+ShiftIM; %390 800;
    % read one image
    i = imread([drs '/' f1{c}]);
    [pathstr,name,ext] = fileparts(f1{c});
    %BL=100; %boundery layer Cylindrical-tin

    % convert image into greyscale
    if b_config_plot_on
    figure(fcr)
    end
    if length(size(i)) == 3
        im = double(i(:,:,2));
    else
        im = double(i);
    end

   c9 = fast9(im, 137.33, 1);      % Data1 127.33 137.33
 %c9 = detectFASTFeatures(i,'MinContrast',.58); %58 Less than 58 result in many landmarks Hard to reco
%    c9 = c9.Location;

    %c9 = corner(rgb2gray(i), 'MinimumEigenvalue');
if b_config_plot_on
    axis image
    colormap(gray)

    subplot(2,2,1)
    hold on
    imshow(im / max(im(:)));
    plot(c9(:,1),c9(:,2),'r.'); % edges

    subplot(2,2,2)
    hold on
    imshow(im / max(im(:)));

    subplot(2,2,3)
    hold on
    imshow(im / max(im(:)));

        subplot(2,2,4)
    hold on
    imshow(im / max(im(:)));
end

    c9 = [c9(:,2),c9(:,1)];     % swap x and y columns
    if fcr == 1
        Size(fcr,1) = numel(c9(:,1));
    else
        Size(fcr,1) = numel(c9(:,1))+Size(fcr-1,1);
    end
    Size(fcr,6) = numel(c9(:,1));
    Edge = c9;
    
  %%  
% kinematic variables (simulated)
Dv = 0.1;
frame = toc;
 
[ggt1]= TimestampToString(str2num(name));
[ggt2]= TimestampToString(ins(Kts,1).timestamp);

while ~((minute(ggt1)<=minute(ggt2)+.01 && minute(ggt1)>=minute(ggt2)-.01)&& (hour(ggt1)<=hour(ggt2)+.005 && hour(ggt1)>=hour(ggt2)-.005)...
        && (second(ggt1)<=second(ggt2)+.01 && second(ggt1)>=second(ggt2)-.01 ))

if minute(ggt1)==minute(ggt2) && hour(ggt1)==hour(ggt2) && second(ggt1)==second(ggt2)
 
else
 Kts=Kts+1;  
[ggt1]= TimestampToString(str2num(name));
[ggt2]= TimestampToString(ins(Kts,1).timestamp);
end


end
 
if c==ShiftIM+1 %initialization of IMU
Av = abs((ins(Kts,10).velocity_north));
Vv = abs(ins(Kts,10).velocity_north);
deltay = (180/pi)*(abs(ins(Kts,15).yaw))/1+5;
deltax =(180/pi)*(abs(ins(Kts,13).roll))/1+5;
DAngz = 0; % No Rotational difference data (straight motion) in default
DAngx = 0;
DAngy=0;
Velk1= abs(ins(Kts,10).velocity_north);
yawk1=ins(Kts,15).yaw;
rollk1=ins(Kts,13).roll;  
pitch1=ins(Kts,14).pitch;
height1=ins(Kts,8).down;
else
Av = abs((ins(Kts,10).velocity_north-Velk1));
Vv = abs( ins(Kts,10).velocity_north) 
deltay = (180/pi)*(abs(ins(Kts,15).yaw-yawk1))/1+5
deltax =(180/pi)*(abs(ins(Kts,13).roll-rollk1))/1+5
DAngy= (ins(Kts,15).yaw-yawk1) % No Rotational difference data (straight motion) in default
DAngx= (ins(Kts,14).pitch-pitch1)
DAngz =(ins(Kts,13).roll-rollk1)
heightY=(ins(Kts,8).down-height1);

height1=ins(Kts,8).down;
Velk1=ins(Kts,10).velocity_north;
yawk1=ins(Kts,15).yaw;
rollk1=ins(Kts,13).roll;  
pitch1=ins(Kts,14).pitch;

end


Fcount=1;   
 %% Kinematic Rotation of Location/Angles Parameters
    %The locations, angles and origins are updated
    if (abs(DAngy)<.0015) && (abs(DAngx)<.0015) && (abs(DAngz)<.0015)% No angular rotation along x and y axes
    
    else
        if En == 0
        else
            Rx=[1 0 0;0 cos(DAngx) -sin(DAngx);0 sin(DAngx) cos(DAngx)];
            Ry=[cos(DAngy) 0 sin(DAngy);0 1 0;-sin(DAngy) 0 cos(DAngy)];
            Rz=[cos(DAngz) -sin(DAngz) 0;sin(DAngz) cos(DAngz) 0; 0 0 1];
            R=Rz*Ry*Rx;
             
            f = 70; %Smaller means smaller distance 25 means 1 pixel
 
            if max(sum(En))~=0 && (~isempty(En))
               % LEn=[En(:,2)-ICX,En(:,1)-ICY,zeros(numel(En(:,1)),1)]*R;
            x=En(:,2)-ICX;
            y=En(:,1)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            Z=1;
            %(-(heightY.*f)./(4*Z))+
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            En(:,1:2)=[ynew xnew]; %Location Update      
              % Z=1;
              %  En(:,1:2)=[(LEn(:,2)/(Z))+ICY (LEn(:,1)/(Z))+ICX]; %Location Update
            end
            if max(sum(Er))~=0 && ~isempty(Er)
            x=Er(:,2)-ICX;
            y=Er(:,1)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            Er(:,1:2)=[ynew xnew]; %Location Update  
            x=Er(:,8)-ICX;
            y=Er(:,7)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            Er(:,7:8)=[ynew xnew]; %Location Update              
%                 LEr=[Er(:,2),Er(:,1),zeros(numel(Er(:,1)),1)]*R;
%                 Er(:,1:2)=[LEr(:,2) LEr(:,1)]; %Location Update
%                 LEr=[Er(:,8),Er(:,7),zeros(numel(Er(:,1)),1)]*R;
%                 Er(:,7:8)=[LEr(:,2) LEr(:,1)];% Origin Update
            end
            if max(sum(C))~=0 && ~isempty(C)
            x=C(:,2)-ICX;
            y=C(:,1)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            C(:,1:2)=[ynew xnew]; %Location Update       
%                 LCn=[C(:,2),C(:,1),zeros(numel(C(:,1)),1)]*R;
%                 C(:,1:2)=[LCn(:,2) LCn(:,1)]; %Location Update
            end
            if max(sum(Cr))~=0 && ~isempty(Cr)
            x=Cr(:,2)-ICX;
            y=Cr(:,1)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            Cr(:,1:2)=[ynew xnew]; %Location Update  
%                 LCr=[Cr(:,2),Cr(:,1),zeros(numel(Cr(:,1)),1)]*R;
%                 Cr(:,1:2)=[LCr(:,2) LCr(:,1)]; %Location Update
            x=Cr(:,8)-ICX;
            y=Cr(:,7)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            Cr(:,7:8)=[ynew xnew]; %Location Update  
%                 LCr=[Cr(:,8),Cr(:,7),zeros(numel(Cr(:,1)),1)]*R;
%                 Cr(:,7:8)=[LCr(:,2) LCr(:,1)];% Origin Update
            end
            if max(sum(S))~=0 && ~isempty(S)
            x=S(:,2)-ICX;
            y=S(:,1)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            S(:,1:2)=[ynew xnew]; %Location Update  
%                 LS=[S(:,2),S(:,1),zeros(numel(S(:,1)),1)]*R;
%                 S(:,1:2)=[LS(:,2) LS(:,1)]; %Location Update
            x=S(:,9)-ICX;
            y=S(:,8)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            S(:,8:9)=[ynew xnew]; %Location Update  
%                 LS=[S(:,9),S(:,8),zeros(numel(S(:,1)),1)]*R;
%                 S(:,8:9)=[LS(:,2) LS(:,1)]; % Origin Update
            end
            if max(sum(alpha))~=0 && ~isempty(alpha) %Alpha (Rebel Landmarks Alignmnet)
            x=alpha(:,2)-ICX;
            y=alpha(:,1)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            alpha(:,1:2)=[ynew xnew]; %Location Update  
            x=alpha(:,4)-ICX;
            y=alpha(:,3)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            alpha(:,3:4)=[ynew xnew]; %Location Update     
            x=alpha(:,6)-ICX;
            y=alpha(:,5)-ICY;
            omx= DAngx; %x-1
            omy=-DAngy; %x-1
            omz= DAngz; 
            V_x=-omy.*f+omz.*y+((omx.*x.*y)./(f))-((omy.*x.^2)./f);
            V_y=+omx.*f-omz.*x-((omy.*x.*y)./(f))+((omx.*y.^2)./f);
            xnew=x-V_x+ICX;
            ynew=y-V_y+ICY;
            alpha(:,5:6)=[ynew xnew]; %Location Update            
%                 Lalp1=[alpha(:,2),alpha(:,1),zeros(numel(alpha(:,1)),1)]*R;
%                 alpha(:,1:2)=[Lalp1(:,2) Lalp1(:,1)];
%                 Lalp2=[alpha(:,4),alpha(:,3),zeros(numel(alpha(:,3)),1)]*R;
%                 alpha(:,3:4)=[Lalp2(:,2) Lalp2(:,1)];
%                 Lalp3=[alpha(:,6),alpha(:,5),zeros(numel(alpha(:,5)),1)]*R;
%                 alpha(:,5:6)=[Lalp3(:,2) Lalp3(:,1)];
            end
            % Update of Angles
            kc=1;
            while kc <= max([numel(En(:,1)) numel(Er(:,1)) numel(C(:,1)) numel(Cr(:,1)) numel(S(:,1))])
                if kc <= numel(En(:,1)) && max(sum(En))~=0 && (~isempty(En))
            En(kc,5) = calculate_vector_angle( En(kc,2), En(kc,1), ICX, ICY ) ;
                end
                if kc <= numel(Er(:,1)) && max(sum(Er))~=0 && (~isempty(Er))
            Er(kc,5) = calculate_vector_angle( Er(kc,2), Er(kc,1), Er(kc,8), Er(kc,7) ) ;
                end
                if kc <= numel(C(:,1)) && max(sum(C))~=0 && (~isempty(C))
            C(kc,5) = calculate_vector_angle( C(kc,2), C(kc,1), ICX, ICY ) ;
                end
                if kc <= numel(Cr(:,1)) && max(sum(Cr))~=0 && (~isempty(Cr))
            Cr(kc,5) = calculate_vector_angle( Cr(kc,2), Cr(kc,1), Cr(kc,8), Cr(kc,7) ) ;
                end
                if kc <= numel(S(:,1)) && max(sum(S))~=0 && (~isempty(S))
            S(kc,6) = calculate_vector_angle( S(kc,2), S(kc,1), S(kc,9), S(kc,8) ) ;
                end
            kc = kc + 1;
            end
        end
    end

%% Line Expert
    hold on
    Edge = Line(lambda,psi,Edge);

    
%%  Circle Expert   

    [En,Er,C,Cr,psi,lambda,alpha,delta] = Circle(Edge,C,Cr,En,Er,psi,delta,Vv,Dv,lambda,alpha);
%     if S==0
%     S=[100 100 40 30 6 60 .1 200 160;112 134 20 20 4 30 .05 100 180];
%     end
%% Square Expert
    [S, psi] = Square(S, C, Cr, delta, Vv, Dv, psi);

    % Square() add square here
    Size(fcr,2) = numel(En(:,1));
    Size(fcr,3) = numel(Er(:,1));
    Size(fcr,4) = numel(C(:,1));
    Size(fcr,5) = numel(Cr(:,1));
    Size(fcr,6) = numel(S(:,1));
  %  Size(c,9) = toc;
    if Fcount<6 % 5 Frame Sum
        Ptemp(Fcount)= numel(c9(:,1));
        Size(fcr,7) = sum(Ptemp);
        Fcount=Fcount+1;
    else
        Fcount=1;
        Ptemp(Fcount)= numel(c9(:,1));
        Size(fcr,7) = sum(Ptemp);
    end
    %***************************************************! Add step k
    %velocity to C and S and Subtract the vel. of k-1*!
    %delta
    %delta=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
  if b_config_plot_on
    hold on
    subplot(2,2,2)
    ploti = plot(En(:,2),En(:,1),'bs');
            xlim([1 SCREEN_X])
            ylim([1 SCREEN_Y])
    hold on
    th = 0:pi/50:2*pi;%for loop for creating circle
    CB = 1;
    hold on
    %plot(ICX,ICY,'- *b','MarkerSize', 18,'LineWidth' , 2.5)
    if C == 0
        % pass
    else
        for i = 1:1:(numel(C(:,1)))
            xunit = (C(i,3) + CB) * cos(th) + C(i,2);%equation of circle :D
            yunit = (C(i,3) + CB) * sin(th) + C(i,1);
            subplot(2,2,3)
            ploti = plot(xunit, yunit,'g','LineWidth' , 1.5);%Plot the boys :v
            xlim([1 SCREEN_X])
            ylim([1 SCREEN_Y])
        end
    end
    CB = 10;
    if Cr == 0
        % pass
    else
        for i = 1:1:(numel(Cr(:,1)))
            xunit = (Cr(i,3) + CB) * cos(th) + Cr(i,2);%equation of circle :D
            yunit = (Cr(i,3) + CB) * sin(th) + Cr(i,1);
            subplot(2,2,3)
            ploti = plot(xunit, yunit,'r','LineWidth' , 2);%Plot the boys :v
            xlim([1 SCREEN_X])
            ylim([1 SCREEN_Y])
        end
    end
       if S == 0
        % pass
        else
        for i = 1:1:(numel(S(:,1)))
            subplot(2,2,4)
            hold on
            TempYPositive= S(i,1)+S(i,4);  
            TempYNegaitive=S(i,1)-S(i,4);  
            TempXPositive=S(i,2)+S(i,3); %
            TempXNegaitive=S(i,2)-S(i,3); %
            %plot(X_o,Y_o,'- *b','MarkerSize', 18,'LineWidth' , 2.5)
            plot([TempXPositive TempXPositive],[TempYNegaitive TempYPositive],'m','LineWidth' , 1.3)
            plot([TempXNegaitive TempXPositive],[TempYPositive TempYPositive],'m','LineWidth' , 1.3)
            plot([TempXNegaitive TempXNegaitive],[TempYNegaitive TempYPositive],'m','LineWidth' , 1.3)
            plot([TempXNegaitive TempXPositive],[TempYNegaitive TempYNegaitive],'m','LineWidth' , 1.3)
            %ploti = plot(xunit, yunit,'r');%Plot the boys :v
            xlim([1 SCREEN_X])
            ylim([1 SCREEN_Y])
        end
       end


    hold on
    subplot(2,2,2)

    hold on
    xlabel('$\bf{E}_n, \bf{E}_r$, $\bf{\tilde{E}}_n$ and $\bf{\tilde{E}}_r$','FontSize',16,'Interpreter','latex')
    hold on
    subplot(2,2,1)
    %txt = ['Frame ',num2str(c)];
    %title(txt,'FontSize',16)
    hold on
    xlabel('\boldmath${\chi}$, \boldmath${\lambda}$ and \boldmath${\psi}$','FontSize',16,'Interpreter','latex')
    subplot(2,2,4)
    %txt = ['Frame ',num2str(c)];
    %title(txt,'FontSize',16)
    hold on
    xlabel('$\bf{S}$ and \boldmath${\psi}_S$','FontSize',16,'Interpreter','latex')
    subplot(2,2,3)
    %txt = ['Frame ',num2str(c)];
    %title(txt,'FontSize',16)
    
    hold on
    xlabel('$\bf{C}_n,\bf{C}_r$ and \boldmath${\psi}_C$','FontSize',16,'Interpreter','latex')

    subplot(2,2,1)
    hold on
    
    txt = ['Frame ',num2str(fcr)];
    text(SCREEN_X+400,SCREEN_Y+200,txt,'FontSize',16)
  %  set(gca,'OuterPosition',[0 0.15 0.31 0.7]);
    subplot(2,2,2)
    hold on
   %  set(gca,'OuterPosition',[0.35 0.15 0.29 0.7]);
    subplot(2,2,3)
    hold on
    % set(gca,'OuterPosition',[0.7 0.15 0.31 0.7]);

    set(gcf,'Units','Inches','renderer','Painters');
   % set(gcf,'Units','Inches');
    pos = get(gcf,'Position');

    %---------- Save Plot
    set(gcf, 'Position',  [100, 100, 1920, 1080])
    set(gcf, 'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(2)*3.3, pos(3)*1.3])
drawnow
    fig_filename = ['./results/fig', num2str(ffr),'.png'];
    saveas(gca, fig_filename);
    close all
    % %--------------------
    % Plot of Rebels layers
    %Crop_Image(ICX,ICY,frame,im / max(im(:)),Cr,S,fcr);
  end
  %  toc
  fcr=fcr+1;
  ffr=ffr+1;
end

clear figure
