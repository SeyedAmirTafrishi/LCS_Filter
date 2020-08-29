%% 
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


% Notes: 
% 1. To use camera in this filter, please install the "MATLAB Support
% Package for USB Webcams".
% 2. The current filter uses the constant values for angular errors
% (deltay and deltax), if you want to change/improve these default values
% please refere to IMU section. Also, the velocity (Vv) and 
% acceleration (Av) are assumed constant. One can change it to desired 
% values or use a real-time IMU sensor. 
 

clc
clear m cam lambda psi En Er C alpha Trs Trcr Cr
%addpath('./helpers/')

%--------------------------- Parameters
CAM_WIDTH = 640;
CAM_HEIGHT = 480;

%--------------------------- Open Camera
%cam = webcam('USB Camera'); %camera name  USB2.0 Camera USB Video Device
%cam = webcam('Logitech HD Pro Webcam C920');
cam = webcam(1);
cam.Resolution = sprintf('%dx%d', CAM_WIDTH, CAM_HEIGHT);
%preview(cam)

% Main intial conditions
%--------------------------- Algorithm Constants
lambda = 0;
psi = 0; %first is x next is y? O_
En = 0;
Er = 0;
C  = 0;
Cr = 0;
S = 0;
delta = zeros(5, 4); %[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
%----------------------------
alpha = [0 0 0 0 0 0 0];
global Trs Trcr Trmax TrsSq TrcrSq TrmaxSq
global frame Av Vv deltay deltax b_config_plot_on
b_config_plot_on = true; %Ploting Graph
Trs = 3; % Circle expert's standard trust factor
Trcr = 2; % Circle expert's critical trust factor
Trmax = 5; % Circle expert's maximum trust factor
delta = 0;
Dv = .1;
Trs   = 3; %Circle Expert Trust
Trcr  = 2; %Circle Expert Critical Trust
Trmax = 5; % Circle Expert Maximum Trust

TrsSq=4; %Square Expert Trust
TrcrSq=3; %Square Expert Critical Trust
TrmaxSq=6; % Square Expert Maximum Trust
SCREEN_X = 640; %Image Dim.
SCREEN_Y = 480; % Image Dim.

global ICX ICY
ICX = SCREEN_X / 2+eps;  %2
ICY = SCREEN_Y / 2+eps;  %1
%-------- Sensors
%m = mobiledev;
pause(5);

%---------
while(1)
    tic
    img = snapshot(cam);% begin snap
    i = img;
    % Make image greyscale
    if length(size(i)) == 3
        im = double(i(:,:,2));
    else
        im = double(i);
    end

    %c9 = fast9(im, 30, 1);      % run fast9 edge detection
    c9 = detectFASTFeatures(rgb2gray(i),'MinContrast',0.2);
    c9 = c9.Location;
    %c9 = corner(rgb2gray(i), 'MinimumEigenvalue');

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


    c9 = [c9(:,2),c9(:,1)];     % swap x and y columns

    Edge = c9;


  %% IMU Sensor Data
%  Instructions: 1. deltay and deltax are parameters that can be found in
%  different ways as follows:
% Option A: Consecetive experiments can be taken place for determining the 
% angular errors (deltay and deltax) that occur. It will be the difference
% between the desired real angles with respect to the obtained IMU angular
% data. 
% Option B: The angular errors (deltay and deltax) can be changing values 
% which are functions of the covariance matrix if the IMU is observed 
% through the Kalman Filter. 
% 2. When there is no IMU sensor, the filter only works for straight motion
% of the camera; Hence, the angular rotations (DAngy and DAngx) are assumed
% approximately zero. 

%-------Trail With IMU
%     [a, t] = accellog(m);
%     [o, t] = orientlog(m);
%     Av = a(end,1) - a(end-10,1);
%     Vv = abs(cumtrapz(a(end-2:end,1) - a(end-10,1)));
%     Vv = Vv(end,1);
%     deltay = abs(180-abs(o(end,1))) + 2;
%     deltax = abs(o(end,3)) + 2;
% 
%     wf = 1;
%     wm = 1;
%     quest = QUEST(fb, mb, fn, mn, wf, wm); % Quest Filter Observer (IMU)
%     [Quest1, Quest2, Quest3] = ... %Angles
%DAngy=Quest2; % Angular Rotation 
%DAngx=Quest1; % Angular Rotation 
%------ Trail without IMU
    Av = 0.0005; % Accelertion in Direction of motion
    Vv = .03; % Velocity in direction of motion
    deltay = 9; % Calculated Angular error along y axis
    deltax = 9; % Calculated Angular error along x axis
    DAngy=0; % No Rotational difference data (straight motion) in default 
    DAngx=0;
    
 %% Kinematic Rotation of Location/Angles Parameters
    %The locations, angles and origins are updated 
    if DAngy<.1 && DAngx<.1 % No angular rotation along x and y axes
    else
    if En==0       
    else
Rx=[1 0 0;0 cos(DAngx) -sin(DAngx);0 sin(DAngx) cos(DAngx)];
Ry=[cos(DAngy) 0 sin(DAngy);0 1 0;-sin(DAngy) 0 cos(DAngy)];
R=Ry*Rx;

if max(sum(En))~=0 && (~isempty(En))
LEn=[En(:,2),En(:,1),zeros(numel(En(:,1)),1)]*R;
En(:,1:2)=[LEn(:,2) LEn(:,1)]; %Location Update
end
if max(sum(Er))~=0 && ~isempty(Er) 
LEr=[Er(:,2),Er(:,1),zeros(numel(Er(:,1)),1)]*R;
Er(:,1:2)=[LEr(:,2) LEr(:,1)]; %Location Update
LEr=[Er(:,8),Er(:,7),zeros(numel(Er(:,1)),1)]*R;
Er(:,7:8)=[LEr(:,2) LEr(:,1)];% Origin Update
end
if max(sum(C))~=0 && ~isempty(C) 
LCn=[C(:,2),C(:,1),zeros(numel(C(:,1)),1)]*R;
C(:,1:2)=[LCn(:,2) LCn(:,1)]; %Location Update
end
if max(sum(Cr))~=0 && ~isempty(Cr) 
LCr=[Cr(:,2),Cr(:,1),zeros(numel(Cr(:,1)),1)]*R;
Cr(:,1:2)=[LCr(:,2) LCr(:,1)]; %Location Update
LCr=[Cr(:,8),Cr(:,7),zeros(numel(Cr(:,1)),1)]*R;
Cr(:,7:8)=[LCr(:,2) LCr(:,1)];% Origin Update
end
if max(sum(S))~=0 && ~isempty(S) 
LS=[S(:,2),S(:,1),zeros(numel(S(:,1)),1)]*R;
S(:,1:2)=[LS(:,2) LS(:,1)]; %Location Update
LS=[S(:,9),S(:,8),zeros(numel(S(:,1)),1)]*R;
S(:,8:9)=[LS(:,2) LS(:,1)]; % Origin Update
end
if max(sum(alpha))~=0 && ~isempty(alpha) %Alpha (Rebel Landmarks Alignmnet)
Lalp1=[alpha(:,2),alpha(:,1),zeros(numel(alpha(:,1)),1)]*R;
alpha(:,1:2)=[Lalp1(:,2) Lalp1(:,1)];
Lalp2=[alpha(:,4),alpha(:,3),zeros(numel(alpha(:,3)),1)]*R;
alpha(:,3:4)=[Lalp2(:,2) Lalp2(:,1)];
Lalp3=[alpha(:,6),alpha(:,5),zeros(numel(alpha(:,5)),1)]*R;
alpha(:,5:6)=[Lalp3(:,2) Lalp3(:,1)];
end
%Update of Angles
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
kc=kc+1;
end
    end
    end
%% Line Expert
    hold on
    Edge = Line(lambda,psi,Edge);
    frame = toc;

%% Circle Expert 
    [En,Er,C,Cr,psi,lambda,alpha,delta] = Circle(Edge,C,Cr,En,Er,psi,delta,Vv,Dv,lambda,alpha);
 
    
%% Square Expert     
    [S, psi] = Square(S, C, Cr, delta, Vv, Dv, psi);
    
%% Plotting
if b_config_plot_on  
    hold on
    subplot(2,2,2)
    ploti = plot(En(:,2),En(:,1),'bs');
    xlim([1 640])
    ylim([1 480])
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
            ploti = plot(xunit, yunit,'g');%Plot the boys :v
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
            ploti = plot(xunit, yunit,'r');%Plot the boys :v
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
    TempYPositive= S(i,1)+S(i,3); %
    TempYNegaitive=S(i,1)-S(i,3); %
    TempXPositive=S(i,2)+S(i,4); %
    TempXNegaitive=S(i,2)-S(i,4); %
    %plot(X_o,Y_o,'- *b','MarkerSize', 18,'LineWidth' , 2.5)
    plot([TempXPositive TempXPositive],[TempYNegaitive TempYPositive],'m','LineWidth' , 2)
    hold on
    plot([TempXNegaitive TempXPositive],[TempYPositive TempYPositive],'m','LineWidth' , 2)
    plot([TempXNegaitive TempXNegaitive],[TempYNegaitive TempYPositive],'m','LineWidth' , 2)
    plot([TempXNegaitive TempXPositive],[TempYNegaitive TempYNegaitive],'m','LineWidth' , 2)
    hold on
            %ploti = plot(xunit, yunit,'r');%Plot the boys :v
            xlim([1 SCREEN_X])
            ylim([1 SCREEN_Y])
        end
       end
    drawnow


    hold on
    subplot(2,2,2)
    hold on
    xlabel('${E}_n, {E}_r$, $\tilde{E}_n$ and $\tilde{E}_r$','FontSize',16,'Interpreter','latex')
    hold on
    subplot(2,2,1)
    %txt = ['Frame ',num2str(c)];
    %title(txt,'FontSize',16)
    hold on
    xlabel('${\chi}, {\lambda}, \psi$','FontSize',16,'Interpreter','latex')
    subplot(2,2,4)
    %txt = ['Frame ',num2str(c)];
    %title(txt,'FontSize',16)
    hold on
    xlabel('${S},\psi_S$','FontSize',16,'Interpreter','latex')
    subplot(2,2,3)
    %txt = ['Frame ',num2str(c)];
    %title(txt,'FontSize',16)
    hold on
    xlabel('${C}_n,{C}_r, \psi_C$','FontSize',16,'Interpreter','latex')
end
    clear figure
    delta = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
end % end of while: dont touch while !!!

%----------
clear(cam)


%
%     hold on
%     subplot(1,2,1)
%     plot(En(:,2),En(:,1),'bs')
%     xlim([1 CAM_WIDTH])
%     ylim([1 CAM_HEIGHT])
%     hold on
%     th = 0:pi/50:2*pi;% for loop for creating circle
%     CB = 1;
%     hold on
%     if C == 0
%         % pass
%     else
%         for i = 1:1:(numel(C(:,1)))
%             xunit = (C(i,3) + CB) * cos(th) + C(i,2);% equation of circle :D
%             yunit = (C(i,3) + CB) * sin(th) + C(i,1);
%             subplot(1,2,2)
%             plot(xunit, yunit,'g');% Plot the boys :v
%             xlim([1 CAM_WIDTH])
%             ylim([1 CAM_HEIGHT])
%         end
%     end
%     hold on
%     CB = 1;
%     if Cr == 0
%         % pass
%     else
%         for i = 1:1:(numel(Cr(:,1)))
%             xunit = (Cr(i,3)+CB) * cos(th) + Cr(i,2);% equation of circle :D
%             yunit = (Cr(i,3)+CB) * sin(th) + Cr(i,1);
%             subplot(1,2,2)
%             plot(xunit, yunit, 'r');% Plot the boys :v
%             xlim([1 CAM_WIDTH])
%             ylim([1 CAM_HEIGHT])
%         end
%     end
%     size(En)
%     size(Er)
