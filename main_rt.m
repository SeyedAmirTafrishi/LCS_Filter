%% Line-Cirle-Square Filter
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

%--------------------------- Camera Parameters
CAM_WIDTH = 640;
CAM_HEIGHT = 480;

%--------------------------- Open the Camera
cam = webcam(1);
% or specify the camera name, for example:
%cam = webcam('USB Camera'); %camera name  USB2.0 Camera USB Video Device
%cam = webcam('Logitech HD Pro Webcam C920');

cam.Resolution = sprintf('%dx%d', CAM_WIDTH, CAM_HEIGHT);
% You can preview the camera with:
%preview(cam)

% Main intial conditions
%--------------------------- Algorithm Constants
lambda = 0;
psi = 0;

En = 0; % Normal edges
Er = 0; % Rebel edges
C  = 0; % Normal circles
Cr = 0; % Rebel circles
S = 0;  % squares

delta = zeros(5, 4); % [0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];

%----------------------------
alpha = [0 0 0 0 0 0 0];
global Trs Trcr Trmax TrsSq TrcrSq TrmaxSq
global frame Av Vv deltay deltax b_config_plot_on

b_config_plot_on = true; % Enable ploting graph

delta = 0;
Dv = .1;

Trs = 3;   % Circle Expert Trust
Trcr = 2;  % Circle Expert Critical Trust
Trmax = 5; % Circle Expert Maximum Trust

TrsSq = 4;   % Square Expert Trust
TrcrSq = 3;  % Square Expert Critical Trust
TrmaxSq = 6; % Square Expert Maximum Trust

SCREEN_X = 640; % Image width
SCREEN_Y = 480; % Image height

global ICX ICY
ICX = SCREEN_X / 2 + eps;
ICY = SCREEN_Y / 2 + eps;

%-------- Sensors
%m = mobiledev;     % uncomment this if a real IMU is used!
pause(5);           % pause for a few seconds for data to be initialized

%---------
while(1)
    tic
    img = snapshot(cam); % begin snap
    i = img;
    % Make image greyscale
    if length(size(i)) == 3
        im = double(i(:,:,2));
    else
        im = double(i);
    end

    c9 = detectFASTFeatures(rgb2gray(i),'MinContrast',0.2);
    % Alternatively, use the fast9 edge detector:
    %c9 = fast9(im, 30, 1);      % run fast9 edge detection
    c9 = c9.Location;

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

%------ Trail with IMU
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
    DAngy = 0; % No Rotational difference data (straight motion) in default
    DAngx = 0;


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
        subplot(2,2,1)
        hold on
        xlabel('${\chi}, {\lambda}, \psi$','FontSize',16,'Interpreter','latex')

        subplot(2,2,2)
        hold on
        xlabel('${E}_n, {E}_r$, $\tilde{E}_n$ and $\tilde{E}_r$','FontSize',16,'Interpreter','latex')
        hold on

        subplot(2,2,3)
        hold on
        xlabel('${C}_n,{C}_r, \psi_C$','FontSize',16,'Interpreter','latex')

        subplot(2,2,4)
        hold on
        xlabel('${S},\psi_S$','FontSize',16,'Interpreter','latex')
    end
        clear figure
        delta = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
end % end of while: dont touch while !!!

%----------
% kill the camera after the code is finished
clear(cam)
