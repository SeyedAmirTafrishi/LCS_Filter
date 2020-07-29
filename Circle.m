function [En,Er,C,Cr,psi,lambda,alpha,delta] = Circle(Edge,C,Cr,En,Er,psi,delta,Vv,Dv,lambda,alpha)
%%
set(0,'DefaultTextInterpreter','Latex');
%---General
global beta
global deltay deltaz Trs Trcr Trmax ploti
global time_diff frame
global ICX ICY

% VeS=2;%Standard Edge Velocity
BLS = 25; %Standard Boundery layer Initialization Constant
L = 3; %The el of rebel edge alignment, Step of accuracy to have rebel edges
options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3]);

NormRows = sqrt(sum(Edge.*Edge,2));
EdgeNorm = bsxfun(@rdivide,abs(Edge),NormRows);
time_diff = 1/frame;%second devided by frame per sec in real activation

%% normal edge
%----------En
k = 1; %counter of En %%REMOVE
if En == 0 %CHANGE!
    %REMOVE
else
    if lambda == 0
        %pass
    else
        %Edge couple counter
        e = 1;
        while e <= numel(En(:,1))%previous steps edges
            ME1 = 0;
            ME2 = 0;
            beta = En(e,5);
            R = (((En(e,1)-ICY)^2) + (En(e,2)-ICX)^2)^(0.5); %The R
            x_0 = R;
            x_1 = abs(En(e,6)+Vv)/2; % CHanged! :D
            [T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated E the 4 space is nutrilized to one since we want just vel
            NEn(1,1) = -(ceil((Y1(end,1)))-R)*sin((pi/180)*beta) + (En(e,1)-ICY); %estimation of En x ceil
            NEn(1,2) = (ceil((Y1(end,1)))-R)*cos((pi/180)*beta) + (En(e,2)-ICX); %estimation of En y
            hold on
            subplot(2,2,2)
            plot(NEn(1,2) + ICX, NEn(1,1) + ICY, 'ys','LineWidth' , 2)
            %hold on
            %plot(NEn(1,2)+ICX,NEn(1,1)+ICY,'y*')
            NVe = Y1(end,2); % Estimated edge velocity
            m = tan((pi/180)*(beta + 90)); % Slope of given angle of En respect to O frame
            deltaT = sqrt(deltay^2 + deltaz^2);
            z = 1; % loop search of suitable edge in Lambda
            FlagleftE=0; %Flag for left out edge 
            
            while z <= (numel(lambda(:,1))) % Lambda column counter
                if ( ((((lambda(z,1)-En(e,1))^2) + ((lambda(z,2)-En(e,2))^2)) ^ 0.5) <= lambda(z,3) ) % this is Lambda Check in En ok?
                    i = 1; %Check the Edge to find related group
                    while (i <= (numel(Edge(1,:)))) %finder of lambda and Edge Match / row counter
                        j = 1;
                        while (j <= (numel(Edge(:,1))))%Column counter
                            if ~(Edge(j,i)==0 && Edge(j,i+1)==0) && (((((lambda(z,1)-Edge(j,i))^2) + ((lambda(z,2)-Edge(j,i+1))^2))^(0.5)) <= lambda(z,3)) % find which Edge is in circle of Lambda
                                %MainEdge(1,1)=Edge(j,i);%no need for this?
                                %MainEdge(1,2)=Edge(j,i+1);
                                %plot(MainEdge(1,2),MainEdge(1,1),'rs')
                                ME1 = j; % Columns of encompassed grouped Edge in Lambda
                                ME2 = i;
                                i = numel(Edge(1,:)) + 1; % if there is a match, break the loop
                                j = numel(Edge(:,1)) + 1;
                                z = numel(lambda(:,1)) + 1;
                                 FlagleftE=1;
                            end
                            j = j + 1;
                        end
                        i = i + 2;
                    end
                end
                z = z + 1;
            end
            if FlagleftE==0 %NO MATCH E_n
            En(e,:)=[];
            if e<2
             e=1;   
            else
             e=e-1;  
            end
            else
            if (ME2==0 && ME1==0) % lonely not a single match?  ([issue] Maybe can be moved)
                En(e,1) = NEn(1,1)+ICY; %This needed a frame center shift (added)
                En(e,2) = NEn(1,2)+ICX;
                En(e,4) = En(e,4) - 1;
            else % Main En and Lambda matching (Lambda 1 3 ...)
                NBL1 = ( (abs(Vv-NVe) / det(corr(EdgeNorm(:,ME2:ME2+1)))) + En(e,3) ) / 2; % Estmated Boundary with using Eq. (4.1) (En and Edge Group)
                NBL = (NBL1 + En(e,3)) / 2;
                j = 1;
                Edgetrans(1,1) = 0; % The nearest Edge to En finder (temporary variable for comparing)
                Edgetrans(1,2) = 0;
                d = -1;   % d is a distance value, if d = -1 that means deactive (En is lambda2,3)
                side = 0; % flag for choosing either Lambda 1 or 2
                while ( j <= (numel(Edge(:,1))) ) % loop check for classifications ME1, Edge,
                    % --- Lambda 2 and 3 condition check
                    if (((((NEn(1,1)-(Edge(j,ME2)-ICY))^2) + ((NEn(1,2)-(Edge(j,ME2+1)-ICX))^2))^(0.5)) <= NBL) && ~(Edge(j,ME2)==0 && Edge(j,ME2+1)==0) && (((abs((-(Edge(j,ME2+1)-ICX))+m*(Edge(j,ME2)-ICY)))/sqrt(1+m^2)) < deltaT) % Circle with NBL size AND deltay delta x error and zero remover,  yEdge=Edge(j,ME2) first parameter
                        if ((NEn(1,1)/abs(NEn(1,1))) == ((Edge(j,ME2)-ICY)/abs((Edge(j,ME2)-ICY)))  && ~(NEn(1,2)/abs(NEn(1,2))) == ((Edge(j,ME2+1)-ICX)/abs((Edge(j,ME2+1)-ICX)))) && (((NEn(1,1)^2+NEn(1,2)^2)^(0.5)) <= (((Edge(j,ME2)-ICY)^2+(Edge(j,ME2+1)-ICX)^2)^(0.5))) %Sign check and magnitude check which is far, lambda_2
                            En(e,1) = (((En(e,4)-Trcr)*(NEn(1,1)+ICY)) + Edge(j,ME2)) / ((En(e,4)-Trcr) + 1);
                            En(e,2) = (((En(e,4)-Trcr)*(NEn(1,2)+ICX)) + Edge(j,ME2+1)) / ((En(e,4)-Trcr) + 1); %Estimation of En, X direction
                            En(e,3) = NBL;
                            %subplot(1,2,1)
                            %hold on
                            %plot(En(e,2),En(e,1),'r*')
                            if (((Edgetrans(1,1)-ICY)^2+(Edgetrans(1,2)-ICX)^2)^(0.5) > (((Edge(j,ME2)-ICY)^2+(Edge(j,ME2+1)-ICX)^2)^(0.5))) && (~(Edgetrans(1,1)==0 && Edgetrans(1,2)==0))%Edge Remmover and modifier
                                d = (abs((-(Edge(j,ME2+1)-ICX)) + m*(Edge(j,ME2)-ICY))) / sqrt(1 + m^2);
                                a = Edgetrans(1,1);
                                b = Edgetrans(1,2);
                                Edgetrans(1,1) = Edge(j,ME2); % Update the nearest point to the estimated En
                                Edgetrans(1,2) = Edge(j,ME2+1);
                                Edge(j,ME2) = a; % return the worst match to the edge again
                                Edge(j,ME2+1) = b;
                                side = 1;
                            elseif (Edgetrans(1,1) == 0 && Edgetrans(1,2) == 0)
                                d = (abs((-(Edge(j,ME2+1)-ICX))+m*(Edge(j,ME2)-ICY)))/sqrt(1+m^2);
                                Edgetrans(1,1) = Edge(j,ME2); % The nearest finder
                                Edgetrans(1,2) = Edge(j,ME2+1);
                                Edge(j,ME2) = 0;
                                Edge(j,ME2+1) = 0;
                                En(e,4) = En(e,4) + 1; % CHECK THIS? right?
                                side = 1; % front
                            end
                        elseif ((NEn(1,1)/abs(NEn(1,1))) == ((Edge(j,ME2)-ICY)/abs((Edge(j,ME2)-ICY)))  && (NEn(1,2)/abs(NEn(1,2))) == ((Edge(j,ME2+1)-ICX)/abs((Edge(j,ME2+1)-ICX)))) && (((NEn(1,1)^2+NEn(1,2)^2)^(0.5)) >= (((Edge(j,ME2)-ICY)^2+(Edge(j,ME2+1)-ICX)^2)^(0.5)))%lambda*_3
                            En(e,1) = ((((En(e,4)-Trcr)*(NEn(1,1)+ICY))+Edge(j,ME2))/((En(e,4)-Trcr)+1));
                            En(e,2) = ((((En(e,4)-Trcr)*(NEn(1,2)+ICX))+Edge(j,ME2+1))/((En(e,4)-Trcr)+1)); %Estimation of En, X direction
                            En(e,3) = NBL;
                            %    subplot(1,2,1)
                            %    hold on
                            %    plot(En(e,2),En(e,1),'r*')
                            if ( ((Edgetrans(1,1)-ICY)^2 + (Edgetrans(1,2)-ICX)^2)^(0.5) < (((Edge(j,ME2)-ICY)^2+(Edge(j,ME2+1)-ICX)^2)^(0.5)) ) && (~(Edgetrans(1,1)==0 && Edgetrans(1,2)==0)) && (~(d==-1)) %Edge Remmover and modifier d is for when there is more than near edge in the boundery to find most fit
                                a = Edgetrans(1,1);
                                b = Edgetrans(1,2);
                                Edgetrans(1,1) = Edge(j,ME2);%The nearest finder
                                Edgetrans(1,2) = Edge(j,ME2+1);
                                Edge(j,ME2) = a;
                                Edge(j,ME2+1) = b;
                                side = -1;%back
                            elseif (Edgetrans(1,1)==0 && Edgetrans(1,2)==0) && (d==-1)
                                d = (abs((-(Edge(j,ME2+1)-ICX)) + m*(Edge(j,ME2)-ICY))) / sqrt(1+m^2);
                                Edgetrans(1,1) = Edge(j,ME2);%The nearest finder
                                Edgetrans(1,2) = Edge(j,ME2+1);
                                Edge(j,ME2) = 0;
                                Edge(j,ME2+1) = 0;
                                En(e,4) = En(e,4) - 1;
                                side = -1;%back
                            end
                        end
                    end
                    j = j + 1;
                end %While of Edge
                   
                
                if (side == 1) %Velocity update of lambda3 and lambda2
                    En(e,6) = abs((En(e,6)+((((((NEn(1,1)+ICY) - Edgetrans(1,1))^2 + (NEn(1,2)+ICX)-Edgetrans(1,2))^2)^(0.5))/(time_diff)))); % lambda 2
                elseif (side == -1)
                    En(e,6) = abs((En(e,6)-((((((NEn(1,1)+ICY) - Edgetrans(1,1))^2 + (NEn(1,2)+ICX)-Edgetrans(1,2))^2)^(0.5))/(time_diff)))); % lambda 3
                else % Failed matching the lambda edges
%                      En(e,:)=[];
%                     if e<2
%                      e=1;   
%                     else
%                      e=e-1;  
%                     end
                end
                %-------------Delta En L
%                 if (d == -1)
%                     % pass
%                 else
%                     if (NEn(1,1) == Inf) || (NEn(1,2) == Inf) || (En(e,1) == Inf) || (En(e,2) == Inf) || (En(e,6) == Inf)
%                         % pass
%                     else
%                         if (delta(2,1) == 0)
%                             delta(2,1) = (((NEn(1,1)+ICY-En(e,1))^2+(NEn(1,2)+ICX-En(e,2))^2)^(.5));
%                             delta(2,2) = abs(En(e,6)-Vv);
%                             delta(2,3) = abs(NBL-NBL1);
%                             delta(2,4) = delta(2,4)+1;
%                         else
%                             delta(2,4) = delta(2,4)+1;
%                             delta(2,1) = ((delta(2,4)-1)/delta(2,4))*delta(2,1)+((1/delta(2,4))*((((NEn(1,1)+ICY)-En(e,1))^2 + ((NEn(1,2)+ICX)-En(e,2))^2)^(.5)));
%                             delta(2,2) = ((delta(2,4)-1)/delta(2,4))*delta(2,1)+((1/delta(2,4))*(abs(En(e,6))-Vv));
%                             delta(2,3) = ((delta(2,4)-1)/delta(2,4))*delta(2,1)+((1/delta(2,4))*(abs(NBL-NBL1)));
%                         end
%                     end
%                 end
                %-----------------Delta En L
                %---------- REBEL EDGES!
                MAINMATCH = 0;
                j = 1;
                match = 0;
                while ( j <= (numel(Edge(:,1))) ) %Loop for failed En in boundry/ lambda_1,lambda_4 and lambda_5
                    if (((abs(((En(e,1)-Edge(j,ME2))^2) + ((En(e,2)-Edge(j,ME2+1))^2)))^(0.5)) <= NBL) && (~(Edge(j,ME2)==0 && Edge(j,ME2+1)==0)) && (((abs((-(Edge(j,ME2+1)-ICX))+m*(Edge(j,ME2)-ICY)))/sqrt(1+m^2)) < deltaT) && (d==-1)%% boundery another % Rebel classification
                        if ( En(e,4) >= Trs ) % for case d=-1 no match, Er
                            En(e,1) = NEn(1,1)+ICY;%estimation of En x
                            En(e,2) = NEn(1,2)+ICX;%estimation of En y%No En trust change
                            En(e,4) = En(e,4)-1; %disipate it by time
                            En(e,3) = NBL;
                            %      hold on
                            %      subplot(1,2,1)
                            %      plot(En(e,2)+ICX,En(e,1)+ICY,'ms')
                        elseif ( En(e,4) < Trs ) && ( En(e,4) >= Trcr)
                            En(e,4) = En(e,4)-1;
                            En(e,1) = NEn(1,1)+ICY;%estimation of En x, DO WE ADD????
                            En(e,2) = NEn(1,2)+ICX;%estimation of En y%No En trust change
                            En(e,3) = NBL;
                            %      hold on
                            %      subplot(1,2,1)
                            %      plot(En(e,2)+ICX,En(e,1)+ICY,'ms')
                            %----    L construction :D
                            el1 = 1;
                            Elkiller = 0;
                            while el1 <=(numel(alpha(:,1))) % Alpha matcher column counter
                                %------- Edge Remover from Alpha
                                if (Elkiller == 1) && ~(Elkillery1 == 0 && Elkillerx1 == 0 && Elkillery2 == 0 && Elkillerx2 == 0 && Elkillery3 == 0 && Elkillerx3 == 0)
                                    if (((alpha(el1,1) ==  Elkillery1) && (alpha(el1,2) == Elkillerx1)) || ((alpha(el1,1) ==  Elkillery2) && (alpha(el1,2) == Elkillerx2)) || ((alpha(el1,1) ==  Elkillery3) && (alpha(el1,2) == Elkillerx3)))
                                        alpha(el1,1) = 0;
                                        alpha(el1,2) = 0;
                                        alpha(el1,7) = alpha(el1,7)-1;%rebel siz
                                    end
                                    if (((alpha(el1,3) ==  Elkillery1) && (alpha(el1,4) == Elkillerx1)) || ((alpha(el1,3) ==  Elkillery2) && (alpha(el1,4) == Elkillerx2)) || ((alpha(el1,3) ==  Elkillery3) && (alpha(el1,4) == Elkillerx3)))
                                        alpha(el1,3) = 0;
                                        alpha(el1,4) = 0;
                                        alpha(el1,7) = alpha(el1,7)-1;%rebel siz
                                    end
                                    if (((alpha(el1,5) ==  Elkillery1) && (alpha(el1,6) == Elkillerx1)) || ((alpha(el1,5) ==  Elkillery2) && (alpha(el1,6) == Elkillerx2)) || ((alpha(el1,5) ==  Elkillery3) && (alpha(el1,6) == Elkillerx3)))
                                        alpha(el1,5) = 0;
                                        alpha(el1,6) = 0;
                                        alpha(el1,7) = alpha(el1,7)-1;%rebel siz
                                    end
                                end
                                %------ Equal omitter :D
                                if (alpha(el1,1) == alpha(el1,3) && alpha(el1,2)==alpha(el1,4)) && ~(alpha(el1,1)==0 && alpha(el1,2) ==0) && ~(alpha(el1,3)==0 && alpha(el1,4) ==0)
                                    alpha(el1,3) = 0;
                                    alpha(el1,4) = 0;
                                    alpha(el1,7) = alpha(el1,7)-1;
                                elseif (alpha(el1,1) == alpha(el1,5) && alpha(el1,2)==alpha(el1,6))&& ~(alpha(el1,1)==0 && alpha(el1,2) ==0) && ~(alpha(el1,5)==0 && alpha(el1,6) ==0)
                                    alpha(el1,5) = 0;
                                    alpha(el1,6) = 0;
                                    alpha(el1,7) = alpha(el1,7)-1;
                                elseif (alpha(el1,3) == alpha(el1,5) && alpha(el1,4)==alpha(el1,6)) && ~(alpha(el1,3)==0 && alpha(el1,4) ==0) && ~(alpha(el1,5)==0 && alpha(el1,6) ==0)
                                    alpha(el1,1) = alpha(el1,5);
                                    alpha(el1,2) = alpha(el1,6);
                                    alpha(el1,3) = 0;
                                    alpha(el1,4) = 0;
                                    alpha(el1,5) = 0;
                                    alpha(el1,6) = 0;
                                    alpha(el1,7) = alpha(el1,7)-1;
                                end
                                %-----
                                %------- Edge Detector in alpha
                                %if (En(e,1) == alpha(el1,el2)) && (En(e,2) == alpha(el1,el2+1)) %BL SAME AS OTHER En
                                if (Elkiller == 1) || (Elkiller == 0)
                                    if alpha(el1,7)==1
                                        MAINMATCH = MAINMATCH + 1;
                                        if (((((alpha(el1,1)-Edge(j,ME2))^2) + ((alpha(el1,2)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (alpha(el1,1)~=0 && alpha(el1,2)~=0)%Edge and alpha in bound  for L=1;  %make it in a way O is clear
                                            alpha(el1,3) = Edge(j,ME2); % KEEP O in first
                                            alpha(el1,4) = Edge(j,ME2+1);
                                            alpha(el1,7) = alpha(el1,7)+1;
                                            match = match+1;
                                        elseif (((((alpha(el1,3)-Edge(j,ME2))^2) + ((alpha(el1,4)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (alpha(el1,3)~=0 && alpha(el1,4)~=0) && alpha(el1,7)==1
                                            alpha(el1,1) = alpha(el1,3);
                                            alpha(el1,2) = alpha(el1,4);
                                            alpha(el1,3) = Edge(j,ME2);
                                            alpha(el1,4) = Edge(j,ME2+1);
                                            alpha(el1,7) = alpha(el1,7)+1;
                                            match = match+1;
                                        elseif (((((alpha(el1,5)-Edge(j,ME2))^2) + ((alpha(el1,6)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (alpha(el1,5)~=0 && alpha(el1,6)~=0) && alpha(el1,7)==1
                                            alpha(el1,1) = alpha(el1,5);
                                            alpha(el1,2) = alpha(el1,6);
                                            alpha(el1,3) = Edge(j,ME2);
                                            alpha(el1,4) = Edge(j,ME2+1);
                                            match=match+1;
                                            alpha(el1,7) = alpha(el1,7)+1;
                                        end
                                    elseif alpha(el1,7) == 2
                                        MAINMATCH=MAINMATCH+1;
                                        ml = ((Edge(j,ME2)-alpha(el1,1))/-(Edge(j,ME2+1)-alpha(el1,2)));
                                        if (((((alpha(el1,3)-Edge(j,ME2))^2) + ((alpha(el1,4)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (((abs((-(alpha(el1,4)-ICX))+ml*(alpha(el1,3)-ICY)))/sqrt(1+ml^2)) < deltaT) && ~(Elkiller==1)
                                            %MOHEM....!!!! far - near finder :\
                                            Er(numel(Er(:,1))+1,1) = Edge(j,ME2);
                                            Er(numel(Er(:,1)),2) = Edge(j,ME2+1);
                                            Er(numel(Er(:,1)),4) = Trs+2;%rebel size  L+1
                                            angle = calculate_vector_angle(Edge(j,ME2+1), Edge(j,ME2), alpha(el1,2), alpha(el1,1)); %[MODIFIED]
                                            %---------------------------- DL calculator!
                                            mk1 = (alpha(el1,3)-alpha(el1,1))/-(alpha(el1,4)-alpha(el1,2));
                                            angle0 = calculate_vector_angle(alpha(el1,4), alpha(el1,3), alpha(el1,2), alpha(el1,1)); %[MODIFIED]
                                            %angle
                                            %angle0
                                            Er(numel(Er(:,1)),3) = angle-angle0; % DL ('-' means clockwise '+' means counter-clockwise)
                                            %----------------------------------------------DL Calculator End
                                            Er(numel(Er(:,1)),5) = angle;%angle
                                            Er(numel(Er(:,1)),6) = (((Edge(j,ME2)-alpha(el1,3))^2+(Edge(j,ME2+1)-alpha(el1,4))^2)^(0.5))/(time_diff); %2 last point velocity
                                            Er(numel(Er(:,1)),7) = alpha(el1,1);%origin
                                            Er(numel(Er(:,1)),8) = alpha(el1,2);
                                            %Prepare the killer of edges!
                                            Elkiller = 1; %all match
                                            Elkillery1 = alpha(el1,1);%O of Y
                                            Elkillerx1 = alpha(el1,2);% O of X
                                            Elkillery2 = alpha(el1,3);
                                            Elkillerx2 = alpha(el1,4);
                                            Elkillery3 = Edge(j,ME2);
                                            Elkillerx3 = Edge(j,ME2+1);%Last Edge
                                            alpha(el1,:) = [];
                                            el1 = 1;
                                            match = match + 1;
                                            %          subplot(1,2,1)
                                            %          hold on
                                            %          plot([Er(:,2) Er(:,8)],[Er(:,1) Er(:,7)],'m')

                                            %maybe add another part to matrix as frame counter !! to remove
                                            %old lines that didnt match in long!
                                            %          else %good bye :D line el1
                                            %          Elkiller=2; %No match
                                            %          Elkillery1=alpha(el1,el2);
                                            %          Elkillerx1=alpha(el1,el2+1);
                                            %          alpha(el1,:)=[];
                                            %         el2=(numel(alpha(1,:)))+1; %breaker
                                            %el2=(numel(alpha(1,:)))+1; %breaker
                                        end
                                    end
                                end
                                %end BL SAME AS OTHER En
                                el1 = el1 + 1;
                            end
                            Edge(j,ME2) = 0;
                            Edge(j,ME2+1) = 0;
                        elseif (En(e,4) <= Trcr-1) %Good bye En

                            %----    L construction :D
                            el1 = 1;
                            Elkiller = 0;
                            while el1 <= (numel(alpha(:,1))) % Alpha matcher column counter
                                %------- Edge Remover from Alpha
                                if (Elkiller == 1) && ~(Elkillery1 == 0 && Elkillerx1 == 0 && Elkillery2 == 0 && Elkillerx2 == 0 && Elkillery3 == 0 && Elkillerx3 == 0)
                                    if (((alpha(el1,1) ==  Elkillery1) && (alpha(el1,2) == Elkillerx1)) || ((alpha(el1,1) ==  Elkillery2) && (alpha(el1,2) == Elkillerx2)) || ((alpha(el1,1) ==  Elkillery3) && (alpha(el1,2) == Elkillerx3)))
                                        alpha(el1,1) = 0;
                                        alpha(el1,2) = 0;
                                        alpha(el1,7) = alpha(el1,7)-1;%rebel siz
                                    end
                                    if (((alpha(el1,3) ==  Elkillery1) && (alpha(el1,4) == Elkillerx1)) || ((alpha(el1,3) ==  Elkillery2) && (alpha(el1,4) == Elkillerx2)) || ((alpha(el1,3) ==  Elkillery3) && (alpha(el1,4) == Elkillerx3)))
                                        alpha(el1,3) = 0;
                                        alpha(el1,4) = 0;
                                        alpha(el1,7) = alpha(el1,7)-1;%rebel siz
                                    end
                                    if (((alpha(el1,5) ==  Elkillery1) && (alpha(el1,6) == Elkillerx1)) || ((alpha(el1,5) ==  Elkillery2) && (alpha(el1,6) == Elkillerx2)) || ((alpha(el1,5) ==  Elkillery3) && (alpha(el1,6) == Elkillerx3)))
                                        alpha(el1,5) = 0;
                                        alpha(el1,6) = 0;
                                        alpha(el1,7) = alpha(el1,7)-1;%rebel siz
                                    end
                                end
                                %------ Equal omitter :D
                                if (alpha(el1,1) == alpha(el1,3) && alpha(el1,2)==alpha(el1,4)) && ~(alpha(el1,1)==0 && alpha(el1,2) ==0) && ~(alpha(el1,3)==0 && alpha(el1,4) ==0)
                                    alpha(el1,3) = 0;
                                    alpha(el1,4) = 0;
                                    alpha(el1,7) = alpha(el1,7)-1;
                                elseif (alpha(el1,1) == alpha(el1,5) && alpha(el1,2)==alpha(el1,6))&& ~(alpha(el1,1)==0 && alpha(el1,2) ==0) && ~(alpha(el1,5)==0 && alpha(el1,6) ==0)
                                    alpha(el1,5) = 0;
                                    alpha(el1,6) = 0;
                                    alpha(el1,7) = alpha(el1,7)-1;
                                elseif (alpha(el1,3) == alpha(el1,5) && alpha(el1,4)==alpha(el1,6)) && ~(alpha(el1,3)==0 && alpha(el1,4) ==0) && ~(alpha(el1,5)==0 && alpha(el1,6) ==0)
                                    alpha(el1,1) = alpha(el1,5);
                                    alpha(el1,2) = alpha(el1,6);
                                    alpha(el1,3) = 0;
                                    alpha(el1,4) = 0;
                                    alpha(el1,5) = 0;
                                    alpha(el1,6) = 0;
                                    alpha(el1,7) = alpha(el1,7)-1;
                                end
                                %-----
                                %------- Edge Detector in alpha
                                %if (En(e,1) == alpha(el1,el2)) && (En(e,2) == alpha(el1,el2+1)) %BL SAME AS OTHER En
                                if (Elkiller == 1) || (Elkiller == 0)
                                    if alpha(el1,7)==1
                                        MAINMATCH = MAINMATCH+1;
                                        if (((((alpha(el1,1)-Edge(j,ME2))^2) + ((alpha(el1,2)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (alpha(el1,1)~=0 && alpha(el1,2)~=0)%Edge and alpha in bound  for L=1;  %make it in a way O is clear
                                            alpha(el1,3) = Edge(j,ME2); % KEEP O in first
                                            alpha(el1,4) = Edge(j,ME2+1);
                                            alpha(el1,7) = alpha(el1,7)+1;
                                            match = match+1;
                                        elseif (((((alpha(el1,3)-Edge(j,ME2))^2) + ((alpha(el1,4)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (alpha(el1,3)~=0 && alpha(el1,4)~=0) && alpha(el1,7)==1
                                            alpha(el1,1) = alpha(el1,3);
                                            alpha(el1,2) = alpha(el1,4);
                                            alpha(el1,3) = Edge(j,ME2);
                                            alpha(el1,4) = Edge(j,ME2+1);
                                            alpha(el1,7) = alpha(el1,7)+1;
                                            match = match+1;
                                        elseif (((((alpha(el1,5)-Edge(j,ME2))^2) + ((alpha(el1,6)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (alpha(el1,5)~=0 && alpha(el1,6)~=0) && alpha(el1,7)==1
                                            alpha(el1,1) = alpha(el1,5);
                                            alpha(el1,2) = alpha(el1,6);
                                            alpha(el1,3) = Edge(j,ME2);
                                            alpha(el1,4) = Edge(j,ME2+1);
                                            match = match+1;
                                            alpha(el1,7) = alpha(el1,7)+1;
                                        end
                                    elseif alpha(el1,7)==2
                                        MAINMATCH = MAINMATCH+1;
                                        ml = ((Edge(j,ME2)-alpha(el1,1))/-(Edge(j,ME2+1)-alpha(el1,2)));
                                        if (((((alpha(el1,3)-Edge(j,ME2))^2) + ((alpha(el1,4)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && (((abs((-(alpha(el1,4)-ICX))+ml*(alpha(el1,3)-ICY)))/sqrt(1+ml^2)) < deltaT) && ~(Elkiller==1)
                                            %MOHEM....!!!! far - near finder :\
                                            Er(numel(Er(:,1))+1,1) = Edge(j,ME2);
                                            Er(numel(Er(:,1)),2) = Edge(j,ME2+1);
                                            Er(numel(Er(:,1)),4) = Trs+2;%rebel size  L+1
                                            angle = calculate_vector_angle(Edge(j,ME2+1), Edge(j,ME2), alpha(el1,2), alpha(el1,1) );
                                            %---------------------------- DL calculator!
                                            mk1 = (alpha(el1,3)-alpha(el1,1))/-(alpha(el1,4)-alpha(el1,2));
                                            angle0 = calculate_vector_angle(alpha(el1,4), alpha(el1,3), alpha(el1,2), alpha(el1,1) );
                                            %angle
                                            %angle0
                                            Er(numel(Er(:,1)),3) = angle-angle0; % DL ('-' means clockwise '+' means counter-clockwise)
                                            %----------------------------------------------DL Calculator End
                                            Er(numel(Er(:,1)),5) = angle;%angle
                                            Er(numel(Er(:,1)),6) = (((Edge(j,ME2)-alpha(el1,3))^2+(Edge(j,ME2+1)-alpha(el1,4))^2)^(0.5))/(time_diff); %2 last point velocity
                                            Er(numel(Er(:,1)),7) = alpha(el1,1);%origin
                                            Er(numel(Er(:,1)),8) = alpha(el1,2);
                                            %Prepare the killer of edges!
                                            Elkiller = 1; %all match
                                            Elkillery1 = alpha(el1,1);%O of Y
                                            Elkillerx1 = alpha(el1,2);% O of X
                                            Elkillery2 = alpha(el1,3);
                                            Elkillerx2 = alpha(el1,4);
                                            Elkillery3 = Edge(j,ME2);
                                            Elkillerx3 = Edge(j,ME2+1);%Last Edge
                                            alpha(el1,:) = [];
                                            el1 = 1;
                                            match = match + 1;
                                            %          subplot(1,2,1)
                                            %          hold on
                                            %          plot([Er(:,2) Er(:,8)],[Er(:,1) Er(:,7)],'m')

                                            %maybe add another part to matrix as frame counter !! to remove
                                            %old lines that didnt match in long!
                                            %          else %good bye :D line el1
                                            %          Elkiller=2; %No match
                                            %          Elkillery1=alpha(el1,el2);
                                            %          Elkillerx1=alpha(el1,el2+1);
                                            %          alpha(el1,:)=[];
                                            %         el2=(numel(alpha(1,:)))+1; %breaker
                                            %el2=(numel(alpha(1,:)))+1; %breaker
                                        end
                                    end
                                end

                                %end BL SAME AS OTHER En
                                el1 = el1 + 1;
                            end
                            Edge(j,ME2) = 0;
                            Edge(j,ME2+1) = 0;
                        end
                    else %En left out?
 
                     %xxxxxx   
                    end
                    j = j + 1;
                end %While of Edge
                if  match == 0 % new comer? :D G o i n t o f i r s t :D
                    j = 1;
                    while (j<=(numel(Edge(:,1)))) %Loop for failed En in boundry/ lambda_1,lambda_4 and lambda_5
                        if (((((En(e,1)-Edge(j,ME2))^2) + ((En(e,2)-Edge(j,ME2+1))^2))^(0.5)) <= NBL) && ((~(Edge(j,ME2)==0 && Edge(j,ME2+1)==0)) && (((abs((-(Edge(j,ME2+1)-ICX)) + m * (Edge(j,ME2)-ICY))) / sqrt(1+m^2)) < deltaT) && (d==-1)) %% boundery another % Rebel classification
                            alpha((numel(alpha(:,1)))+1,1) = Edge(j,ME2);
                            alpha((numel(alpha(:,1))),2) = Edge(j,ME2+1);
                            alpha((numel(alpha(:,1))),7) = 1;
                            Edge(j,ME2) = 0;
                            Edge(j,ME2+1) = 0;
                        end
                    j = j + 1;
                    end
                             %% En(e,:) = []; %left out en negative update
%                             e = e - 1;
%                             if e < 1
%                             e = 1;
%                             end
                            En(e,1) = NEn(1,1)+ICY;%estimation of En x
                            En(e,2) = NEn(1,2)+ICX;%estimation of En y%No En trust change
                            En(e,4) = En(e,4)-1; %disipate it by time
                            En(e,3) = NBL;
                end
            end % PUT HERE
            end
            e = e + 1;
        end%end of Story for En
 
    end %? HERE?
    %alpha
end

 %% rebel edge
 %-------------- Er Estimator
 %REBELIONNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNSSSSSSS! :))
if (Er == 0) % no way for youngesters :V
     % pass
else
     r = 1;
     while r <= numel(Er(:,1))
        MatchR = 0;
        if ~(Er(r,4) == -1) %not a currently calculated rebel edge L
            betar = Er(r,5);
            R = (((Er(r,1)-Er(r,7))^2)+(Er(r,2)-Er(r,8))^2)^(0.5);%The R
            x_0 = R;
            x_1 = Er(r,6);%check velocity
            [T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated E the 4 space is nutrilized to one since we want just vel
            NEr(1,1) = -(ceil(Y1(end,1))-R)*sin((pi/180)*(betar+Er(r,3)))+(Er(r,1));%estimation of En x Without removal of center ICX and ICY
            NEr(1,2) = (ceil(Y1(end,1))-R)*cos((pi/180)*(betar+Er(r,3)))+(Er(r,2));%estimation of En y
            NBL = BLS; % WILL CHANGE
            mr = (NEr(1,1)-Er(r,7))/-(Er(r,2)-Er(r,8));
            NVe = Y1(end,2);%Estimated edge velocity
            deltaTr = sqrt(deltay^2+deltaz^2);
            i = 1; %Check the Edge to find related group
            while (i<=(numel(Edge(1,:)))) %finder of lambda and Edge Match / row counter
                j = 1;
                while (j<=(numel(Edge(:,1))))%Column counter
                    %----------- Matching of REbel?
                    if (((((NEr(1,1)-(Edge(j,i)))^2) + ((NEr(1,2)-(Edge(j,i+1)))^2))^(0.5)) <= NBL) && ~(Edge(j,i)==0 && Edge(j,i+1)==0) && (((abs((-(Edge(j,i+1)))+mr*(Edge(j,i))))/sqrt(1+mr^2)) < deltaTr) && (MatchR==0)% I DO NOT PUT !!!! The zeroes in the Er omittion so! NO!!! ~(Er(r,1)==0 && Er(r,2)==0)
                        MatchR = MatchR+1;
                        me = ((Edge(j,i)-Er(r,7))/-(Edge(j,i+1)-Er(r,8)));
                        angle = calculate_vector_angle( Edge(j,i+1), Edge(j,i), Er(r,8), Er(r,7) );
                        %-------------Delta Er (NOT INCLUDED IN THIS PAPER)
                        % MAYBE REMOVE IT?
                        if (NEr(1,1)==Inf) || (NEr(1,2)==Inf) || (Er(r,1)==Inf) || (Er(r,2)==Inf) || (Er(r,6)==Inf)
                            % pass
                        else
                            if delta(3,1)==0
                                delta(3,1) = (((NEr(1,1)-Er(r,1))^2+(NEr(1,2)-Er(r,2))^2)^(.5));
                                delta(3,2) = abs(Er(r,6)-((((Edge(j,i)-Er(r,1))^2+(Edge(j,i+1)-Er(r,2))^2)^(0.5))/(time_diff)));
                                delta(3,3) = abs(Er(r,3)-((angle-(betar+Er(r,3)))));
                                delta(3,4) = delta(3,4)+1;
                            else
                                delta(3,4) = delta(3,4)+1;
                                delta(3,1) = ((delta(3,4)-1)/delta(3,4))*delta(3,1)+((1/delta(3,4))*(((NEr(1,1)-Er(r,1))^2+(NEr(1,2)-Er(r,2))^2)^(.5)));
                                delta(3,2) = ((delta(3,4)-1)/delta(3,4))*delta(3,1)+((1/delta(3,4))*(abs(Er(r,6)-((((Edge(j,i)-Er(r,1))^2+(Edge(j,i+1)-Er(r,2))^2)^(0.5))/(time_diff)))));
                                delta(3,3) = ((delta(3,4)-1)/delta(3,4))*delta(3,1)+((1/delta(3,4))*(abs(Er(r,3)-((angle-(betar+Er(r,3)))))));
                            end
                        end
                        %-----------------Delta Er
                        Er(r,6) = (((Edge(j,i)-Er(r,1))^2+(Edge(j,i+1)-Er(r,2))^2)^(0.5))/(time_diff);
                        Er(r,1) = ((((Er(r,4)-Trcr)*(NEr(1,1)))+Edge(j,i))/((Er(r,4)-Trcr)+1));
                        Er(r,2) = ((((Er(r,4)-Trcr)*(NEr(1,2)))+Edge(j,i+1))/((Er(r,4)-Trcr)+1));
                        Er(r,3) = (Er(r,3)+(angle-(betar+Er(r,3)))); %DL - Error of Edge
                        Er(r,4) = Er(r,4)+1;
                        Er(r,5) = angle;
                        hold on
                        subplot(2,2,2)
                        ploti = plot(NEr(1,2),NEr(1,1),'ms');
                        hold on
                        subplot(2,2,2)
                        ploti = plot(Edge(j,i+1),Edge(j,i),'y*');
                        hold on
                        subplot(2,2,2)
                        ploti = plot(Er(r,2),Er(r,1),'rs','LineWidth' , 2.5);
                    end
                    j = j + 1;
                end
                i = i + 2;
            end
            if MatchR==0 % No match for Er :\
                % pass
                if (Er(r,4) >= Trcr ) % update what? last point?
                    Er(r,1) = NEr(1,1);
                    Er(r,2) = NEr(1,2);%DL - Error of Edge
                    Er(r,4) = Er(r,4)-1;
                    hold on
                    subplot(2,2,2)
                    ploti = plot(Er(r,2),Er(r,1),'rs','LineWidth' , 2.5);
                    %elseif ( Er(e,4) < Trs ) && ( Er(e,4) >= Trcr)
                elseif (Er(r,4) <= Trcr-1)
                    Er(r,:) = [];
                    r = r - 1;
                end
            end
        else
            Er(r,4) = round((Trcr+Trs)/2);
        end
        r = r + 1;
    end %total Er counter 'r'
end

%-------------------------------------------------------------------------------------
%% LEFT EDGEs with En + Initiation of En
        k = (numel(En(:,1)))+1; % The size of latest En matrix (REMOVE)
        i = 1; %Check the Edge to find related group
        while (i<=(numel(Edge(1,:)))) %finder of lambda and Edge Match / row counter
            j = 1;
            while (j<=(numel(Edge(:,1))))%Column counter
                if  ~(Edge(j,i)==0 && Edge(j,i+1)==0)
                    En(k,1) = Edge(j,i); %Y
                    En(k,2) = Edge(j,i+1); %X
                    En(k,3) = BLS; %BL not good (((abs(Vv-VeS)/det(corr(EdgeNorm(:,j:j+1))))+BLS)/2)
                    En(k,4) = round((Trcr+Trs)/2);
                    %angle = calculate_vector_angle(Edge(j,i+1), Edge(j,i), ICX, ICY);
                    m = (Edge(j,i)-ICY)/-(Edge(j,i+1)-ICX);
                    angle = calculate_vector_angle( Edge(j,i+1), Edge(j,i), ICX, ICY);
                    En(k,5) = angle;
                    En(k,6) = Vv;
                    k = k + 1;
                end
                Edge(j,i) = 0;
                Edge(j,i+1) = 0;
                j = j + 1;
            end
            i = i + 2;
        end
%-----------------------------------------------------------------------------
%% En infinity ones remover
u = 1;
while u <= (numel(En(:,1)))
    if En(u,4) > Trmax
        En(u,4) = Trmax;
    end
    if En(u,4) < 0          %if trust of exsting edges goes to minus, has to be removed
        En(u,:) = [];
        u = u-1;
        if u < 1
            u = 1;
        end
    end
    if ((En(u,1) > (2*ICY)) || (En(u,2) > (2*ICX)) || (En(u,1) < 0) || (En(u,2) < 0)) %remove edges that are out of screen
        En(u,:) = [];
        u = u - 1;
    end
    u = u + 1;
end

%% match circle with rebel edges (estimating the rebel circle)
%------------Ciculing Er This must be first :))
PIN = 20; %In percentage
BetaDev = 50;
TEr = Er;
TErM = Er;
action = 0;
if (Cr==0)
    % pass
else
    ci = 1;
    while ci <= (numel(Cr(:,1))) % proper C finderd
        beta = Cr(ci,5);
        R = (((Cr(ci,1)-Cr(ci,7))^2)+(Cr(ci,2)-Cr(ci,8))^2)^(0.5);%The R
        x_0 = R;
        x_1 = Cr(ci,6);
        [T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated C the 4 space is nutrilized to one since we want just vel
        NCn(1,1) = -(ceil(Y1(end,1))-R)*sin((pi/180)*(beta))+(Cr(ci,1));%estimation of Cn x
        NCn(1,2) = (ceil(Y1(end,1))-R)*cos((pi/180)*(beta))+(Cr(ci,2));%estimation of Cn y

        u = 1;
        if Er==0
            % pass
        else
            while u<=(numel(TErM(:,1))) %Take first En
                M = 0;  %template En group
                Mk = 1;
                i = 1;
                EXTEr = TErM(u,:);
                while i<=(numel(TEr(:,1)))
                    if TEr(i,5)+TEr(i,3) > 360 || TEr(i,5)+TEr(i,3) < -360
                        AngleTEr = abs(TEr(i,5)+TEr(i,3))-360;
                    else
                        AngleTEr = TEr(i,5)+TEr(i,3);
                    end
                    if EXTEr(1,5)+EXTEr(1,3) > 360 || EXTEr(1,5)+EXTEr(1,3) < -360
                        AngleEXTEr = abs(EXTEr(1,5)+EXTEr(1,3))-360;
                    else
                        AngleEXTEr = EXTEr(1,5)+EXTEr(1,3);
                    end
                    if (AngleTEr>= AngleEXTEr-BetaDev) && ((AngleTEr <= AngleEXTEr+BetaDev)) && (abs(EXTEr(1,6))>= abs(abs(TEr(i,6))-40*abs(Vv))) && (abs(EXTEr(1,6))<= abs(abs(TEr(i,6))+40*abs(Vv))) %kick out wrong velovity jumps and angle level. may velocity be removed && (abs(TEr(i,6)) < abs(3*Vv))
                        M(Mk,1:8) = TEr(i,:);%Our Circule Mother! :D
                        Mk = Mk + 1;
                        TEr(i,:) = [];
                        i = i - 1;
                    end
                    i = i + 1;
                end
                if Mk>1
                    countin = 0;
                    R = (((M(:,1)-NCn(1,1)).^2)+((M(:,2)-NCn(1,2)).^2)).^(.5);
                    for i = 1:1:(numel(M(:,1))) %Number of in Edge
                        if R(i,1) <= Cr(ci,3)
                            countin = countin + 1;
                        end
                    end
                end
                if Mk>1
                    MY = round(mean(M(:,1))); %center of Y max(A)
                    MX = round(mean(M(:,2))); %center of X max(A)
                    DLAVE = round(mean(M(:,3)));
                    R1 = (((M(:,1)-MY).^2)+((M(:,2)-MX).^2)).^(.5);
                    MR = max(R1(:,1));
                    OY = round(mean(M(:,7)));
                    OX = round(mean(M(:,8)));

                    if ((countin/(numel(M(:,1)))) > PIN/100)
                        if mean(M(:,5)+M(:,3)) > 360
                            Mangle = abs(mean(M(:,5)+M(:,3)))-360;
                        else
                            Mangle = (mean(M(:,5)+M(:,3)));
                        end
                        if ((Mangle-(BetaDev/5) < Cr(ci,5)) && (Mangle+(BetaDev/5) > Cr(ci,5))) && (Cr(ci,6) >= mean(M(:,6))-100*abs(Vv)) && (Cr(ci,6) <= mean(M(:,6))+100*abs(Vv))%Proportion Match PIN and angular similarity and velocity alighnment &&
                            %update C!
                            action = 1;
                            Cr(ci,1) = ((((Cr(ci,4)-Trcr)*(NCn(1,1)))+MY)/((Cr(ci,4)-Trcr)+1));
                            Cr(ci,2) = ((((Cr(ci,4)-Trcr)*(NCn(1,2)))+MX)/((Cr(ci,4)-Trcr)+1)); %Estimation of En, X direction
                            Cr(ci,4) = Cr(ci,4)+1; %Trust High
                            Cr(ci,6) = (Cr(ci,6)+mean(M(:,6)))/2;
                            u = (numel(TErM(:,1)))+1;
                            %%%----- Killer of Mother Er
                            TErM = TEr;
                            %%%
                        elseif  ((Mangle-(BetaDev) < Cr(ci,5)) && (Mangle+(BetaDev) > Cr(ci,5))) && (Cr(ci,6) >= mean(M(:,6))-100*abs(Vv)) && (Cr(ci,6) <= mean(M(:,6))+100*abs(Vv)) %Somehow Match
                            action = 1;
                            Cr(ci,1) = ((((Cr(ci,4)-Trcr)*(NCn(1,1)))+MY)/((Cr(ci,4)-Trcr)+1));
                            Cr(ci,2) = ((((Cr(ci,4)-Trcr)*(NCn(1,2)))+MX)/((Cr(ci,4)-Trcr)+1)); %Estimation of En, X direction
                            Cr(ci,3) = ((((Cr(ci,4)-Trcr)*(Cr(ci,3)))+MR)/((Cr(ci,4)-Trcr)+1));
                            Cr(ci,4) = Cr(ci,4)-1; %Trust Low
                            Cr(ci,5) = ((((Cr(ci,4)-Trcr)*(Cr(ci,5)))+mean(M(:,5)+M(:,3)))/((Cr(ci,4)-Trcr)+1));
                            Cr(ci,6) = (Cr(ci,6)+mean(M(:,6)))/2;
                            Cr(ci,7) = OY;
                            Cr(ci,8) = OX;
                            u = (numel(TErM(:,1)))+1;
                            %%%----- The TEr that lost edges which grouped is placed to mother
                            %%%TErM
                            TErM = TEr;
                            %%%
                        end
                    else
                        TEr = TErM;
                    end
                end
                u = u + 1;
            end %Er counter
        end
        %--------No match to Edges? Lets weakens!
        if action == 0
            Cr(ci,1) = NCn(1,1);
            Cr(ci,2) = NCn(1,2);
            Cr(ci,4) = Cr(ci,4)-1; %Trust Low
        end
        %---------------------------
        ci = ci + 1;
    end
end
%Left Er circle new
u = 1;
if  TErM==0
    % pass
else
    while u<=(numel(TErM(:,1))) %Take first En
        M = 0;  %template En group
        Mk = 1;
        i = 1;
        EXTEr = TErM(u,:);
        while i<=(numel(TEr(:,1)))
            %        A=TEn(i,5)
            %        B=EXTEn(1,5)-BetaDev
            %        C=abs(TEn(i,6))
            %        D=abs(3*Vv)
            if TEr(i,5)+TEr(i,3) > 360 || TEr(i,5)+TEr(i,3) < -360
                AngleTEr = abs(TEr(i,5)+TEr(i,3))-360;
            else
                AngleTEr = TEr(i,5)+TEr(i,3);
            end
            if EXTEr(1,5)+EXTEr(1,3) > 360 || EXTEr(1,5)+EXTEr(1,3) < -360
                AngleEXTEr = abs(EXTEr(1,5)+EXTEr(1,3))-360;
            else
                AngleEXTEr = EXTEr(1,5)+EXTEr(1,3);
            end
            if (AngleTEr>= AngleEXTEr-BetaDev) && ((AngleTEr <= AngleEXTEr+BetaDev)) && (abs(EXTEr(1,6))>= abs(abs(TEr(i,6))-40*abs(Vv))) && (abs(EXTEr(1,6))<= abs(abs(TEr(i,6))+40*abs(Vv))) %kick out wrong velovity jumps and angle level. may velocity be removed && (abs(TEr(i,6)) < abs(3*Vv))
                M(Mk,1:8) = TEr(i,:);%Our Circule Mother! :D
                Mk = Mk + 1;
                TEr(i,:) = [];
                i = i - 1;
            end
            i = i + 1;
        end
        if Mk>1
            TErM = TEr;
            MY = round(mean(M(:,1))); %center of Y max(A)
            MX = round(mean(M(:,2))); %center of X max(A)
            DLAVE = round(mean(M(:,3)));
            R1 = (((M(:,1)-MY).^2)+((M(:,2)-MX).^2)).^(.5);
            MR = max(R1(:,1));
            OY = round(mean(M(:,7)));
            OX = round(mean(M(:,8)));
            Cr(numel(Cr(:,1))+1,1) = MY;
            Cr(numel(Cr(:,1)),2) = MX; %Estimation of En, X direction
            Cr(numel(Cr(:,1)),3) = MR;
            Cr(numel(Cr(:,1)),4) =Trs+2; %Trust Low
            angleC = 0;
            mC = (MY-OY)/-(MX-OX);
            angleC = calculate_vector_angle( MX, MY, OX, OY );
            Cr(numel(Cr(:,1)),5) = angleC+mean(M(:,3));
            Cr(numel(Cr(:,1)),6) = Vv;
            Cr(numel(Cr(:,1)),7) = OY;
            Cr(numel(Cr(:,1)),8) = OX;
            u = 1;
        else
            TEr = TErM;
        end
        u = u + 1;
    end
end



%% normal circles
%-----------Circuling En
PIN = 40; %In percentage
BetaDev = 20;
TEn = En;
TEnM = En;
if (C == 0)
    % pass
else
    ci = 1;
    while ci <= (numel(C(:,1))) % proper C finderd
        beta = C(ci,5);
        R = (((C(ci,1)-ICY)^2)+(C(ci,2)-ICX)^2)^(0.5); %The R
        x_0 = R;
        x_1 = C(ci,6);
        [T1, Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated C the 4 space is nutrilized to one since we want just vel
        NCn(1,1) = -(ceil(Y1(end,1))-R)*sin((pi/180)*(beta))+(C(ci,1)); %estimation of Cn x
        NCn(1,2) = (ceil(Y1(end,1))-R)*cos((pi/180)*(beta))+(C(ci,2)); %estimation of Cn y

        u = 1;
        if En == 0
            % pass
        else
            while u <= (numel(TEnM(:,1))) %Take first En
                M = 0;  %template En group
                Mk = 1;
                i = 1;
                action = 0;
                EXTEn = TEnM(u,:);
                while i<=(numel(TEn(:,1)))
                    %        A=TEn(i,5)
                    %        B=EXTEn(1,5)-BetaDev
                    %        C=abs(TEn(i,6))
                    %        D=abs(3*Vv)
                    if (TEn(i,5)>= EXTEn(1,5)-BetaDev) && ((TEn(i,5) <= EXTEn(1,5)+BetaDev)) && (abs(TEn(i,6)) <= abs(10*Vv)) %kick out wrong velovity jumps and angle level. may velocity be removed && (abs(TEr(i,6)) < abs(3*Vv))
                        M(Mk,1:6) = TEn(i,:);%Our Circule Main
                        Mk = Mk + 1;
                        TEn(i,:) = [];
                        i = i - 1;
                    end
                    i = i + 1;
                end
                if Mk>1
                    countin = 0;
                    R=(((M(:,1)-NCn(1,1)).^2)+((M(:,2)-NCn(1,2)).^2)).^(.5);
                    for i = 1:1:(numel(M(:,1))) %Number of in Edge
                        if R(i,1) <= C(ci,3)
                            countin = countin + 1;
                        end
                    end
                end
                if Mk>1
                    MY = round(mean(M(:,1))); %center of Y max(A)
                    MX = round(mean(M(:,2))); %center of X max(A)
                    R1 = (((M(:,1)-MY).^2)+((M(:,2)-MX).^2)).^(.5);
                    MR = max(R1(:,1));
                    if ((countin/(numel(M(:,1)))) > PIN/100)
                        if ((mean(M(:,5))-(BetaDev/5) < C(ci,5)) && (mean(M(:,5))+(BetaDev/5) > C(ci,5))) && (C(ci,6) >= mean(M(:,6))-100*abs(Vv)) &&  (C(ci,6) <= mean(M(:,6))+100*abs(Vv)) %update C!
                            action = 1;
                            C(ci,1) = ((((C(ci,4)-Trcr)*(NCn(1,1)))+MY)/((C(ci,4)-Trcr)+1));
                            C(ci,2) = ((((C(ci,4)-Trcr)*(NCn(1,2)))+MX)/((C(ci,4)-Trcr)+1)); %Estimation of En, X direction
                            C(ci,4) = C(ci,4)+1; %Trust High
                            C(ci,6) = (C(ci,6)+mean(M(:,6)))/2;
                            u = (numel(TEnM(:,1)))+1;
                            %%%----- Killer of Mother Er
                            TEnM = TEn;
                            %%%
                        elseif  ((mean(M(:,5))-(BetaDev) < C(ci,5)) && (mean(M(:,5))+(BetaDev) > C(ci,5))) && (C(ci,6) >= mean(M(:,6))-100*abs(Vv)) &&  (C(ci,6) <= mean(M(:,6))+100*abs(Vv)) %update C!%Somehow Match       action=1;
                            C(ci,1) = ((((C(ci,4)-Trcr)*(NCn(1,1)))+MY)/((C(ci,4)-Trcr)+1));
                            C(ci,2) = ((((C(ci,4)-Trcr)*(NCn(1,2)))+MX)/((C(ci,4)-Trcr)+1)); %Estimation of En, X direction
                            C(ci,3) = ((((C(ci,4)-Trcr)*(C(ci,3)))+MR)/((C(ci,4)-Trcr)+1));
                            C(ci,4) = C(ci,4)-1; %Trust Low
                            C(ci,5) = ((((C(ci,4)-Trcr)*(C(ci,5)))+mean(M(:,5)+M(:,3)))/((C(ci,4)-Trcr)+1));
                            C(ci,6) = (C(ci,6)+mean(M(:,6)))/2;
                            u = (numel(TEnM(:,1)))+1;
                            %%%----- The TEr that lost edges which grouped is placed to mother
                            %%%TErM
                            TEnM = TEn;
                            %%%
                        end
                    else
                        TEn = TEnM;
                    end
                end
                u=u+1;
            end %Er counter
        end
        %--------No match to Edges? Lets weakens!
        if action == 0
            C(ci,1) = NCn(1,1);
            C(ci,2) = NCn(1,2);
            C(ci,4) = C(ci,4)-1; %Trust Low
        end
        %---------------------------
        ci = ci+1;
    end
end
%Left Er circle new
u = 1;
if  TEnM==0
    % pass
else
    while u <= (numel(TEnM(:,1))) %Take first En
        M = 0;  %template En group
        Mk = 1;
        i = 1;
        EXTEn = TEnM(u,:);
        while i<=(numel(TEn(:,1)))
            %A=TEn(i,5)
            %B=EXTEn(1,5)-BetaDev
            %C=abs(TEn(i,6))
            %D=abs(3*Vv)
            if (TEn(i,5) >= EXTEn(1,5)-BetaDev) && ((TEn(i,5) <= EXTEn(1,5)+BetaDev)) && (abs(TEn(i,6)) < abs(3*Vv)) %kick out wrong velovity jumps and angle level. may velocity be removed && (abs(TEr(i,6)) < abs(3*Vv))
                M(Mk,1:6)=TEn(i,:);%Our Circule Mother! :D
                Mk = Mk + 1;
                TEn(i,:) = [];
                i = i - 1;
            end
            i = i + 1;
        end
        if Mk > 1
            TEnM=TEn;
            u = 1;
            MY = round(mean(M(:,1))); %center of Y max(A)
            MX = round(mean(M(:,2))); %center of X max(A)
            R1 = (((M(:,1)-MY).^2)+((M(:,2)-MX).^2)).^(.5);
            MR = max(R1(:,1));
            C(numel(C(:,1))+1,1) = MY;
            C(numel(C(:,1)),2) = MX; %Estimation of En, X direction
            C(numel(C(:,1)),3) = MR;
            C(numel(C(:,1)),4) = round((Trcr+Trs)/2); %Trust Low
            mC = (MY-ICY)/-(MX-ICX);
            angleC = calculate_vector_angle( MX, MY, ICX, ICY );
            C(numel(C(:,1)),5) = angleC;
            C(numel(C(:,1)),6) = Vv;
        else
            TEn = TEnM;
        end
        u = u + 1;
    end
end

%% updating psi and lambda
%----------Psi and lambda
u = 1;
psi = 0;

%  if C(1,4) == 0
%      C(u,:)=[];
%  end
if C==0
    C = 0;
else
    %   a=lambda;
    lambda = 0;
    while u <= (numel(C(:,1)))
        L1 = 0;
        L2 = 0;
        if C(u,4) >= Trmax
            C(u,4) = Trmax-2;
            psi(numel(psi(:,1))+1,1) = C(u,1);
            psi(numel(psi(:,1)),2) = C(u,2);
            psi(numel(psi(:,1)),3) = C(u,3);
            psi(numel(psi(:,1)),4) = 1;
        end
        if ((C(u,1) > (2*ICY)) || (C(u,2) > (2*ICX)) || (C(u,1) < 0) || (C(u,2) < 0)) %Kill more than that :D
            C(u,:) = [];
            u = u - 1;
            L1 = 1;
        end
        if (L1==1) && (u==0)
            u = u + 1;
        end
        if C(u,4) < Trcr
            C(u,:)=[];
            u = u - 1;
            L2 = 1;
        end
        if (L1==1 || L2==1) && (u==0) %either or both active and intial 0 or -1 make it 1
            u = 1;
        elseif (L1==1) && (L2==1) %both active make it one
            u = u + 1;
        end
        if C(u,4) > Trcr
            lambda(numel(lambda(:,1))+1,1) = C(u,1);
            lambda(numel(lambda(:,1)),2) = C(u,2);
            lambda(numel(lambda(:,1)),3) = C(u,3);
        end
        u = u + 1;
    end
end

r = 1;
if Cr == 0
    Cr = 0;
else
    while r <= (numel(Cr(:,1)))
        if Cr(r,4) < Trcr
            Cr(r,:) = [];
            r = r-1;
            if r < 1
                r = 1;
            end
        end
        r = r + 1;
    end
end
size(En)
end % end of function
