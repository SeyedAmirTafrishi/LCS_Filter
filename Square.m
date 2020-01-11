function [S] = Square(S, C, Cr, delta, Vv, Dv, psi)
% Subsitute C and Cr to new Ct (Temprery matrix) *Done
global ICX ICY Trs time_diff Trcr Trmax
C(:,7) = ICY; % Make the normal circle in same dimension with rebel circle (Center is the center of image)
C(:,8) = ICX;
Stem = [S zeros(numel(S(:,1)),1)];
Sup=[];
countSup=1;
if Cr==0
 Ctem=C;   
else
 Ctem = [C;Cr];   
end
u_m = 1;
e_v = .4; % Deviation for circle velocity,
DeltaBeta = 16; % Deviation for beta angle
Betaconstant = 90; % Angle of two circle from each other
Betaconsame = 10;  % Angle offset for case B when there is no 90+/- angle matches of couple circles
PercntSqComp= 50; %Minimum Overlap percentage of two squares
% Take a normal circle *Done
global xd_1 yd_2 % temporory global variables in solution of tangential points on circles
%%


while u_m <= (numel(Ctem(:,1)))  %*Done Circle Counter
    %% Find a furthest distance of circle at same velocity
    u_mn = 1; %Internal counter
    ushift_mn = 0; % number of shifts in u_mn after adding back the C_A and C_B excluded by d_m
    d_tem = ICY*2; %temprorary distance, large  (Must be the maximum pixel)
    Cmain = Ctem(u_m,:); %Main Circle C_A of our Loop
    Ctem(u_m,:) = []; % Remove the chosen circle
    CanswerA = []; %Case A Matrix
    CanswerB = []; %CAse B Matrix
    CtT = []; %Stupid matrix
    FlagBReak = 0;
    %% Reorder the matrix
    % Expand another column in Ctem
    Ctem = [Ctem zeros(numel(Ctem(:,1)),1)];
    % Calculate all the distance
    u_mn = 1;
    while u_mn <= (numel(Ctem(:,1)))
        d_m = sqrt((Cmain(1,1)-Ctem(u_mn,1))^2+((Cmain(1,2)-Ctem(u_mn,2))^2));
        Ctem(u_mn,9) = d_m;
        u_mn = u_mn + 1;
    end
    
    % Now we have the circle matrix argumented with a distance (as the final column)
    % let's sort it:
    Ctem  = sortrows(Ctem, 9, 'descend'); %from far to near
    %%
    u_mn = 1;
    while  u_mn <= (numel(Ctem(:,1))) % Check all circles Case A Certain Range Angle
        %----- Angle refinment
        A = cos(Cmain(1,5)*(pi/180));
        B = sin(Cmain(1,5)*(pi/180));
        Cmain(1,5) = atan2(B,A)*(180/pi);% Note: Gives back the angle in degree
        A = cos(Ctem(u_mn,5)*(pi/180));
        B = sin(Ctem(u_mn,5)*(pi/180));
        Ctem(u_mn,5) = atan2(B,A)*(180/pi);
        %----- Angle refinment
        d_m = Ctem(u_mn,9); % Check !
       
        if ~isempty(CanswerA) || ~isempty(CanswerB) % if we have potential C_B from previous iterations (Case A and B of C_B)
            %% Calculate the whether D_m removes cooresponding C_B and add it back to loop for all previous Circles
            %The order, i change, update CanswerA/B and Ctem!!!% make it in line for easiness
            if (d_m<d_tem && (( Ctem(u_mn,6) < Cmain(1,6))))
                
                if ~isempty(CanswerA)
                    u_dm=1;
                    while u_dm<=(numel(CanswerA(:,1))) %X_B Y_B
                        x0 = [0 0 0 0];
                        xd_1=0;
                        yd_2=0;
                        if Ctem(u_mn,3)<.3
                        Ctem(u_mn,3)=2; 
                        end
                        REF = [Ctem(u_mn,2),Ctem(u_mn,1),Cmain(1,2),Cmain(1,1),Ctem(u_mn,3)];%Cmain(1,1)=C_A X_d,Y_d,X_A,Y_A,R
                        f = @(x) FindTangenfordm(x,REF); % function of dummy variable y
                        %fsolve doesnt give multiple solutons
                        F = fsolve(f,x0);
                        Point1(1,1) = real(xd_1(1,1));
                        Point2(1,1) = real(F(1,1));
                        Point1(2,1) = real(F(1,2));
                        Point2(2,1) = real(yd_2(1,1));
                        AngleCtangL=asin(sqrt((Point2(1,1)-Ctem(u_mn,1))^2+(Point1(1,1)-Ctem(u_mn,2))^2)/sqrt((Ctem(u_mn,1)-Cmain(1,1))^2+(Ctem(u_mn,2)-Cmain(1,2))^2)); %Due to circular form lower and uper has same angle but you can use only one
                        AngleCtangU=asin(sqrt((Point2(2,1)-Ctem(u_mn,1))^2+(Point1(2,1)-Ctem(u_mn,2))^2)/sqrt((Ctem(u_mn,1)-Cmain(1,1))^2+(Ctem(u_mn,2)-Cmain(1,2))^2));
                        D_Ad=sqrt((Cmain(1,2)-Ctem(u_mn,2))^2+(Cmain(1,1)-Ctem(u_mn,1))^2);
                        D_AB=sqrt((Cmain(1,2)-CanswerA(u_dm,2))^2+(Cmain(1,1)-CanswerA(u_dm,1))^2);
                        D_Bd=sqrt((CanswerA(u_dm,2)-Ctem(u_mn,2))^2+(CanswerA(u_dm,1)-Ctem(u_mn,1))^2);
                        Ahh=acos((D_AB^2+D_Ad^2-D_Bd^2)/(2*D_AB*D_Ad));
                        if (abs(AngleCtangL)>=abs(Ahh) || abs(AngleCtangU)>=abs(Ahh))  %+/- Deltae can be added due to pixal relations
                            d_tem=d_tem-d_m; %decrease the distance
                            Ctem = cat(1,CanswerA(u_dm,:),Ctem); %add CanswerA to begining of Ctem all arrays! we add it to begining because they are farthest circles and we dont need to re-do searching by u_nm of main loop
                            CanswerA(u_dm,:)=[];
                            ushift_mn=ushift_mn+1;
                            FlagBReak=1;
                        end
                        
                        u_dm=u_dm+1;
                    end
                end
                
                
                if ~isempty(CanswerB)
                    u_dm=1;
                    while u_dm<=(numel(CanswerB(:,1))) %X_B Y_B
                        x0 = [0 0 0 0];
                        xd_1=0;
                        yd_2=0;
                        if Ctem(u_mn,3)<.3 %For edge=circle with zero radius
                        Ctem(u_mn,3)=2; 
                        end
                        REF = [Ctem(u_mn,2),Ctem(u_mn,1),Cmain(1,2),Cmain(1,1),Ctem(u_mn,3)] %Cmain(1,1)=C_A X_d,Y_d,X_A,Y_A,R
                        f = @(x) FindTangenfordm(x,REF)  % function of dummy variable y
                        %fsolve doesnt give multiple solutons
                        F = fsolve(f,x0) 
                        Point1(1,1) = real(xd_1(1,1)) 
                        Point2(1,1) = real(F(1,1))  
                        Point1(2,1) = real(F(1,2)) 
                        Point2(2,1) = real(yd_2(1,1)) 
                        AngleCtangL=asin(sqrt((Point2(1,1)-Ctem(u_mn,1))^2+(Point1(1,1)-Ctem(u_mn,2))^2)/sqrt((Ctem(u_mn,1)-Cmain(1,1))^2+(Ctem(u_mn,2)-Cmain(1,2))^2)); %Due to circular form lower and uper has same angle but you can use only one
                        AngleCtangU=asin(sqrt((Point2(2,1)-Ctem(u_mn,1))^2+(Point1(2,1)-Ctem(u_mn,2))^2)/sqrt((Ctem(u_mn,1)-Cmain(1,1))^2+(Ctem(u_mn,2)-Cmain(1,2))^2));
                        D_Ad=sqrt((Cmain(1,2)-Ctem(u_mn,2))^2+(Cmain(1,1)-Ctem(u_mn,1))^2);
                        D_AB=sqrt((Cmain(1,2)-CanswerB(u_dm,2))^2+(Cmain(1,1)-CanswerB(u_dm,1))^2);
                        D_Bd=sqrt((CanswerB(u_dm,2)-Ctem(u_mn,2))^2+(CanswerB(u_dm,1)-Ctem(u_mn,1))^2);
                        Ahh=acos((D_AB^2+D_Ad^2-D_Bd^2)/(2*D_AB*D_Ad));
                        if (abs(AngleCtangL)>=abs(Ahh) || abs(AngleCtangU)>=abs(Ahh))  %+/- Deltae can be added due to pixal relations
                            d_tem=d_tem-d_m; %decrease the distance
                            Ctem = cat(1,CanswerB(u_dm,:),Ctem); %%add CanswerA to begining of Ctem all arrays!
                            CanswerB(u_dm,:)=[];
                            ushift_mn=ushift_mn+1;
                            FlagBReak=1;
                        end
                        
                        u_dm=u_dm+1;
                    end
                    
                end
                % if d_m<d_tem && (( Ctem(u_mn,6) < Cmain(1,6))) %Case for exceptional circles with low velocity between two potentially matched circles
                % % Failiur OBject exists with in C_A and C_B
                % d_tem=d_tem-d_m; %decrease the distance
                % if CanswerA ~=[]
                % % Ctem(numel(Ctem(:,1))+1,:) = CanswerA;
                % %CanswerA=[];
                % end
                % if CanswerB~= []
                % %Ctem(numel(Ctem(:,1))+1,:) = CanswerB;
                % %CanswerB=[];
                % end
            end
        end
        %%
        if FlagBReak==0    %***!!!!! ATTENTION Depending Results the sum of Case A and B is not yet decided! We may assume them together if both exist!!!**
            if ((Cmain(1,6)< Ctem(u_mn,6)+e_v) ... %CASE A
                    && (Cmain(1,6)> Ctem(u_mn,6)-e_v)) ...
                    && ((abs(Cmain(1,5)) < abs(Ctem(u_mn,5))+Betaconstant+DeltaBeta) ...
                    && (abs(Cmain(1,5)) > abs(Ctem(u_mn,5))+Betaconstant-DeltaBeta)) ...
                    && d_m<d_tem %Find the match of 90^o angle and same velocity threshold with furthest distance,NOte: it checks both +/- 90
                %---------- Remove and ADD matched circles from main matrix
                %Ctem
                if isempty(CanswerA)
                    CanswerA = Ctem(u_mn,:); %Update temporary Circle
                    d_tem = d_m;% update distance
                    Ctem(u_mn,:) = [];
                else %Canswer is not empty
                    CtT = Ctem(u_mn,:); %Goes to temporary 0
                    Ctem(u_mn,:) = [];
                    Ctem = cat(1,CanswerA,Ctem);
                    %Ctem(numel(Ctem(:,1))+1,:) = CanswerA; % PLease verify
                    CanswerA = CtT;
                    d_tem = d_m;% update distance
                end
            elseif ((Cmain(1,6) < Ctem(u_mn,6)+e_v) ... %CASE B
                    && (Cmain(1,6) > Ctem(u_mn,6)-e_v)) ...
                    && ((abs(Cmain(1,5)) < abs(Ctem(u_mn,5))+Betaconsame+DeltaBeta) ...
                    && (abs(Cmain(1,5)) > abs(Ctem(u_mn,5))+Betaconsame-DeltaBeta)) %Find same velocity threshold and small angle offset
                
                if isempty(CanswerB)
                    CanswerB = Ctem(u_mn,:); %Update temporary Circle
                    d_tem = d_m;% update distance
                    Ctem(u_mn,:) = [];
                else %Canswer is not empty
                    CtT = Ctem(u_mn,:); %Goes to temporary 0
                    Ctem(u_mn,:) = [];
                    Ctem = cat(1,CanswerB,Ctem);
                    %Ctem(numel(Ctem(:,1))+1,:) = CanswerB; % PLease verify
                    CanswerB = CtT;
                    d_tem = d_m;% update distance
                end
            else %if Case A/B fails, check whether circle can be a minor
                %% if circle u_m fails to be C_B, is it minor between C_A and C_B
                %!!!! Check if u_m be the CanswerA/B do we have to check this condition
                %again?
               % Ctem
                if ~isempty(CanswerA) % Check if u)m is a Minor circles for correspondin C_B CASE A
                    MeanYO=(mean(CanswerA(:,7))+Cmain(1,7))/2; % Cmain(u_mn) is candidate, C_A=Cmain C_B=CanswerA
                    MeanXO=(mean(CanswerA(:,8))+Cmain(1,8))/2;
                    NbetaO=calculate_vector_angle(Ctem(u_mn,2), Ctem(u_mn,1), MeanYO, MeanXO);%[MODIFIED]
                    Vmean=(mean(CanswerA(:,6))+Cmain(1,6))/2;
%                     u_m
%                     Ctem
                    SQYPositive=max([(CanswerA(:,1)+CanswerA(:,3));(Cmain(1,1)+Cmain(1,3))]); %Sqaure boundaries are determined to see whether u_mn is inside this square
                    SQYNegaitive=min([(CanswerA(:,1)-CanswerA(:,3));(Cmain(1,1)-Cmain(1,3))]); %Y
                    SQXPositive=max([(CanswerA(:,2)+CanswerA(:,3));(Cmain(1,2)+Cmain(1,3))]); %Y
                    SQXNegaitive=min([(CanswerA(:,2)-CanswerA(:,3));(Cmain(1,2)-Cmain(1,3))]);%Y
                    if Ctem(u_mn,1) > SQYNegaitive && Ctem(u_mn,1) < SQYPositive && Ctem(u_mn,2)<SQXPositive && Ctem(u_mn,2)>SQXNegaitive
                        if abs(NbetaO)+Betaconsame> abs(Cmain(1,5)) && abs(NbetaO)-Betaconsame< abs(Cmain(1,5)) && ((Vmean < Ctem(u_mn,6)+e_v) ... %CASE B
                                && (Vmean > Ctem(u_mn,6)-e_v)) ...  % add the minor circle if it follows a circular array with
                                %collected couple C_A \SumC_B
                            % there are two conditions 1: angle match 2: the circle be inside the regions of C_A and C_B
                            %IMPORTANT: CHeck after running Square whether condition cathes
                            CtT = Ctem(u_mn,:); %Goes to temporary 0
                            Ctem(u_mn,:) = [];
                            CanswerA = cat(1,CanswerA,CtT);
                            %CanswerA(numel(CanswerA(:,1))+1,:) = CtT; %works
                        end
                    end
                end
                if ~isempty(CanswerB) % The Minor circles for correspondin C_B CASE B
                    MeanYO=(mean(CanswerB(:,7))+Cmain(1,7))/2;
                    MeanXO=(mean(CanswerB(:,8))+Cmain(1,8))/2;
                    NbetaO=calculate_vector_angle(Ctem(u_mn,2), Ctem(u_mn,1), MeanYO, MeanXO);%[MODIFIED]
                    Vmean=(mean(CanswerB(:,6))+Cmain(1,6))/2;
                    SQYPositive=max([(CanswerB(:,1)+CanswerB(:,3));(Cmain(1,1)+Cmain(1,3))]) %Sqaure boundaries are determined to see whether u_mn is inside this square
                    SQYNegaitive=min([(CanswerB(:,1)-CanswerB(:,3));(Cmain(1,1)-Cmain(1,3))])  %Y
                    SQXPositive=max([(CanswerB(:,2)+CanswerB(:,3));(Cmain(1,2)+Cmain(1,3))])  %Y
                    SQXNegaitive=min([(CanswerB(:,2)-CanswerB(:,3));(Cmain(1,2)-Cmain(1,3))]) %Y
                    Ctem(u_mn,:)
%                     Cmain
%                     u_mn
%                     3
                    if Ctem(u_mn,1) > SQYNegaitive && Ctem(u_mn,1) < SQYPositive && Ctem(u_mn,2)<SQXPositive && Ctem(u_mn,2)>SQXNegaitive
                       % 1
                        if  abs(NbetaO)+Betaconsame> abs(Cmain(1,5)) && abs(NbetaO)-Betaconsame< abs(Cmain(1,5)) && ((Vmean < Ctem(u_mn,6)+e_v) ... %CASE B
                                && (Vmean > Ctem(u_mn,6)-e_v)) ...  % add the minor circle if it follows a circular array with
                            % collected couple C_A \SumC_B
                            % there are two conditions 1: angle match 2: the circle be inside the regions of C_A and C_B
                            %IMPORTANT: CHeck after running Square whether condition cathes
                           % 2
                            CtT = Ctem(u_mn,:); %Goes to temporary 0
                            Ctem(u_mn,:) = [];
                            CanswerB = cat(1,CanswerB,CtT);
                            %CanswerB(numel(CanswerB(:,1))+1,:) = CtT; %works
                        end
                    end
                end
            end
            u_mn = u_mn+1;
        else % Flag = 1 Reset!
            u_mn = u_mn + ushift_mn;
            FlagBReak = 0;
        end
    end %NOTE: Matrix Computation is possible for future work
    
    % D_D = 0 distance of further (Maybe maximum value initiation)
    % Same while loop for Ct (1)
    % if condition (furthest distance) && (Angle and velocity) Match CASE A
    % WE have couple C(A) and C(B)
    % Apply Unified Circle to the two circles
SA=[];
SB=[];
    %% Construcst the Square from Collected Circles
    if ~isempty(CanswerA) || ~isempty(CanswerB)
        if ~isempty(CanswerA) %Case A
            % Cv_Y = round(mean(CanswerA(:,1))); % Virtual Center of Y for Square
            % Cv_X = round(mean(CanswerA(:,2))); %Virtual Center of Y for Square
            TempYPositive = max(CanswerA(:,1)+CanswerA(:,3)); % Position of summed Circles
            TempYNegaitive = min(CanswerA(:,1)-CanswerA(:,3)); %Y
            TempXPositive = max(CanswerA(:,2)+CanswerA(:,3)); %Y
            TempXNegaitive = min(CanswerA(:,2)-CanswerA(:,3)); %Y
            Y_o = ((TempYPositive+TempYNegaitive)/2);% The new center of construsted square by sum of circles
            X_o = ((TempXPositive+TempXNegaitive)/2);
            a = abs(TempYPositive-TempYNegaitive)/2; %Y direction major
            b = abs(TempXPositive-TempXNegaitive)/2; %X direction major
            SA(1,1) = Y_o; % The location
            SA(1,2) = X_o;
            SA(1,3) = a;% R of grouped Circles Y dis
            SA(1,5) = Trs;%Standard Trust factor for new co
            SA(1,4) = b;% R of grouped circles X dis
            SA(1,6) = calculate_vector_angle(((TempYPositive+TempYNegaitive)/2),((TempXPositive+TempXNegaitive)/2), mean(CanswerA(:,8)), mean(CanswerA(:,7)));% beta angle of square
            SA(1,7) = mean(CanswerA(:,6));
            SA(1,8) = mean(CanswerA(:,7)); % Center of Frame for moving Square
            SA(1,9) = mean(CanswerA(:,8)); % Center of Frame for moving Square
        end
        if  ~isempty(CanswerB) %Case B
            TempYPositive = max(CanswerB(:,1)+CanswerB(:,3)); % Position of summed Circles
            TempYNegaitive = min(CanswerB(:,1)-CanswerB(:,3)); %Y
            TempXPositive = max(CanswerB(:,2)+CanswerB(:,3)); %Y
            TempXNegaitive = min(CanswerB(:,2)-CanswerB(:,3)); %Y
            Y_o = ((TempYPositive+TempYNegaitive)/2);% The new center of construsted square by sum of circles
            X_o = ((TempXPositive+TempXNegaitive)/2);
            a = abs(TempYPositive-TempYNegaitive)/2; %Y direction major
            b = abs(TempXPositive-TempXNegaitive)/2; %X direction major
            SB(1,1) = Y_o; % The location
            SB(1,2) = X_o;
            SB(1,3) = a;% R of grouped Circles Y dis
            SB(1,5) = Trs;%Standard Trust factor for new co
            SB(1,4) = b;% R of grouped circles X dis
            SB(1,6) = calculate_vector_angle(((TempYPositive+TempYNegaitive)/2),((TempXPositive+TempXNegaitive)/2), mean(CanswerB(:,8)), mean(CanswerB(:,7)));% beta angle of square
            SB(1,7) = mean(CanswerB(:,6));
            SB(1,8) = mean(CanswerB(:,7)); % Center of Frame for moving Square
            SB(1,9) = mean(CanswerB(:,8)); % Center of Frame for moving Square
        end
    else % Case Lonely
        %Lonely Square
        SB(1,1) = Cmain(1,1); % The location
        SB(1,2) = Cmain(1,2);
        SB(1,3) = Cmain(1,3);% R of grouped Circles Y dis
        SB(1,5) = Trs;%Standard Trust factor for new co
        SB(1,4) = Cmain(1,3);% R of grouped circles X dis
        SB(1,6) = Cmain(1,5);% beta angle of square
        SB(1,7)= Cmain(1,6);
        SB(1,8) = Cmain(1,7); % Center of Frame for moving Square
        SB(1,9) = Cmain(1,8); % Center of Frame for moving Square
    end
    
    %%    Create Square S(k) and estimated with any matched S'(k-1)
  if S==0
      if ~isempty(SA)
                 Sup(countSup,1) = SA(1,1);
                 Sup(countSup,2) = SA(1,2); %Estimation of Square, X direction
                 Sup(countSup,3)=  SA(1,3);
                 Sup(countSup,4)=  SA(1,4);
                 Sup(countSup,5) = SA(1,5);
                 Sup(countSup,6)= SA(1,6);
                 Sup(countSup,7)= SA(1,7);
                 Sup(countSup,8)= SA(1,8);
                 Sup(countSup,9)= SA(1,9);
                 countSup=countSup+1;
      elseif ~isempty(SB)
                 Sup(countSup,1) = SB(1,1);
                 Sup(countSup,2) = SB(1,2); %Estimation of Square, X direction
                 Sup(countSup,3)=  SB(1,3);
                 Sup(countSup,4)=  SB(1,4);
                 Sup(countSup,5) = SB(1,5);
                 Sup(countSup,6)= SB(1,6);
                 Sup(countSup,7)= SB(1,7);
                 Sup(countSup,8)= SB(1,8);
                 Sup(countSup,9)= SB(1,9);
                 countSup=countSup+1;
      end
  else
    if ~isempty(SA)
        % S(K),find S'(K) to match with our current square from real
        % data
        %1) sortrows the Complete S(.) respect to our S(K)
        %2) every step find S'(K) of each square and compare with S(K)
        %3) Move the Estimated New Square to Temp matrix and Remove it
        %from S(.)
        % FINAL step after the main while loop of Cmain, update the
        %remining sqaures in Temp. Matrix with T-1 Trust, estimate them and construct Main Square :D
        FlagSA=0;
        IK = 1; %Reordering the squares
        while IK <= (numel(Stem(:,1)))
            d_m = sqrt((SA(1,1)-Stem(IK,1))^2+((SA(1,2)-Stem(IK,2))^2));
            Stem(IK,10) = d_m;
            IK = IK + 1;
        end
        Stem = sortrows(Stem, 10, 'descend'); %from far to near
        u_sm = 1;
        while u_sm <= (numel(Stem(:,1)))
            %% Square Estimator
            %Estimate the Squares from Stem
            betaS = Stem(u_sm,6);
            R = (((Stem(u_sm,1)-Stem(u_sm,7))^2)+(Stem(u_sm,2)-Stem(u_sm,8))^2)^(0.5);%The R
            x_0 = R;
            x_1 = Stem(u_sm,6);
            options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3]);
            [T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated S the 4 space is nutrilized to one since we want just vel
            NSn(1,1) = -(Y1(end,1)-R)*sin((pi/180)*(betaS))+(Stem(u_sm,1));%estimation of Sn x
            NSn(1,2) = (Y1(end,1)-R)*cos((pi/180)*(betaS))+(Stem(u_sm,2));%estimation of Sn y
            Y_e = Stem(u_sm,1);
            X_e = Stem(u_sm,2);
            Y_o = Stem(u_sm,7);
            X_o = Stem(u_sm,8);
            a_1 = Stem(u_sm,3);
            b_1 = Stem(u_sm,4);
            Delta_r = sqrt((NSn(1,1)-Y_e)^2+(NSn(1,2)-X_e)^2);
            % Estimation of Square from Elipse_
            t_al = calculate_vector_angle( X_o, Y_o, X_e, Y_e );
            R_al = sqrt((X_o-X_e)^2+(Y_o-Y_e)^2);
            r_eE = sqrt((a_1^2*(cos(t_al*(pi/180)))^2)+(b_1^2*(sin(t_al*(pi/180)))^2));
            if a_1==0 && b_1==0 % Very small edges that transformed to circles and then to square
            a_1=2;    
            b_1=2;   
            end
            if R_al > r_eE % Case the (X_o,Y_o) is out of the ellipse
                [Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r); %A_n on X axis B_n on Y Axis
            elseif R_al <= r_eE % Case the (X_o,Y_o) is in of the ellipse
                %betaang = calculate_vector_angle(X_e,Y_e ,X_o,Y_o);    %Degree unit
                A_n = a_1+Delta_r;% Estimated Square
                B_n = b_1+Delta_r;
                Y_ef1 = NSn(1,1);
                X_ef1 = NSn(1,2);
            end
            %% Overlap computation
            [flagG, overlapPrec] = Square_Intersects(X_ef1,Y_ef1, A_n,B_n, SA(1,2),SA(1,1),SA(1,3),SA(1,4));%First is estimated square, second is 
            %  obtained square by combined circles
            if flagG==1 
            if ((SA(1,7) < Stem(u_sm,7)+e_v) ...  % Comapring the Angle and Velocity 
                    && (SA(1,7) > Stem(u_sm,7)-e_v)) ...
                    && ((abs(SA(1,6)) < abs(Stem(u_sm,6))+Betaconsame+DeltaBeta) ...
                    && (abs(SA(1,6)) > abs(Stem(u_sm,6))+Betaconsame-DeltaBeta)) ...
                    && (overlapPrec> PercntSqComp)
              
                 Sup(countSup,1) = ((((Stem(u_sm,5)-Trcr)*(Y_ef1))+SA(1,1))/((Stem(u_sm,5)-Trcr)+1));
                 Sup(countSup,2)  = ((((Stem(u_sm,5)-Trcr)*(X_ef1))+SA(1,2))/((Stem(u_sm,5)-Trcr)+1)); %Estimation of Square, X direction
                 Sup(countSup,3)= ((((Stem(u_sm,5)-Trcr)*(B_n))+SA(1,3))/((Stem(u_sm,5)-Trcr)+1));
                 Sup(countSup,4)= ((((Stem(u_sm,5)-Trcr)*(A_n))+SA(1,4))/((Stem(u_sm,5)-Trcr)+1));
                 Sup(countSup,5) = Stem(u_sm,5)+1; %Trust Low
                 Sup(countSup,6)=Stem(u_sm,6);
                 Sup(countSup,7)= ((((Stem(u_sm,5)-Trcr)*(Stem(u_sm,7)))+SA(1,7))/((Stem(u_sm,5)-Trcr)+1));% Check Velocity Formula maybe better?
                 Sup(countSup,8)=Stem(u_sm,8);
                 Sup(countSup,9)=Stem(u_sm,9);
                 Stem(u_sm,:)=[];
                 FlagSA=1;
                countSup=countSup+1;
                %Update the final Square and put it to the Sready and remove it from Stemp
                u_sm = (numel(Stem(:,1))); %break the looop! :D
                %T+1
            end
            end
            u_sm = u_sm+1;
        end
        if FlagSA == 0 % No match with the Squares of Stem Then put it as new Square 
        Sup(countSup,:)= SA; % The location
        countSup=countSup+1;
        end
    elseif ~isempty(SB)
        % Same operation for this part as SA
        FlagSB=0;
                IK = 1; %Reordering the squares
        while IK <= (numel(Stem(:,1)))
            d_m = sqrt((SB(1,1)-Stem(IK,1))^2+((SB(1,2)-Stem(IK,2))^2));
            Stem(IK,10) = d_m;
            IK = IK + 1;
        end
        Stem = sortrows(Stem, 10, 'descend'); %from far to near
        u_sm = 1;
        while u_sm <= (numel(Stem(:,1)))
            %% Square Estimator
            %Estimate the Squares from Stem
            betaS = Stem(u_sm,6);
            R = (((Stem(u_sm,1)-Stem(u_sm,7))^2)+(Stem(u_sm,2)-Stem(u_sm,8))^2)^(0.5);%The R
            x_0 = R;
            x_1 = Stem(u_sm,6);
            options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3]);
            [T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated S the 4 space is nutrilized to one since we want just vel
            NSn(1,1) = -(Y1(end,1)-R)*sin((pi/180)*(betaS))+(Stem(u_sm,1));%estimation of Sn x
            NSn(1,2) = (Y1(end,1)-R)*cos((pi/180)*(betaS))+(Stem(u_sm,2));%estimation of Sn y
            Y_e = Stem(u_sm,1);
            X_e = Stem(u_sm,2);
            Y_o = Stem(u_sm,7);
            X_o = Stem(u_sm,8);
            a_1 = Stem(u_sm,3);
            b_1 = Stem(u_sm,4);
            Delta_r = sqrt((NSn(1,1)-Y_e)^2+(NSn(1,2)-X_e)^2);
            % Estimation of Square from Elipse_
            t_al = calculate_vector_angle( X_o, Y_o, X_e, Y_e );
            R_al = sqrt((X_o-X_e)^2+(Y_o-Y_e)^2);
            r_eE = sqrt((a_1^2*(cos(t_al*(pi/180)))^2)+(b_1^2*(sin(t_al*(pi/180)))^2));
            if a_1==0 && b_1==0 % Very small edges that transformed to circles and then to square
            a_1=2;    
            b_1=2;   
            end
            if R_al > r_eE % Case the (X_o,Y_o) is out of the ellipse
                [Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r); %A_n on X axis B_n on Y Axis
            elseif R_al <= r_eE % Case the (X_o,Y_o) is in of the ellipse
                %betaang = calculate_vector_angle(X_e,Y_e ,X_o,Y_o);    %Degree unit
                A_n = a_1+Delta_r;% Estimated Square
                B_n = b_1+Delta_r;
                Y_ef1 = NSn(1,1);
                X_ef1 = NSn(1,2);
            end
            %% Overlap computation
            [flagG, overlapPrec] = Square_Intersects(X_ef1,Y_ef1, A_n,B_n, SB(1,2),SB(1,1),SB(1,3),SB(1,4));%First is estimated square, second is 
            %  obtained square by combined circles
            if flagG==1 
            if ((SB(1,7) < Stem(u_sm,7)+e_v) ...  % Comapring the Angle and Velocity 
                    && (SB(1,7) > Stem(u_sm,7)-e_v)) ...
                    && ((abs(SB(1,6)) < abs(Stem(u_sm,6))+Betaconsame+DeltaBeta) ...
                    && (abs(SB(1,6)) > abs(Stem(u_sm,6))+Betaconsame-DeltaBeta)) ...
                    && (overlapPrec> PercntSqComp)
              
                 Sup(countSup,1) = ((((Stem(u_sm,5)-Trcr)*(Y_ef1))+SB(1,1))/((Stem(u_sm,5)-Trcr)+1));
                 Sup(countSup,2)  = ((((Stem(u_sm,5)-Trcr)*(X_ef1))+SB(1,2))/((Stem(u_sm,5)-Trcr)+1)); %Estimation of Square, X direction
                 Sup(countSup,3)= ((((Stem(u_sm,5)-Trcr)*(B_n))+SB(1,3))/((Stem(u_sm,5)-Trcr)+1));
                 Sup(countSup,4)= ((((Stem(u_sm,5)-Trcr)*(A_n))+SB(1,4))/((Stem(u_sm,5)-Trcr)+1));
                 Sup(countSup,5) = Stem(u_sm,5)+1; %Trust Low
                 Sup(countSup,6)=Stem(u_sm,6);
                 Sup(countSup,7)= ((((Stem(u_sm,5)-Trcr)*(Stem(u_sm,7)))+SB(1,7))/((Stem(u_sm,5)-Trcr)+1));% Check Velocity Formula maybe better?
                 Sup(countSup,8)=Stem(u_sm,8);
                 Sup(countSup,9)=Stem(u_sm,9);
                 Stem(u_sm,:)=[];
                countSup=countSup+1;
                 FlagSB=1;
                %Update the final Square and put it to the Sready and remove it from Stemp
                u_sm = (numel(Stem(:,1))); %break the looop! :D
                %T+1
            end
            end
            u_sm = u_sm+1;
        end
        if FlagSB == 0 % No match with the Squares of Stem Then put it as new Square 
        Sup(countSup,:)= SB; % The location
        countSup=countSup+1;
        end
    end
  end
    %%
    u_m = u_m + 1;%*Done
end
if S==0
S=Sup
else
Stem(:,5)=Stem(:,5)-1;
Stem(:,10)=[];
Sup = cat(1,Stem,Sup)   
S=Sup
end
%After checking all the circles in Cmain, The remaining squares in Stemp
%are estimated with T-1 if it is less than the critical trust remove them.
%Psi and trust cleaner! :D 
u = 1;
if S==0
    S = 0;
else
    %   a=lambda;
    %lambda = 0;
    while u <= (numel(S(:,1)))
        L1 = 0;
        L2 = 0;
        if S(u,5) > Trmax
            S(u,5) = Trmax-2;
            psi(numel(psi(:,1))+1,1) = S(u,1);
            psi(numel(psi(:,1)),2) = S(u,2);
            psi(numel(psi(:,1)),3) = S(u,3);
            psi(numel(psi(:,1)),4) = 1; %--- 2 for square 1 circle, 4 to 5th (because a and b)
        end
        if ((S(u,1) > (2*ICY)) || (S(u,2) > (2*ICX)) || (S(u,1) < 0) || (S(u,2) < 0)) %Kill more than that :D
            S(u,:) = [];
            u = u - 1;
            L1 = 1;
        end
        if (L1==1) && (u==0)
            u = u + 1;
        end
        if S(u,5) < Trcr
            S(u,:)=[];
            u = u - 1;
            L2 = 1;
        end
        if (L1==1 || L2==1) && (u==0) %either or both active and intial 0 or -1 make it 1
            u = 1;
        elseif (L1==1) && (L2==1) %both active make it one
            u = u + 1;
        end
        if S(u,5) > Trcr
%            lambda(numel(lambda(:,1))+1,1) = S(u,1); % change all the psi parts to 5 array square and circle! 
%            lambda(numel(lambda(:,1)),2) = S(u,2);
%            lambda(numel(lambda(:,1)),3) = S(u,3);
        end
        u = u + 1;
    end
end


% return S delta Psi

end
