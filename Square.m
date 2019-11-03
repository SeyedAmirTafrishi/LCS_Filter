function [S, C, Cr] = Square(S, C, Cr, delta, Vv, Dv)
% Subsitute C and Cr to new Ct (Temprery matrix) *Done
global ICX ICY Trs time_diff
C(:,7) = ICX; % Make the normal circle in same dimension with rebel circle (Center is the center of image) 
C(:,8) = ICY;
Stem = [S zeros(numel(S(:,1)),1)];
Ctem = [C;Cr];
u_m = 1;
e_v = .05; % Deviation for circle velocity,
DeltaBeta = .05; % Deviation for beta angle
Betaconstant = 90; % Angle of two circle from each other
Betaconsame = 6;  % Angle offset for case B when there is no 90+/- angle matches of couple circles
% Take a normal circle *Done
global xd_1 yd_2 % temporory global variables in solution of tangential points on circles
%%
    while u_m <= (numel(Ctem(:,1)))  %*Done Circle Counter
    %% Find a furthest distance of circle at same velocity
    u_mn = 1; %Internal counter
    ushift_mn=0; % number of shifts in u_mn after adding back the C_A and C_B excluded by d_m
    d_tem = ICY*2; %temprorary distance, large  (Must be the maximum pixel)
    Cmain = Ctem(u_m,:); %Main Circle C_A of our Loop
    Ctem(u_m,:) = []; % Remove the chosen circle
    CanswerA = []; %Case A Matrix
    CanswerB = []; %CAse B Matrix
    CtT = []; %Stupid matrix
    FlagBReak=0;
    %% Reorder the matrix
    % Expand another column in Ctem
    Ctem = [Ctem zeros(numel(Ctem(:,1)),1)];
    % Calculate all the distance
    u_mn = 1;
    while u_mn <= (numel(Ctem(:,1)))
        d_m = sqrt((Cmain(1,1)-Ctem(u_mn,1))^2+((Cmain(1,2)-Ctem(u_mn,2))^2));
        Ctem(u_mn,7) = d_m;
        u_mn = u_mn + 1;
    end
    
    % Now we have the circle matrix argumented with a distance (as the final column)
    % let's sort it:
    Ctem  = sortrows(Ctem, 7, 'descend'); %from far to near
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
            d_m=Ctem(u_mn,7); % Check !
            
            if CanswerA ~= [] || CanswerB ~= [] % if we have potential C_B from previous iterations (Case A and B of C_B) 
                %% Calculate the whether D_m removes cooresponding C_B and add it back to loop for all previous Circles
                %The order, i change, update CanswerA/B and Ctem!!!% make it in line for easiness
             if (d_m<d_tem && (( Ctem(u_mn,6) < Cmain(1,6))))    
                  
                if CanswerA ~=[]
                 u_dm=1;
                    while u_dm<=(numel(CanswerA(:,1))) %X_B Y_B
                 x0 = [0 0 0 0]; 
                 xd_1=0;
                 yd_2=0;
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
                
                
                if CanswerB~= [] 
                 u_dm=1;
                    while u_dm<=(numel(CanswerB(:,1))) %X_B Y_B
                 x0 = [0 0 0 0]; 
                 xd_1=0;
                 yd_2=0;
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
%                 if d_m<d_tem && (( Ctem(u_mn,6) < Cmain(1,6))) %Case for exceptional circles with low velocity between two potentially matched circles
%                 % Failiur OBject exists with in C_A and C_B
%                 d_tem=d_tem-d_m; %decrease the distance
%                 if CanswerA ~=[]
%                % Ctem(numel(Ctem(:,1))+1,:) = CanswerA;  
%                 %CanswerA=[];
%                 end  
%                 if CanswerB~= [] 
%                 %Ctem(numel(Ctem(:,1))+1,:) = CanswerB;  
%                 %CanswerB=[];    
%                 end
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
                    if CanswerA == []
                        CanswerA = Ctem(u_mn,:); %Update temporary Circle
                        d_tem = d_m;% update distance   
                        Ctem(u_mn,:) = [];
                    else %Canswer is not empty
                        CtT = Ctem(u_mn,:); %Goes to temporary 0
                        Ctem(u_mn,:) = [];
                        Ctem(numel(Ctem(:,1))+1,:) = CanswerA; % PLease verify
                        CanswerA = CtT;
                        d_tem = d_m;% update distance  
                    end
                elseif ((Cmain(1,6) < Ctem(u_mn,6)+e_v) ... %CASE B
                    && (Cmain(1,6) > Ctem(u_mn,6)-e_v)) ... 
                    && ((abs(Cmain(1,5)) < abs(Ctem(u_mn,5))+Betaconsame+DeltaBeta) ...
                    && (abs(Cmain(1,5)) > abs(Ctem(u_mn,5))+Betaconsame-DeltaBeta)) %Find same velocity threshold and small angle offset 
                    
                    if CanswerB == []
                        CanswerB = Ctem(u_mn,:); %Update temporary Circle
                        d_tem = d_m;% update distance   
                        Ctem(u_mn,:) = [];
                    else %Canswer is not empty
                        CtT = Ctem(u_mn,:); %Goes to temporary 0
                        Ctem(u_mn,:) = [];
                        Ctem(numel(Ctem(:,1))+1,:) = CanswerB; % PLease verify
                        CanswerB = CtT;
                        d_tem = d_m;% update distance  
                    end
                else %if Case A/B fails, check whether circle can be a minor
                 %%         
 %!!!! Check if u_m be the CanswerA/B do we have to check this condition
 %again?
               if CanswerA ~= [] % Check if u)m is a Minor circles for correspondin C_B CASE A
               MeanYO=(mean(CanswerA(:,7))+Cmain(1,7))/2;
               MeanXO=(mean(CanswerA(:,8))+Cmain(1,8))/2;
               NbetaO=calculate_vector_angle(Ctem(u_mn,2), Ctem(u_mn,1), MeanYO, MeanXO);%[MODIFIED]
               SQYPositive=max([(CanswerA(:,1)+CanswerA(:,3));(Ctem(u_m,1)+Ctem(u_m,3))]); %Sqaure boundaries are determined to see whether u_mn is inside this square
               SQYNegaitive=min([(CanswerA(:,1)-CanswerA(:,3));(Ctem(u_m,1)-Ctem(u_m,3))]); %Y 
               SQXPositive=max([(CanswerA(:,2)+CanswerA(:,3));(Ctem(u_m,2)+Ctem(u_m,3))]); %Y 
               SQXNegaitive=min([(CanswerA(:,2)-CanswerA(:,3));(Ctem(u_m,2)-Ctem(u_m,3))]);%Y  
                if Ctem(u_m,1) > SQYNegaitive && Ctem(u_m,1) < SQYPositive && Ctem(u_m,2)<SQXPositive && Ctem(u_m,2)>SQXNegaitive
                    if abs(NbetaO)+Betaconsame> abs(Ctem(u_mn,5)) && abs(NbetaO)-Betaconsame< abs(Ctem(u_mn,5)) && ((Cmain(1,6) < Ctem(u_mn,6)+e_v) ... %CASE B
                    && (Cmain(1,6) > Ctem(u_mn,6)-e_v)) ...  % add the minor circle if it follows a circular array with 
                %collected couple C_A \SumC_B
                % there are two conditions 1: angle match 2: the circle be inside the regions of C_A and C_B    
                %IMPORTANT: CHeck after running Square whether condition cathes
                        CtT = Ctem(u_mn,:); %Goes to temporary 0
                        Ctem(u_mn,:) = [];
                        CanswerA(numel(CanswerA(:,1))+1,:) = CtT; %works
                    end
                end
               end
               if CanswerB ~= [] % The Minor circles for correspondin C_B CASE B
               MeanYO=(mean(CanswerA(:,7))+Cmain(1,7))/2;
               MeanXO=(mean(CanswerA(:,8))+Cmain(1,8))/2;
               NbetaO=calculate_vector_angle(Ctem(u_mn,2), Ctem(u_mn,1), MeanYO, MeanXO);%[MODIFIED]
               SQYPositive=max([(CanswerB(:,1)+CanswerB(:,3));(Ctem(u_m,1)+Ctem(u_m,3))]); %Sqaure boundaries are determined to see whether u_mn is inside this square
               SQYNegaitive=min([(CanswerB(:,1)-CanswerB(:,3));(Ctem(u_m,1)-Ctem(u_m,3))]); %Y 
               SQXPositive=max([(CanswerB(:,2)+CanswerB(:,3));(Ctem(u_m,2)+Ctem(u_m,3))]); %Y 
               SQXNegaitive=min([(CanswerB(:,2)-CanswerB(:,3));(Ctem(u_m,2)-Ctem(u_m,3))]);%Y  
                if Ctem(u_m,1) > SQYNegaitive && Ctem(u_m,1) < SQYPositive && Ctem(u_m,2)<SQXPositive && Ctem(u_m,2)>SQXNegaitive
                    if abs(NbetaO)+Betaconsame> abs(Ctem(u_mn,5)) && abs(NbetaO)-Betaconsame< abs(Ctem(u_mn,5)) && ((Cmain(1,6) < Ctem(u_mn,6)+e_v) ... %CASE B
                    && (Cmain(1,6) > Ctem(u_mn,6)-e_v)) ...  % add the minor circle if it follows a circular array with 
                %collected couple C_A \SumC_B
                % there are two conditions 1: angle match 2: the circle be inside the regions of C_A and C_B    
                %IMPORTANT: CHeck after running Square whether condition cathes
                        CtT = Ctem(u_mn,:); %Goes to temporary 0
                        Ctem(u_mn,:) = [];
                        CanswerB(numel(CanswerB(:,1))+1,:) = CtT; %works
                    end
                end    
               end  
            
                    
                end

              
        u_mn = u_mn+1;    
        else % Flag = 1 Reset!
        u_mn =u_mn+ushift_mn; 
        FlagBReak=0;
       end           
        end %NOTE: Matrix Computation is possible for future work 
   
% D_D = 0 distance of further (Maybe maximum value initiation)
% Same while loop for Ct (1) 
% if condition (furthest distance) && (Angle and velocity) Match CASE A
% WE have couple C(A) and C(B)
% Apply Unified Circle to the two circles

%% Construcst the Square from Collected Circles 
   if CanswerA ~= [] || CanswerB ~= []
    if CanswerA ~= [] %Case A 
       % Cv_Y = round(mean(CanswerA(:,1))); % Virtual Center of Y for Square 
       % Cv_X = round(mean(CanswerA(:,2))); %Virtual Center of Y for Square     
        TempYPositive=max(CanswerA(:,1)+CanswerA(:,3)); % Position of summed Circles
        TempYNegaitive=min(CanswerA(:,1)-CanswerA(:,3)); %Y 
        TempXPositive=max(CanswerA(:,2)+CanswerA(:,3)); %Y 
        TempXNegaitive=min(CanswerA(:,2)-CanswerA(:,3)); %Y    
        Y_o=((TempYPositive+TempYNegaitive)/2);% The new center of construsted square by sum of circles
        X_o=((TempXPositive+TempXNegaitive)/2);
        a=abs(TempYPositive-TempYNegaitive)/2; %Y direction major
        b=abs(TempXPositive-TempXNegaitive)/2; %X direction major
        SA(1,1)=Y_o; % The location 
        SA(1,2)=X_o;
        SA(1,3) = a;% R of grouped Circles Y dis     
        SA(1,5) = Trs;%Standard Trust factor for new co 
        SA(1,4) = b;% R of grouped circles X dis  
        SA(1,6) = calculate_vector_angle(((TempYPositive+TempYNegaitive)/2),((TempXPositive+TempXNegaitive)/2), mean(CanswerA(:,8)), mean(CanswerA(:,7)));% beta angle of square 
        SA(1,7) = mean(CanswerA(:,7)); % Center of Frame for moving Square
        SA(1,8) = mean(CanswerA(:,8)); % Center of Frame for moving Square  
    end   
     if  CanswerB ~= [] %Case B
        TempYPositive=max(CanswerB(:,1)+CanswerB(:,3)); % Position of summed Circles
        TempYNegaitive=min(CanswerB(:,1)-CanswerB(:,3)); %Y 
        TempXPositive=max(CanswerB(:,2)+CanswerB(:,3)); %Y 
        TempXNegaitive=min(CanswerB(:,2)-CanswerB(:,3)); %Y    
        Y_o=((TempYPositive+TempYNegaitive)/2);% The new center of construsted square by sum of circles
        X_o=((TempXPositive+TempXNegaitive)/2);
        a=abs(TempYPositive-TempYNegaitive)/2; %Y direction major
        b=abs(TempXPositive-TempXNegaitive)/2; %X direction major
        SB(1,1)=Y_o; % The location 
        SB(1,2)=X_o;
        SB(1,3) = a;% R of grouped Circles Y dis     
        SB(1,5) = Trs;%Standard Trust factor for new co 
        SB(1,4) = b;% R of grouped circles X dis  
        SB(1,6) = calculate_vector_angle(((TempYPositive+TempYNegaitive)/2),((TempXPositive+TempXNegaitive)/2), mean(CanswerB(:,8)), mean(CanswerB(:,7)));% beta angle of square 
        SB(1,7) = mean(CanswerB(:,7)); % Center of Frame for moving Square
        SB(1,8) = mean(CanswerB(:,8)); % Center of Frame for moving Square  
     end
    else % Case Lonely
        %Lonely Square
        SB(1,1)= Cmain(1,1); % The location 
        SB(1,2)= Cmain(1,2);
        SB(1,3) = Cmain(1,3);% R of grouped Circles Y dis     
        SB(1,5) = Trs;%Standard Trust factor for new co 
        SB(1,4) = Cmain(1,4);% R of grouped circles X dis  
        SB(1,6) = Cmain(1,5);% beta angle of square 
        SB(1,7) = Cmain(1,7); % Center of Frame for moving Square
        SB(1,8) = Cmain(1,8); % Center of Frame for moving Square 
    end

%%    Create Square S(k) and estimated with any matched S'(k-1)   
 
        if SA ~= [] 
            % S(K),find S'(K) to match with our current square from real
            % data
            %1) sortrows the Complete S(.) respect to our S(K)  
            %2) every step find S'(K) of each square and compare with S(K)
            %3) Move the Estimated New Square to Temp matrix and Remove it
            %from S(.) 
            % FINAL step after the main while loop of Cmain, update the
            %remining sqaures in Temp. Matrix with T-1 Trust, estimate them and construct Main Square :D  
             IK = 1; %Reordering the squares
         while IK <= (numel(Stem(:,1)))
        d_m = sqrt((SA(1,1)-Stem(IK,1))^2+((SA(1,2)-Stem(IK,2))^2));
        Stem(IK,7) = d_m;
        IK = IK + 1;
         end
        Stem  = sortrows(Stem, 7, 'descend'); %from far to near 
        u_sm=1;
         while u_sm <= (numel(Stem(:,1)))
             %% Square Estimator
           %Estimate the Squares from Stem  
        betaS = Stem(u_sm,6);
        R = (((Stem(u_sm,1)-Stem(u_sm,7))^2)+(Stem(u_sm,2)-Stem(u_sm,8))^2)^(0.5);%The R
        x_0 = R;
        x_1 = Stem(u_sm,6);
        [T1,Y1] = ode45(@EdgeTR,[0 time_diff],[x_0 x_1],options); %location of estimated S the 4 space is nutrilized to one since we want just vel
        NSn(1,1) = -(Y1(end,1)-R)*sin((pi/180)*(betaS))+(Stem(u_sm,1));%estimation of Sn x
        NSn(1,2) = (Y1(end,1)-R)*cos((pi/180)*(betaS))+(Stem(u_sm,2));%estimation of Sn y   
        Y_e=Stem(u_sm,1);
        X_e=Stem(u_sm,2); 
        Y_o=Stem(u_sm,7);
        X_o=Stem(u_sm,8);
        a_1=Stem(u_sm,3);
        b_1=Stem(u_sm,4);
        Delta_r=sqrt((NSn(1,1)-Y_e)^2+(NSn(1,2)-X_e)^2);
        % Estimation of Square from Elipse_
            t_al = calculate_vector_angle( X_o, Y_o, X_e, Y_e );
            R_al=sqrt((X_o-X_e)^2+(Y_o-Y_e)^2);
            r_eE= sqrt((a_1^2*(cos(t_al*(pi/180)))^2)+(b_1^2*(sin(t_al*(pi/180)))^2));
            if R_al > r_eE % Case the (X_o,Y_o) is out of the ellipse
            [Y_ef1,X_ef1,A_n,B_n] = EstimSquare(X_e,Y_e,X_o,Y_o,a_1,b_1,Delta_r); %A_n on X axis B_n on Y Axis
            elseif R_al <= r_eE % Case the (X_o,Y_o) is in of the ellipse
            %betaang = calculate_vector_angle(X_e,Y_e ,X_o,Y_o);    %Degree unit
            A_n=a_1+Delta_r;
            B_n=b_1+Delta_r;
            Y_ef1=NSn(1,1);
            X_ef1=NSn(1,2);
            end
        
            %%
            if (1)%1.Beta_angle respect to approx. origin of SA 2.check velocity 3. The Percentage of involvement if 60% of estimated square is in SA/SB we are done
          %Update the final Square and put it to the Sready and remove it from Stemp
          u_sm=(numel(Stem(:,1))); %break the looop! :D 
          %T+1
            end
          u_sm=u_sm+1;   
         end  
        elseif  SB ~= []  
         % Same operation for this part as SA
    
        end
%%
u_m = u_m + 1;%*Done
    end
 %After checking all the circles in Cmain, The remaining squares in Stemp
 %are estimated with T-1 if it is less than the critical trust remove them.
    
% return S delta Psi 

 end
