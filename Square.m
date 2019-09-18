function [S, C, Cr] = Square(S, C, Cr, delta, Vv, Dv)
% Subsitute C and Cr to new Ct (Temprery matrix) *Done
global ICX ICY
C(:,7) = ICX; % Make the normal circle in same dimension with rebel circle (Center is the center of image) 
C(:,8) = ICY;
Ctem = [C;Cr];
u_m = 1;
e_v = .05; % Deviation for circle velocity,
DeltaBeta = .05; % Deviation for beta angle
Betaconstant = 90; % Angle of two circle from each other
Betaconsame = 10;  % Angle offset for case B when there is no 90+/- angle matches of couple circles
% Take a normal circle *Done

%%
while u_m <= (numel(Ctem(:,1)))  %*Done Circle Counter
    %% Find a furthest distance of circle at same velocity
    u_mn = 1; %Internal counter
    d_tem = 10000; %temprorary distance, large 
    Cmain = Ctem(u_m,:); %Main Circle
    Ctem(u_m,:) = []; % Remove the chosen circle
    CanswerA = []; %Case A Matrix
    CanswerB = []; %CAse B Matrix
    CtT = []; %Stupid matrix
        while  u_mn <= (numel(Ctem(:,1))) % Check all circles Case A Certain Range Angle 
            %----- Angle refinment 
            A = cos(Cmain(1,5)*(pi/180)); 
            B = sin(Cmain(1,5)*(pi/180));
            Cmain(1,5) = atan2(B,A)*(180/pi);% Note: Gives back the angle in degree
            A = cos(Ctem(u_mn,5)*(pi/180));
            B = sin(Ctem(u_mn,5)*(pi/180));
            Ctem(u_mn,5) = atan2(B,A)*(180/pi);
            %----- Angle refinment
            d_m=sqrt((Cmain(1,1)-Ctem(u_mn,1))^2+((Cmain(1,2)-Ctem(u_mn,2))^2)); % Check !
                if ((Cmain(1,6)< Ctem(u_mn,6)+e_v) ...
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
                elseif ((Cmain(1,6) < Ctem(u_mn,6)+e_v) ...
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
                end
        u_mn = u_mn+1;
        end %NOTE: Matrix Computation is possible for future work 
        
% D_D = 0 distance of further (Maybe maximum value initiation)
% Same while loop for Ct (1) 
% if condition (furthest distance) && (Angle and velocity) Match CASE A
% WE have couple C(A) and C(B)
% Apply Unified Circle to the two circles

%% Construcst the Square from Collected Circles 
    if CanswerA ~= [] %Case A 
        Cv_Y = round(mean(CanswerA(:,1))); % Virtual Center of Y for Square 
        Cv_X = round(mean(CanswerA(:,2))); %Virtual Center of Y for Square     
        TempD = max(CanswerA(:,1)-Cv_Y+CanswerA(:,3)) %Y 
        
        
        S(numel(S)+1,3) = % R of grouped Circles Y dis     
        S(numel(S)+1,4) = % R of grouped circles X dis    
    elseif  CanswerB ~= [] %Case B
        % need coding
    
    else % Case Lonely
    
    end



%% Does couple and circle found? Case B Certain Range Angle 
        if Canswer == [] 
%             MS=[]; %Collector temporary matrix
%             nMS=0; % Counter of collector temp. matrix
%                 while  u_mn <= (numel(Ctem(:,1))) % Check all circles 
%             %----- Angle refinment 
%             A=cos(Cmain(1,5)*(pi/180)); 
%             B=sin(Cmain(1,5)*(pi/180));
%             Cmain(1,5)=atan2(B,A)*(180/pi);% Note: Gives back the angle in degree
%             A=cos(Ctem(u_mn,5)*(pi/180));
%             B=sin(Ctem(u_mn,5)*(pi/180));
%             Ctem(u_mn,5)=atan2(B,A)*(180/pi);
%             %----- Angle refinment
%                 if ((Cmain(1,6)< Ctem(u_mn,6)+e_v) && (Cmain(1,6)> Ctem(u_mn,6)-e_v)) && ((abs(Cmain(1,5))< abs(Ctem(u_mn,5))+Betaconsame+DeltaBeta) && (abs(Cmain(1,5))> abs(Ctem(u_mn,5))+Betaconsame-DeltaBeta)) %Find same velocity threshold and small angle offset 
%                     nMS=nMS+1;
%                     MS(nMS,:)=Ctem(u_mn,:); %Update temporary Circle
%                 end
%                     u_mn=u_mn+1;
%                 end %NOTE: Matrix Computation is possible for future work 
%%    Create Square S(k) and estimated with any matched S'(k-1)        
            if MS == [] %Lonely Circle :D
                %------- Estimation and Matching Squares 
                
            else
            
            
            end
                      
        else % Remove
    
    
        end



% State: Search loop for verifying any circle exsitis that makes two C(A) and C(B)
% as seperate circle 
% IF condition (2) V*_C < V_C(B), V_C(A)
% Same while loop excluding C(A) and C(B) 
% Find matches of exstiing circles
% if matches are over
% State Remove the circles from Ct
% Square While loop
% Find Square matching the create square from C(A) and C(B) couples
% Estimate the Square 
% end all loops 
% Else of If CASE A (No match) (1) CASE B
% Find C(A) where it is align with other circles in delta Beta
% Loop for C_t 
% Construct the unified circle
% Loop for Square search
% Estimate for Square
% Remove the Circles From C_t
%ELSE of IF (2) 
%decrease the D_D
% Return to while loop of 1 to fin

% end

% 
    % Does couple C(A) and C(B) found?

    % Any normal/rebel circle found?

    % Include any circle that matches same velocity and certain angle

    % Create square and estimated with any matched S'(k-1)

    % Does just unmatched circles left (no match of case A and B)?
u_m = u_m + 1;%*Done
end


% Match remaining circles as a squre
% S = []

% return S delta Psi 

end
