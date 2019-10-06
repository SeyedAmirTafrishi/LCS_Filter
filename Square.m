function [S, C, Cr] = Square(S, C, Cr, delta, Vv, Dv)
% Subsitute C and Cr to new Ct (Temprery matrix) *Done
global ICX ICY Trs
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
    Ctem_sorted = sortrows(Ctem, 7, 'descend'); %from far to near
    %%
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
            if d_m<d_tem && (( Ctem(u_mn,6) < Cmain(1,6))) %Case for exceptional circles with low velocity between two potentially matched circles
            % Failiur OBject exists with in C_A and C_B
            d_tem=d_tem-d_m; %decrease the distance
            if CanswerA ~=[]
            Ctem(numel(Ctem(:,1))+1,:) = CanswerA;  
            CanswerA=[];
            end  
            if CanswerB~= [] 
            Ctem(numel(Ctem(:,1))+1,:) = CanswerB;  
            CanswerB=[];    
            end
            FlagBReak=1;
            end
            end
                   if FlagBReak==0
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
                   else % Flag = 1 Reset!
                    u_mn =1; 
                    FlagBReak=0;
                   end
        end %NOTE: Matrix Computation is possible for future work 
        
% D_D = 0 distance of further (Maybe maximum value initiation)
% Same while loop for Ct (1) 
% if condition (furthest distance) && (Angle and velocity) Match CASE A
% WE have couple C(A) and C(B)
% Apply Unified Circle to the two circles

%% Construcst the Square from Collected Circles 
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
        S(numel(S)+1,1)=Y_o; % The location 
        S(numel(S)+1,2)=X_o;
        S(numel(S)+1,3) = a;% R of grouped Circles Y dis     
        S(numel(S)+1,5) = Trs;%Standard Trust factor for new co 
        S(numel(S)+1,4) = b;% R of grouped circles X dis  
        S(numel(S)+1,6) = calculate_vector_angle(((TempYPositive+TempYNegaitive)/2),((TempXPositive+TempXNegaitive)/2), mean(CanswerA(:,8)), mean(CanswerA(:,7)));% beta angle of square 
        S(numel(S)+1,7) = mean(CanswerA(:,7)); % Center of Frame for moving Square
        S(numel(S)+1,8) = mean(CanswerA(:,8)); % Center of Frame for moving Square  
         
    elseif  CanswerB ~= [] %Case B
        TempYPositive=max(CanswerB(:,1)+CanswerB(:,3)); % Position of summed Circles
        TempYNegaitive=min(CanswerB(:,1)-CanswerB(:,3)); %Y 
        TempXPositive=max(CanswerB(:,2)+CanswerB(:,3)); %Y 
        TempXNegaitive=min(CanswerB(:,2)-CanswerB(:,3)); %Y    
        Y_o=((TempYPositive+TempYNegaitive)/2);% The new center of construsted square by sum of circles
        X_o=((TempXPositive+TempXNegaitive)/2);
        a=abs(TempYPositive-TempYNegaitive)/2; %Y direction major
        b=abs(TempXPositive-TempXNegaitive)/2; %X direction major
        S(numel(S)+1,1)=Y_o; % The location 
        S(numel(S)+1,2)=X_o;
        S(numel(S)+1,3) = a;% R of grouped Circles Y dis     
        S(numel(S)+1,5) = Trs;%Standard Trust factor for new co 
        S(numel(S)+1,4) = b;% R of grouped circles X dis  
        S(numel(S)+1,6) = calculate_vector_angle(((TempYPositive+TempYNegaitive)/2),((TempXPositive+TempXNegaitive)/2), mean(CanswerB(:,8)), mean(CanswerB(:,7)));% beta angle of square 
        S(numel(S)+1,7) = mean(CanswerB(:,7)); % Center of Frame for moving Square
        S(numel(S)+1,8) = mean(CanswerB(:,8)); % Center of Frame for moving Square  
    else % Case Lonely
    %Lonely Square
    end



%% Does couple and circle found? Case B Certain Range Angle 
        if Canswer == [] 
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
