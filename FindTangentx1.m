function F = FindTangentx1(x, REF)
X_e = REF(1,1); %Elipse center coordinate
Y_e = REF(1,2);
X_o = REF(1,3);
Y_o = REF(1,4); %Original coordinate
a_1 = REF(1,5); %elipse minor past
b_1 = REF(1,6); %elipse major past
global x_1 y_2 y_1 x_2
% x_2=Y_e-((b_1/a_1)*sqrt((a_1)^2-(x(1)-X_e)^2))
% x_1=X_e-((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2))
% F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
% F(2) = ((x(2)-Y_e)^2/(b_1^2))+((x(1)-X_e)^2/(a_1^2))-1;

% Solution 1 
% x_1=X_e+((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2))
% F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));

% Solution 2

% if (x(1)-Y_e)^2>(b_1)^2 && (x(2)-X_e)^2>(a_1)^2
%  y_1 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(1)-X_e)^2));    
%  F(1) = (a_1^2*(Y_o-y_1)*(y_1-Y_e)) + (b_1^2*(X_o-x(1))*(x(1)-X_e));   
%  
%  
% elseif (x(1)-Y_e)^2>(b_1)^2 && (x(1)-Y_e)^2<(b_1)^2
%     
% elseif (x(1)-Y_e)^2<(b_1)^2 && (x(2)-X_e)^2>(a_1)^2
%     
% else
%     
% end

%%
 % ADD SAME COMPARISON
 
   
    x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
    F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
    y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  
    F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
    if isreal(x_1)==1 && isreal(y_2)==1
    elseif isreal(x_1)==1 && isreal(y_2)~=1
        if isreal(x_1)==1 
        else
        x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
        F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e)); 
        end
        if isreal(y_2)==1 
        else
         y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  
         F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
        end    
        if isreal(x_1)==1 
        else
        x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
        F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e)); 
        end
        if isreal(y_2)==1 
        else
        y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  
        F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
        end
    elseif isreal(x_1)~=1 && isreal(y_2)~=1
        if isreal(x_1)==1 
        else
        x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
        F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e)); 
        end
        if isreal(y_2)==1 
        else
         y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  
         F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
        end    
        if isreal(x_1)==1 
        else
        x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
        F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e)); 
        end
        if isreal(y_2)==1 
        else
        y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  
        F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
        end  
    else

    x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
    F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
    y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  
    F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));    
    end
%% SAVEd




% if (X_o>abs(X_e-a_1) && Y_o>abs(Y_e-b_1))%+ + 
%     1
%    if Y_o<Y_e
%    if abs(Y_o-Y_e)<b_1
%        11
%      x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));     
%    else
%        12
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
%    end
%    else
%     if abs(Y_o-Y_e)<b_1
%        13
%      x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));     
%    else
%        14
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2));
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
%    end     
%    end
% elseif (X_o<abs(X_e-a_1) && Y_o>abs(Y_e-b_1))
%     2
%    if Y_o<Y_e
%    if abs(Y_o-Y_e)<b_1
%        21
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));     
%    else
%        22
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
%    end
%    else
%     if abs(Y_o-Y_e)<b_1
%        23
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));     
%    else
%        24
%      x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
%    end     
%    end
%     
% elseif (X_o<X_e-a_1 && Y_o<Y_e-b_1)
%     3
%     x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2));
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e))+(b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
% 
% else%- +
%    4
%    if X_o<X_e
%    if abs(X_o-X_e)<a_1
%        41
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));     
%    else
%        42
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
%    end
%    else
%     if abs(X_o-X_e)<a_1
%        43
%      x_1 = X_e + ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e - ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));     
%    else
%        44
%      x_1 = X_e - ((a_1/b_1)*sqrt((b_1)^2-(x(1)-Y_e)^2)); 
%     % save x_1
%     F(1) = (a_1^2*(Y_o-x(1))*(x(1)-Y_e)) + (b_1^2*(X_o-x_1)*(x_1-X_e));
%     y_2 = Y_e + ((b_1/a_1)*sqrt((a_1)^2-(x(2)-X_e)^2));  %change this guy
%     % save y_2
%     F(2) = (a_1^2*(Y_o-y_2)*(y_2-Y_e)) + (b_1^2*(X_o-x(2))*(x(2)-X_e));  
%    end     
%    end
%    
%     
%     %elseif if X_o=X_e    
% end




 
end
