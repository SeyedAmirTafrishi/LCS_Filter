%Problem Solved :D center of Square

CanswerA = [1 2 5; 3 10 4; 12 3 5];
for i=1:1:numel(CanswerA(:,1))
    xunit = 0;
    yunit = 0;
    th = 0:pi/50:2*pi; %for loop for creating circle
    CB = 1;
    xunit = (CanswerA(i,3) ) * cos(th) + CanswerA(i,2); %equation of circle :D
    yunit = (CanswerA(i,3) ) * sin(th) + CanswerA(i,1);
    ploti = plot(xunit, yunit,'g');% Ellipse
    hold on    
end

Cv_Y=(mean(CanswerA(:,1))) %Virtual Center of Y for Square 
Cv_X=(mean(CanswerA(:,2))) %Virtual Center of Y for Square 
plot(Cv_X,Cv_Y,'- xb','MarkerSize', 18,'LineWidth' , 2.5)  

% CanswerAN=CanswerA(:,1)-Cv_Y
% CanswerAN=CanswerA(:,2)-Cv_X

TempYPositive=max(CanswerA(:,1)+CanswerA(:,3)) %Y It find the max of array!!!!
TempYNegaitive=min(CanswerA(:,1)-CanswerA(:,3)) %Y 
TempXPositive=max(CanswerA(:,2)+CanswerA(:,3)) %Y 
TempXNegaitive=min(CanswerA(:,2)-CanswerA(:,3)) %Y    

Y_o=((TempYPositive-TempYNegaitive)/2)
X_o=((TempXPositive-TempXNegaitive)/2)
a=abs(TempYPositive-TempYNegaitive)/2 %height
b=abs(TempXPositive-TempXNegaitive)/2
plot(X_o,Y_o,'- *b','MarkerSize', 18,'LineWidth' , 2.5)  

plot([TempXPositive TempXPositive],[TempYNegaitive TempYPositive],'r')
plot([TempXNegaitive TempXPositive],[TempYPositive TempYPositive],'r')
plot([TempXNegaitive TempXNegaitive],[TempYNegaitive TempYPositive],'r')
plot([TempXNegaitive TempXPositive],[TempYNegaitive TempYNegaitive],'r')
