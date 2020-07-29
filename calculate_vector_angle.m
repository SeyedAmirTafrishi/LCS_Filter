function [ angle ] = calculate_vector_angle( x1, y1, ICX, ICY )
%calculate_vector_angle
%   Calculate the angle between two pixels
%   [ICX, ICY] is the origin
%   [x1, y1] is the new location

    ml = (y1-ICY) / -(x1-ICX);

    if (((y1 - ICY)/abs((y1 - ICY)))>=0 && ((x1 - ICX)/abs((x1 - ICX)))>=0)
        angle=(180/pi)*atan(ml);
    elseif (((y1 - ICY)/abs((y1 - ICY)))<0 && ((x1 - ICX)/abs((x1 - ICX)))>0)
        angle=(180/pi)*atan(ml);
    elseif (((y1 - ICY)/abs((y1 - ICY)))<0 && ((x1 - ICX)/abs((x1 - ICX)))<0)
        angle=(180/pi)*atan(ml)+180;
    elseif (((y1 - ICY)/abs((y1 - ICY)))>0 && ((x1 - ICX)/abs((x1 - ICX)))<0)
        angle=(180/pi)*atan(ml)+180;
    elseif (x1 - ICX)==0
        angle=-((y1 - ICY)/abs((y1 - ICY)))*90;
    elseif (y1 - ICY)==0 && ((x1 - ICX)/abs((x1 - ICX)))>0
        angle=0;
    elseif (y1 - ICY)==0 && ((x1 - ICX)/abs((x1 - ICX)))<0
        angle=180;
    else
     angle=  0;   
    end
angle=-angle;
  %  angle2 = -atan2(y1 - ICY, (x1 - ICX)) * (180/pi);
%   if isempty(angle)
%     angle=  0;
%   else
%    angle=-angle;   
 % end

    %[angle, angle2, x1, y1, xo, yo];
end
