function [ angle ] = calculate_vector_angle( x2, y2, x1, y1 )
%calculate_vector_angle
%   Calculate the angle between two pixels

    ml = (y1-y2) / -(x1-x2);

    if (((y1 - y2)/abs((y1 - y2)))>=0 && ((x1 - x2)/abs((x1 - x2)))>=0)
        angle=(180/pi)*atan(ml);
    elseif (((y1 - y2)/abs((y1 - y2)))<0 && ((x1 - x2)/abs((x1 - x2)))>0)
        angle=(180/pi)*atan(ml);
    elseif (((y1 - y2)/abs((y1 - y2)))<0 && ((x1 - x2)/abs((x1 - x2)))<0)
        angle=(180/pi)*atan(ml)+180;
    elseif (((y1 - y2)/abs((y1 - y2)))>0 && ((x1 - x2)/abs((x1 - x2)))<0)
        angle=(180/pi)*atan(ml)+180;
    elseif (x1 - x2)==0
        angle=-((y1 - y2)/abs((y1 - y2)))*90;
    elseif (y1 - y2)==0 && ((x1 - x2)/abs((x1 - x2)))>0
        angle=0;
    elseif (y1 - y2)==0 && ((x1 - x2)/abs((x1 - x2)))<0
        angle=180;
    end

    angle2 = atan2(y2 - y1, -(x2 - x1)) * (180/pi);

    if (abs(angle - angle2) >= 10^-5)
        [angle, angle2, x1, y1, x2, y2]
    end
end
