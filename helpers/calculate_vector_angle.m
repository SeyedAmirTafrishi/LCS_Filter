function [ angle ] = calculate_vector_angle( x1, y1, x2, y2 )
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

    angle2 = atan2(y1 - y2, x1 - x2) * (180/pi);

    angle
    angle2
end
