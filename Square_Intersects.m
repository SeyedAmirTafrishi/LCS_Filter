function [flag, overlapArea] = Square_Intersects(Lx_1,Ly_1, a_1,b_1, Lx_2,Ly_2,a_2,b_2) 
% Square_Intersects()
% Test case: [flag, overlapArea] = Square_Intersects(60, 15, 5, 10, 50, 30, 11, 14)

%Find where Sq 1 is located around 2
if (Lx_1+a_1<Lx_2-a_2 || Lx_1-a_1>Lx_2+a_2 ) ||  (Ly_1-b_1>Ly_2+b_2 || Ly_1+b_1<Ly_2-b_2 )% Sq_1 is out of Sq_2
    flag = 0; % If out then good bye
    overlapArea = 0;
else %Always in
    flag = 1; % if in who is bigger?
	% Find intersecting area!
	% need to know direction

    %The overlap area can be computed as follows:
    Left1 = Lx_1-a_1;
    Top1 = Ly_1+b_1;
    Bottom1 = Ly_1-b_1;
    Right1 = Lx_1+a_1;

    Left2 = Lx_2-a_2;
    Top2 = Ly_2+b_2;
    Bottom2 = Ly_2-b_2;
    Right2 = Lx_2+a_2;

    x_overlap = abs(min([Right1, Right2]) - max([Left1, Left2]));
    y_overlap =abs(max([Bottom1, Bottom2]) - min([Top1, Top2]));
    overlapArea = x_overlap * y_overlap;
end

% hold on
% plot(Lx_1,Ly_1,'- *r','MarkerSize', 18,'LineWidth' , 2.5)
% hold on
% plot([Lx_1+a_1 Lx_1+a_1],[Ly_1-b_1 Ly_1+b_1],'b')
% plot([Lx_1-a_1 Lx_1+a_1],[Ly_1+b_1 Ly_1+b_1],'b')
% plot([Lx_1-a_1 Lx_1-a_1],[Ly_1-b_1 Ly_1+b_1],'b')
% plot([Lx_1-a_1 Lx_1+a_1],[Ly_1-b_1 Ly_1-b_1],'b')
% 
% 
% hold on
% plot(Lx_2,Ly_2,'- *g','MarkerSize', 18,'LineWidth' , 2.5)
% hold on
% plot([Lx_2+a_2 Lx_2+a_2],[Ly_2-b_2 Ly_2+b_2],'m')
% plot([Lx_2-a_2 Lx_2+a_2],[Ly_2+b_2 Ly_2+b_2],'m')
% plot([Lx_2-a_2 Lx_2-a_2],[Ly_2-b_2 Ly_2+b_2],'m')
% plot([Lx_2-a_2 Lx_2+a_2],[Ly_2-b_2 Ly_2-b_2],'m')
end