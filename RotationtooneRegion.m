function [X_e1,Y_e1,AngleB] = RotationtooneRegion(X_e,Y_e,X_o,Y_o,AngleB)
if AngleB == -1
    Angle=calculate_vector_angle(X_e,Y_e, X_o, Y_o); %[MODIFIED]
    SingAngle=Angle/abs(Angle);

    if (Angle<90 && SingAngle>=0 )|| (Angle<-270 && SingAngle<0)% POINNT PPOINT You can't make this kind of comparison!!! URGENT
      %   1
     X_e1=X_e;
     Y_e1=Y_e;
     %No rotation X_e , Y_e
     elseif (Angle>=90 && Angle<=180 && SingAngle>=0) || (Angle>=-270 && Angle<=-180 && SingAngle<0)
         if SingAngle>=0
             AngleB=90-155;
         else
             AngleB=-90-155;
         end
         X_e1=(X_e-X_o)*cos((pi/180)*(+AngleB))-(Y_e-Y_o)*sin((pi/180)*(+AngleB))+X_o;
         Y_e1=(X_e-X_o)*sin((pi/180)*(+AngleB))+(Y_e-Y_o)*cos((pi/180)*(+AngleB))+Y_o;
     elseif (Angle>180 && Angle<=270 && SingAngle>=0) || (Angle>-180 && Angle<=-90 && SingAngle<0)
         if SingAngle>=0
             AngleB=180-155;
         else
             AngleB=-180-155;
         end

         X_e1=(X_e-X_o)*cos((pi/180)*(+AngleB))-(Y_e-Y_o)*sin((pi/180)*(+AngleB))+X_o;
         Y_e1=(X_e-X_o)*sin((pi/180)*(+AngleB))+(Y_e-Y_o)*cos((pi/180)*(+AngleB))+Y_o;
     else
         if SingAngle>=0
             AngleB=270-155;
         else
             AngleB=-270-155;
         end
         X_e1=(X_e-X_o)*cos((pi/180)*(+AngleB))-(Y_e-Y_o)*sin((pi/180)*(+AngleB))+X_o;
         Y_e1=(X_e-X_o)*sin((pi/180)*(+AngleB))+(Y_e-Y_o)*cos((pi/180)*(+AngleB))+Y_o;
     end
     
    else %Reverse the Operations
        Angle=calculate_vector_angle(X_e,Y_e, X_o, Y_o); %[MODIFIED]
        SingAngle=Angle/abs(Angle);
        AngleB=-AngleB;
        X_e1=(X_e-X_o)*cos((pi/180)*(+AngleB))-(Y_e-Y_o)*sin((pi/180)*(+AngleB))+X_o;
        Y_e1=(X_e-X_o)*sin((pi/180)*(+AngleB))+(Y_e-Y_o)*cos((pi/180)*(+AngleB))+Y_o;
    end
end
