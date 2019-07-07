classdef Vector
   properties
      x
      y
      icx
      icy
   end
   methods
      function obj = Vector(x, y, icx, icy)
         if nargin == 2
            obj.icx = 0;
            obj.icy = 0;
         elseif nargin == 4
            obj.icx = icx;
            obj.icy = icy;     
         end
         obj.x = x;
         obj.y = y;
      end
   end
end

