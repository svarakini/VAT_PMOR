function theta = VAT_fiber_ply_angle_1D(T0,T1,x,center,width)

% use coordinate to identify the fiber ply orientation
%
%

% size(T0)
% size(T1)


xbar = (x-center(1))/(width/2);
% size(xbar)

if xbar>=0
    
    theta = T0 + xbar.*(T1-T0);
    
else
    
    
    theta = T0 - xbar.*(T1-T0);
end