function [trd, plg] = my_Cart2Sph(cn,ce,cd)
% Converts from cartesion to spherical coordinates
%   [trd,plg] = my_Cart2Sph(cn,ce,cd) returns the trend and plunge of a
%   line for input north, east, and down direction cosines

plg = asin(cd);
if cn == 0
    if ce < 0
        trd = 3/2*pi;
    else
        trd = pi/2;
    end
    
else
    trd = atan(ce./cn);
    if cn < 0
        trd = trd+pi;
    end
    trd = ZeroTwoPi(trd);
    trd= rad2deg(trd);
    
end
