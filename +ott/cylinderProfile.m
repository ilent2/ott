function profile = cylinderProfile(length, diameter, corner, cradius, cres)
%cylinderProfile Summary of this function goes here
%   Detailed explanation goes here

profile = {};

if nargin == 2
    profile.rho = [0, diameter/2, diameter/2, 0];
    profile.z = [-length/2, -length/2, length/2, length/2];
elseif nargin >= 4
    if cradius > diameter/2 || cradius <= 0
        error('Invalud value, make sure 0 < cradius <= diameter/2');
    elseif cradius > length/2
        error('Invalud value, make sure 0 < cradius <= length/2');
    end
    
    if corner == 'bevel'
        profile.rho = [0, diameter/2-cradius, diameter/2, ...
            diameter/2, diameter/2-cradius, 0];
        profile.z = [-length/2, -length/2, -length/2 + cradius, ...
            length/2 - cradius, length/2, length/2];
    elseif corner == 'round'
        profile.rho = [0, diameter/2-cradius];
        profile.z = [-length/2, -length/2];
        for ii = 1:cres
            theta = pi/2*ii/cres;
            profile.rho(2+ii) = diameter/2-cradius + cradius*sin(theta);
            profile.z(2+ii) = -length/2 + cradius - cradius*cos(theta);
        end
        profile.rho = [profile.rho, flip(profile.rho)];
        profile.z = [profile.z, -flip(profile.z)];
    else
        error('Unrecognized corner type');
    end
else
    error('Too few arguments');
end

end

