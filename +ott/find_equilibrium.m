function eq = find_equilibrium(z, fz)
%FIND_EQUILIBRIUM estimates equilibrium positions from position-force data
%
% Based on the code in example_gaussian from the original toolbox.
%
% TODO: Generalize the code to find multiple equilibriums.
% TODO: z need not be a vector of scalars, we could have a array of
%   position vectors representing some path we want to find the
%   equilibrium on.  We could do a similar thing for fz.

% This function is not directly concerned with force/torque calculation
warning('This function will move in a future release');

zeroindex=find(fz<0,1);

if length(zeroindex)~=0
    %fit to third order polynomial the local points. (only works when dz
    %sufficiently small)
    zrange = max([zeroindex-2,1]):min([zeroindex+2,length(z)]);
    pz=polyfit(z(zrange), fz(zrange), 3);
    root_z=roots(pz); %find roots of 3rd order poly.
    dpz=[3*pz(1),2*pz(2),1*pz(3)]; %derivative of 3rd order poly.

    real_z=root_z(imag(root_z)==0); % finds real roots only.

    rootsofsign=polyval(dpz,real_z); %roots that are stable
    eq=real_z(rootsofsign<0); %there is at most 1 stable root. critical roots give error.
    try
        eq=zeq(abs(zeq-z(zeroindex))==min(abs(zeq-z(zeroindex))));
    end
else
    eq=[];
end

end
