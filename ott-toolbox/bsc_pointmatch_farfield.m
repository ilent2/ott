function [nn,mm,a,b] = bsc_pointmatch_farfield( nmax, beam_type, ...
    parameters, varargin )
%BSC_POINTMATCH_FARFIELD calculate beam coefficients for Gaussian beams
%
% [nn,mm,a,b] = BSC_POINTMATCH_FARFIELD(Nmax, type, parameters[, optional])
% Uses an overdetermined point-matching method to find
% spherical harmonic expansion of a laser beam.
%
% Currently available types of beams [parameters]:
% 0 Gauss-Hermite beam
%   [ m n beam_angle P xcomponent ycomponent truncation_angle beam_offset ]
% 1 Laguerre-Gauss beam
%   [ p l beam_angle P xcomponent ycomponent truncation_angle beam_offset ]
% 2 Ince-Gauss beam
%   [ o m p xi beam_angle P xcomponent ycomponent truncation_angle beam_offset ]
%
% beam_angle is the angle of the incoming beam waist.
% beam_offset is a three-column vector [x,y,z] of new beam origin.
%
% Optional parameters:
%
% 'radial' - makes radial component with weighting xcomponent.
% 'azimuthal' - makes azimuthal component with weighting ycomponent. (note:
%   azimuthal and radial components are not mutually exclusive.)
% 'sintheta' - angular scaling function is the same as the one present in
%   standard microscope objectives. Preserves high order mode shape!
% 'tantheta' - default angular scaling function, "small angle
%   approximation" which is valid for thin lenses ONLY. Does not preserve
%   high order mode shape at large angles.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

axisymmetry = 1;
%axisymmetry = 0;

zero_rejection_level = 1e-8;

speed_of_light = 3.00e8;
%medium_refractive_index = 1.33;
medium_refractive_index = 1;
%beam_wavelength0 = 1064e-9;
beam_wavelength0 = 1;

%radial and azimuthal polarisation.
radial=0;
azimuthal=0;

%% these parameters aren't used
beam_power = 1e-3;
speed_in_medium = speed_of_light/medium_refractive_index;
epsilon0 = 8.854187817e-12;
kappa = medium_refractive_index^2;
%%
radial_mode = parameters(1);
azimuthal_mode = parameters(2);

ott_warning('internal');

%% mode selection
switch beam_type
    case 0
        m=radial_mode;
        n=azimuthal_mode;
        paraxial_order=n+m;
        [modeweights,initial_mode,final_mode]=paraxial_transformation_matrix(paraxial_order,0,1,0);
        [row]=find(final_mode(:,1)==m,1);
    case 1
        paraxial_order=2*radial_mode+abs(azimuthal_mode);
        modeweights=eye(paraxial_order+1);
        row=(azimuthal_mode+paraxial_order)/2+1;
        
        i2_out=[-paraxial_order:2:paraxial_order].';
        i1_out=floor((paraxial_order-abs(i2_out))/2);
        
        initial_mode=[i1_out,i2_out];
    case 2
        paraxial_order = radial_mode;
        parity = parameters(3);
        elipticity = parameters(4);
        parameters=parameters(3:end); %pop off the extra parameters.
        
        [modeweights,initial_mode,final_mode]=paraxial_transformation_matrix(paraxial_order,0,[2,elipticity],0);
        
        [row]=find(and(final_mode(:,2)==azimuthal_mode,final_mode(:,3)==parity),1);
        
        if and(paraxial_order>1,isempty(row))
            ott_warning('external');
            error('Observe parity convensions!')
        end
end
% find the mode columns:
keepz=(abs(modeweights(row,:))>0);
initial_mode=initial_mode(keepz,:);
c=modeweights(row,keepz);

beam_angle = parameters(3);
k = 2*pi * medium_refractive_index / beam_wavelength0;
xcomponent = parameters(5);
ycomponent = parameters(6);

% Truncation angle
if length(parameters) > 6
    truncation_angle = parameters(7);
else
    truncation_angle = 90;
end

ott_warning('external');

% Offset of focal point from coordinate origin
if length(parameters) > 9
    offset = parameters(8:10);
    if any(abs(parameters(8:9))>0)
        ott_warning('ott:bsc_pointmatch_farfield:offsets', ...
            ['Beam offsets with x and y components cannot be ' ...
             'axi-symmetric, beam symmetry is now off, and the ' ...
             'calculation will be much slower. It is highly recommended ' ...
             'that a combination of rotations and translations are ' ...
             'used on BSCs instead.']);
        axisymmetry=0;
    end
else
    offset = [];
end

aperture_function=0;
if numel(varargin)>0
    for ii=1:length(varargin)
        
        switch class(varargin{ii})
            case 'char'
                switch lower(varargin{ii})
                    case {'sin','sintheta'}
                        aperture_function=2;
                    case {'tan','tantheta'}
                        aperture_function=0;
                    case {'radial'}
                        radial=1;
                    case {'azimuthal'}
                        azimuthal=1;
                end
            otherwise
                ott_warning('ott:bsc_pointmatch_farfield:input', ...
                    ['Unrecognised input: ' varargin{ii} '.'])
        end
    end
end

ott_warning('internal');

% Grid of points over sphere
ntheta = (nmax + 1);
nphi = 2*(nmax + 1);
if axisymmetry
    ntheta = 2*(nmax+1);
    nphi = 3;
    if beam_type~=1
        nphi = paraxial_order+3-rem(paraxial_order,2);
    end
end

[theta,phi] = angulargrid(ntheta,nphi);

np = length(theta);

% Find electric field at all points
% In the far-field, we have:
% w = 2|z|/(k w0)     (cylindrical coords)
% r/w = kr w0 / 2 |z| (cylindrical coords)
% r = z tan(theta)    (cylindrical -> spherical conversion)
% r/w = k w0 |tan(theta)|/2 (spherical)

%central_irradiance = 2*beam_power / (pi*w0^2);
%central_amplitude = sqrt(2*central_irradiance / ...
%   (speed_in_medium*kappa));

central_amplitude = 1;

w0 = paraxial_beam_waist(paraxial_order);

wscaling=1/tan(abs(beam_angle/180*pi));

rw = 2*(wscaling * w0)^2 * tan(theta).^2 ;
dr = (wscaling * w0) * (sec(theta)).^2 ;

if aperture_function==2
    
    wscaling=1/sin(abs(beam_angle/180*pi));
    
    rw = 2*(wscaling * w0)^2 * sin(theta).^2 ;
    dr = (wscaling * w0) * abs(cos(theta)) ;
end

% degree and order of all modes
total_modes = nmax^2 + 2*nmax;
[nn,mm] = combined_index((1:total_modes)');

mode_index_vector=[];
beam_envelope = zeros(np,length(c));
for ii=1:length(c)
    radial_mode=initial_mode(ii,1);
    azimuthal_mode=initial_mode(ii,2);
    
    norm_paraxial=sqrt(2*factorial(radial_mode)/(pi*factorial(radial_mode+abs(azimuthal_mode))));
    L = laguerre(radial_mode,abs(azimuthal_mode),rw);
    beam_envelope(:,ii) = norm_paraxial.*rw.^abs(azimuthal_mode/2) .* L .* exp(-rw/2 + 1i*azimuthal_mode*phi+1i*pi/2*(radial_mode*2+abs(azimuthal_mode)+1));
    mode_input_power=sqrt(sum(2*pi*abs(beam_envelope(:,ii)).^2.*sqrt(rw/2).*abs(dr)));
    aperture_power_normalization=sqrt(sum(2*pi*abs(beam_envelope(:,ii)).^2.*sin(theta)));
    
    beam_envelope(:,ii)=c(ii)*beam_envelope(:,ii)/aperture_power_normalization*mode_input_power;
    
    mode_index_vector=[mode_index_vector;find(mm==azimuthal_mode+1-max([azimuthal,radial])|mm==azimuthal_mode-1+max([azimuthal,radial]))];

end
mode_index_vector=unique(mode_index_vector);

beam_envelope=sum(beam_envelope,2);
outbeam = find(theta<pi*(180-truncation_angle)/180);
beam_envelope(outbeam) = 0;

if ~isempty(offset)
    rhat = rtpv2xyzv( ones(size(theta)), zeros(size(theta)), zeros(size(theta)), ones(size(theta)), theta, phi );
    [offset,rhat] = matchsize(offset,rhat);
    phase_shift = exp( -1i * k * dot(offset,rhat,2) );
    beam_envelope = beam_envelope .* phase_shift;
end
Ex = xcomponent * beam_envelope * central_amplitude;
Ey = ycomponent * beam_envelope * central_amplitude;

if any(azimuthal|radial)
    Etheta=-radial*xcomponent*beam_envelope * central_amplitude;
    Ephi=azimuthal*ycomponent*beam_envelope * central_amplitude;
else
    Etheta = - Ex .* cos(phi) - Ey .* sin(phi);
    Ephi = - Ex .* sin(phi) + Ey .* cos(phi);
end

e_field = [[ Etheta(:); Ephi(:) ]];

if axisymmetry
    nn=nn(mode_index_vector);
    mm=mm(mode_index_vector);
    
    removeels=find(abs(mm)>paraxial_order+1);
    nn(removeels)=[];
    mm(removeels)=[];
end

coefficient_matrix = zeros(2*np,2*length(nn));

for n = 1:max(nn)
    ci=find(nn==n);
    
    [~,dtY,dpY]= spharm(n,mm(ci),theta,phi);
    
    coefficient_matrix(:,ci) = [dpY;-dtY] * 1i^(n+1)/sqrt(n*(n+1));
    coefficient_matrix(:,ci+length(nn)) = [dtY;dpY] * 1i^(n)/sqrt(n*(n+1));
    
end

a=zeros(size(nn));
b=zeros(size(nn));

expansion_coefficients = coefficient_matrix \ e_field;
%fprintf(1,'done!\n');
%toc
a = expansion_coefficients(1:end/2,:);
b = expansion_coefficients(1+end/2:end,:);

p=abs(a).^2+abs(b).^2;
binaryvector=(p>zero_rejection_level*max(p));

if nargout>2
    nn=nn(binaryvector);
    mm=mm(binaryvector);
    a=a(binaryvector);
    b=b(binaryvector);
end

if nargout==2
    ci=combined_index(nn,mm);
    
    mm=sparse(ci,1,b,nmax*(nmax+2),1);
    nn=sparse(ci,1,a,nmax*(nmax+2),1);
end

ott_warning('external');

return

