function T = cylinderTMatrix(Nmax, k_medium, k_particle, profile, varargin)
%cylinderTMatrix Calculate T-matrix for a cylinder.
%
%   Includes options to calculate T-matrix for elongated sphere,
%   cylinder and cylinder with bevelled/rounded edges.
%
%   Uses either PM, EBCM or DDA.  If method not explicitly specified,
%   will use results of Qi et al., 2014 to calculate the T-matrix at
%   the specified accuracy using the fastest method.
%
%   Requires: Optical tweezers toolbox.
%
%   Copyright Isaac Lenton 2018

% Load optical tweezers path
addpath('/home/ilenton/Documents/Study/OtToolbox');

%% Parse inputs

p = inputParser;
p.addParameter('method', 'qi2014');
p.addParameter('tol', 0.1);
p.parse(varargin{:});
pResults = p.Results;

pResults.Nmax = Nmax;
pResults.k_medium = k_medium;
pResults.k_particle = k_particle;

%% Determine which method to use

if strcmp(pResults.method, 'ebcm') ...
        || strcmp(pResults.method, 'pm') || strcmp(pResults.method, 'dda')
    % Nothing to do
elseif strcmp(pResults.method, 'qi2014')
    
    % EBCM 1% data
    ebcm1 = {};
    ebcm1.x = [73, 169, 198, 228, 261, 391, 586, 718, 718, 657, 523, 457, 262, 73];
    ebcm1.y = [409, 406, 418, 423, 397, 412, 400, 375, 223, 193, 195, 165, 204, 390];
    ebcm1.x = (ebcm1.x - 73) * 2.0 / (718 - 73);
    ebcm1.y = -(ebcm1.y - 438) * 6.0 / (438 - 9);

    % PM 1% data
    pm1 = {};
    pm1.x = [297, 355, 394, 718, 718, 591, 525, 391, 361, 297];
    pm1.y = [943, 933, 946, 894, 868, 846, 874, 864, 913, 913];
    pm1.x = (pm1.x - 73) * 2.0 / (718 - 73);
    pm1.y = -(pm1.y - 985) * 6.0 / (985 - 555);

    % EBCM 10% data
    ebcm10 = {};
    ebcm10.x = [73, 193, 718, 718, 525, 328, 229, 160, 73];
    ebcm10.y = [430, 426, 381, 37, 94, 177, 214, 274, 375];
    ebcm10.x = (ebcm10.x - 73) * 2.0 / (718 - 73);
    ebcm10.y = -(ebcm10.y - 438) * 6.0 / (438 - 9);

    % PM 10% data
    pm10 = {};
    pm10.x = [130, 160, 328, 397, 462, 522, 589, 718, 718, 654, 589, 522, 328, 265, 130];
    pm10.y = [961, 970, 967, 951, 946, 946, 925, 912, 753, 784, 798, 798, 865, 874, 948];
    pm10.x = (pm10.x - 73) * 2.0 / (718 - 73);
    pm10.y = -(pm10.y - 985) * 6.0 / (985 - 555);
    
    % Conversion factor, paper uses 1064nm illumination
    k = pResults.k_medium * 1.064 / 2.0 / pi;
    diameter = 2.0 * max(profile.rho) * k;
    len = 2.0 * max(profile.z) * k;
    
    if pResults.tol == 0.1
        if inpolygon(diameter, len, pm10.x, pm10.y)
            pResults.method = 'pm';
        elseif inpolygon(diameter, len, ebcm10.x, ebcm10.y)
            pResults.method = 'ebcm';
        else
            pResults.method = 'dda';
        end
    elseif pResults.tol == 0.01
        if inpolygon(diameter, len, pm1.x, pm1.y)
            pResults.method = 'pm';
        elseif inpolygon(diameter, len, ebcm1.x, ebcm1.y)
            pResults.method = 'ebcm';
        else
            pResults.method = 'dda';
        end
    else
        error('method=qi2014 only accepts tol=0.1 or tol=0.01');
    end
else
    error('Unknown method selected');
end

%% Calculate the T-matrix

% Assume the object is a cylinder and calculate a suitable Nmax

if strcmp(pResults.method, 'ebcm')
    T = tmatrix_ebcm_axisym(pResults.Nmax, pResults.k_medium, ...
            pResults.k_particle, profile.rho, profile.z);
elseif strcmp(pResults.method, 'pm')
    T = tmatrix_pm_axisym(pResults.Nmax, pResults.k_medium, ...
            pResults.k_particle, profile.rho, profile.z);
elseif strcmp(pResults.method, 'dda')
    % TODO
    error('dda not yet implemented');
end

end