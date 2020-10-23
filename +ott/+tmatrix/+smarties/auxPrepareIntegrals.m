function stRtfunc = auxPrepareIntegrals(nNint, sInt)
%% auxPrepareIntegrals
% Calculates points and weights for integral quadrature
%
%	auxPrepareIntegrals(nNint, sInt) calculates the points and weights for
%	integral quadrature, using the given method.
%
%	Input:
%	- nNint: the number of integration points required
%	- sInt: the integration scheme, either 'gauss' for Gaussian quadrature, or
%	'rectangle' for Simpson.
%
%	Output:
%	A structure with fields nNbTheta (number of theta), theta, and wTheta
%	(weights)
%  For gauss quadrature, precalculated values stored in the file
%  Utils/quadTable.mat will be used if available. Those values include
%  nNint in the range 50 to 500 by steps of 5, and from 600 to 2000 by
%  steps of 100 with pairs 600 605 700 705 ... 2000 2005 to facilitate
%  convergence tests
%
% Dependency:
% auxInitLegendreQuad

import ott.tmatrix.smarties.*;

% Defines integration types
switch lower(sInt)
    case 'gauss'
        % For Gauss-Legendre quadrature
        stRtfunc.nNbTheta = nNint;

        if(nNint < 50)
            % need to calculate those
            % disp 'case 1'
            [xi,wi]=auxInitLegendreQuad(nNint);
            stRtfunc.theta=acos(xi); % [T x 1]
            stRtfunc.wTheta=wi; % [T x 1]
        elseif(exist([fileparts(which('storeGLquadrature.m')), '/quadTable.mat'], 'file') == 2)
            load([fileparts(which('storeGLquadrature.m')), '/quadTable.mat'])
            Nt = quadTable.Nt;
            if(any(nNint == Nt))
                % disp 'case 2'
                % use tabulated values
                ind = Nt == nNint;
                xw = quadTable.values{ind};
                stRtfunc.theta=xw(:,1); % [T x 1]
                stRtfunc.wTheta=xw(:,2); % [T x 1]
            else % no tabulated value for this N
                % disp 'case 3'
                [xi,wi]=auxInitLegendreQuad(nNint);
                stRtfunc.theta=acos(xi); % [T x 1]
                stRtfunc.wTheta=wi; % [T x 1]
            end
        else  % file failed to load
            % disp 'case 4'
            [xi,wi]=auxInitLegendreQuad(nNint);
            stRtfunc.theta=acos(xi); % [T x 1]
            stRtfunc.wTheta=wi; % [T x 1]
        end

        % faster algorithm using asymptotic expansion
        % but needs to be checked (slightly different values)
        %  elseif(nNint >= 100)
        %  [~,wi,~,ti]=legpts(nNint, [-1, 1]);
        %  stRtfunc.theta=ti; % [T x 1]
        %  stRtfunc.wTheta=wi'; % [T x 1]
        %  end

    case 'rectangle'
        % For Simpson integration
        nNbTheta=nNint;
        stRtfunc.nNbTheta = nNbTheta;
        stRtfunc.theta=transpose(linspace(0,pi,nNbTheta)); % [T x 1]
        dtheta=pi/(nNbTheta-1);
        stRtfunc.wTheta=dtheta*sin(stRtfunc.theta); % [T x 1]

    otherwise
        disp 'Integration type not recognized'
end
