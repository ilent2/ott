classdef BscPmGauss < ott.BscPointmatch
%BscPmGauss provides HG, LG and IG beams using point matching method
%
% Properties
%  - gtype          --  Type of beam ('gaussian', 'lg', 'hg', or 'ig')
%  - mode           --  Beam modes (2 or 4 element vector)
%  - polarisation   --  Beam polarisation
%  - truncation_angle -- Truncation angle for beam [rad]
%  - offset         --  Offset for original beam calculation
%  - angle          --  Angle of incoming beam waist
%  - angular_scaling -- Angular scaling function (tantheta | sintheta)
%  See :class:`+ott.Bsc` for inherited properties.
%
% This class is based on ``bsc_pointmatch_farfield.m`` and
% ``bsc_pointmatch_focalplane.m`` from OTT (version 1).
%
% See also BscPmGauss.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    gtype              % Type of beam ('gaussian', 'lg', 'hg', or 'ig')
    mode               % Beam modes (2 or 4 element vector)
    polarisation       % Beam polarisation
    truncation_angle   % Truncation angle for beam [rad]
    offset             % Offset for original beam calculation
    angle              % Angle of incoming beam waist
    translation_method % Translation method to use
    angular_scaling    % Angular scaling function (tantheta | sintheta)
  end

  % TODO: Incorperate bsc_pointmatch_focalplane option for gaussian beams.

  methods (Static)
    function l = supported_beam_type(s)
      l = strcmp(s, 'lg') || strcmp(s, 'hg') || strcmp(s, 'ig');
    end
    
    function validate_translation_method(method)
      assert(any(strcmpi(method, {'Default', 'NewBeamOffset'})), ...
        'Translation method must be Default or NewBeamOffset');
    end
  end

  methods
    function beam = BscPmGauss(varargin)
      % Construct a new IG, HG, LG or Gaussian beam.
      %
      % Usage
      %   BscPmGauss(...) constructs a new Gassian beam (LG00).
      %
      %   BscPmGauss(type, mode, ...) constructs a new beam with the given type.
      %   Supported types [mode]:
      %    - 'lg' -- Laguarre-Gauss  [ radial azimuthal ]
      %    - 'hg' -- Hermite-Gauss   [ m n ]
      %    - 'ig' -- Ince-Gauss      [ paraxial azimuthal parity elipticity ]
      %
      % Optional named parameters
      %   - 'Nmax'  --          Truncation number for beam shape coefficients.
      %     If omitted, Nmax is initially set to 100, the beam is
      %     calculated and Nmax is reduced so that the power does
      %     not drop significantly.
      %   - 'zero_rejection_level' -- Level used to determine non-zero
      %     beam coefficients in far-field point matching.  Default: 1e-8.
      %
      %   - 'NA'     --         Numerical aperture of objective
      %   - 'polarisation'  --  Polarisation of the beam
      %   - 'power'         --  Rescale the power of the beam (default: [])
      %
      %   - 'omega'         --  Optical angular frequency (default: 2*pi)
      %   - 'k_medium'      --  Wave number in medium
      %   - 'index_medium'  --  Refractive index of medium
      %   - 'wavelength_medium' -- Wavelength in medium
      %
      %   - 'wavelength0'   --  Wavelength in vacuum
      %   - 'offset'        --  Offset of the beam from origin
      %
      %   - translation_method -- Method to use when calculating translations.
      %     Can either be 'Default' or 'NewBeamOffset', the latter calculates
      %     new beam shape coefficients for every new position.
      %
      %   - angular_scaling (enum) -- Angular scaling function.
      %     For a discussion of this parameter, see Documentation
      %     (:ref:`conception-angular-scaling`).
      %      - 'sintheta' -- angular scaling function is the same as the
      %        one present in standard microscope objectives.
      %        Preserves high order mode shape!
      %      - 'tantheta' -- default angular scaling function,
      %        "small angle approximation" which is valid for thin
      %        lenses ONLY. Does not preserve high order mode shape
      %        at large angles.
      %
      %   - truncation_angle (numeric) -- Adds a hard edge to the
      %     beam, this can be useful for simulating the back-aperture
      %     of a microscope objective.  Default: ``pi/2`` (i.e. no edge).
      %
      %   - truncation_angle_deg -- Same as `truncation_angle` but
      %     with degrees instead of radians.

      beam = beam@ott.BscPointmatch(varargin{:});
      beam.type = 'incident';
      beam.basis = 'regular';

      % Parse inputs
      p = inputParser;
      p.addOptional('type', 'lg', @ott.BscPmGauss.supported_beam_type);
      p.addOptional('mode', [ 0 0 ]);

      p.addParameter('Nmax', []);
      p.addParameter('zero_rejection_level', 1e-8);
      p.addParameter('offset', []);
      p.addParameter('polarisation', [ 1 1i ]);
      p.addParameter('wavelength0', 1);
      p.addParameter('power', []);
      p.addParameter('progress_callback', []);
      p.addParameter('translation_method', 'Default', ...
        @ott.BscPmGauss.validate_translation_method);

      p.addParameter('omega', 2*pi);
      p.addParameter('k_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('wavelength_medium', []);

      p.addParameter('NA', []);
      p.addParameter('angle_deg', []);
      p.addParameter('angle', []);

      p.addParameter('truncation_angle_deg', []);
      p.addParameter('truncation_angle', []);
      p.addParameter('angular_scaling', 'tantheta');

      p.parse(varargin{:});

      % Store parameters
      beam.gtype = p.Results.type;
      beam.mode = p.Results.mode;
      beam.polarisation = p.Results.polarisation;
      beam.offset = p.Results.offset;
      beam.k_medium = ott.Bsc.parser_k_medium(p, 2*pi);
      beam.omega = p.Results.omega;
      
      % Ensure beam offset is not empty
      if isempty(beam.offset)
        beam.offset = [0;0;0];
      end

      % Store truncation angle
      if isempty(p.Results.truncation_angle_deg) &&  ...
          isempty(p.Results.truncation_angle)
        beam.truncation_angle = pi/2;
      elseif ~isempty(p.Results.truncation_angle)
        beam.truncation_angle = p.Results.truncation_angle;
      elseif ~isempty(p.Results.truncation_angle_deg)
        beam.truncation_angle = p.Results.truncation_angle_deg * pi/180;
      else
        ott.warning('external');
        error('Truncation angle given in degrees and radians');
      end
      
      beam.translation_method = p.Results.translation_method;

      % TODO: bsc_pointmatch_farfield.m had other arguments
      % optional parameters:
      %
      % 'radial' - makes radial component with weighting xcomponent.
      % 'azimuthal' - makes azimuthal component with weighting ycomponent.
      % (note: azimuthal and radial components are not mutually exclusive.)

      import ott.utils.*

      axisymmetry = 1;

      ott.warning('internal');

      %radial and azimuthal polarisation.
      radial=0;
      azimuthal=0;

      % mode selection
      switch p.Results.type
        case 'hg'
          assert(numel(p.Results.mode) == 2, ...
            'ott:BscPmGauss:wrong_mode_length', ...
            'mode must be 2 element vector');

          m = p.Results.mode(1);
          n = p.Results.mode(2);
          paraxial_order=n+m;
          [modeweights,initial_mode,final_mode] = ...
              paraxial_transformation_matrix(paraxial_order,0,1,0);
          [row]=find(final_mode(:,1)==m,1);
        case 'lg'

          assert(numel(p.Results.mode) == 2, ...
            'ott:BscPmGauss:wrong_mode_length', ...
            'mode must be 2 element vector');

          radial_mode = p.Results.mode(1);
          azimuthal_mode = p.Results.mode(2);
          assert(radial_mode - floor(radial_mode) == 0 && radial_mode >= 0, ...
            'ott:BscPmGauss:invalid_radial_mode', ...
            'Radial mode index must be positive integer');
          assert(azimuthal_mode - floor(azimuthal_mode) == 0, ...
            'ott:BscPmGauss:invalid_azimuthal_mode', ...
            'Azimuthal mode index must be integer');

          paraxial_order=2*radial_mode+abs(azimuthal_mode);
          modeweights=eye(paraxial_order+1);
          row=(azimuthal_mode+paraxial_order)/2+1;
          
          i2_out= (-paraxial_order:2:paraxial_order).';
          i1_out=floor((paraxial_order-abs(i2_out))/2);
          
          initial_mode=[i1_out,i2_out];
        case 'ig'
          assert(numel(p.Results.mode) == 4, ...
            'ott:BscPmGauss:wrong_mode_length', ...
            'mode must be 4 element vector');

          paraxial_order = p.Results.mode(1);
          azimuthal_mode = p.Results.mode(2);
          parity = p.Results.mode(3);
          elipticity = p.Results.mode(4);
          
          [modeweights,initial_mode,final_mode] = ...
             paraxial_transformation_matrix(paraxial_order,0,[2,elipticity],0);
          
          [row]=find(and(final_mode(:,2)==azimuthal_mode, ...
              final_mode(:,3)==parity),1);
          
          if and(paraxial_order>1,isempty(row))
            ott.warning('external');
            error('Observe parity convensions!')
          end
      end
      % find the mode columns:
      keepz=(abs(modeweights(row,:))>0);
      initial_mode=initial_mode(keepz,:);
      c=modeweights(row,keepz);

      beam_angle_specified = ~isempty(p.Results.angle) ...
          + ~isempty(p.Results.angle_deg) + ~isempty(p.Results.NA);

      % Find beam_angle_deg
      if beam_angle_specified > 1
        ott.warning('external');
        error('Too many inputs.  Only specify NA, angle or angle_deg');
      elseif isempty(p.Results.angle) && isempty(p.Results.angle_deg)
        if isempty(p.Results.NA)
          NA = 1.02;
          index = 1.33;
        else
          NA = p.Results.NA;
          if ~isempty(p.Results.index_medium)
            index = p.Results.index_medium;
          else
            ott.warning('external');
            error('Need to specify index_medium with NA');
          end
        end
        beam_angle_deg = asin(NA/index)*180.0/pi;
      elseif ~isempty(p.Results.angle_deg)
        beam_angle_deg = p.Results.angle_deg;
      elseif ~isempty(p.Results.angle)
        beam_angle_deg = p.Results.angle*180/pi;
      end
      beam.angle = beam_angle_deg * pi/180;

      xcomponent = p.Results.polarisation(1);
      ycomponent = p.Results.polarisation(2);
      offset = p.Results.offset;

      if numel(offset) == 3 && any(abs(offset(1:2))>0)
        
        % Only warn if using beams that support matrix translations
        if strcmpi(p.Results.translation_method, 'Default')
          ott.warning('external');
          ott.warning('ott:bsc_pointmatch_farfield:offsets', ...
              ['Beam offsets with x and y components cannot be ' ...
               'axi-symmetric, beam symmetry is now off, and the ' ...
               'calculation will be much slower. It is highly recommended ' ...
               'that a combination of rotations and translations are ' ...
               'used on BSCs instead.']);
            ott.warning('internal');
        end
        
        % Turn off axissymmetry
        axisymmetry=0;
      end
      
      w0 = ott.utils.paraxial_beam_waist(paraxial_order);
      wscaling=1/tan(abs(beam_angle_deg/180*pi));

      % Store or calculate Nmax
      if isempty(p.Results.Nmax)
        nmax = 100;
      else
        nmax = p.Results.Nmax;
      end

      % Grid of points over sphere
      ntheta = (nmax + 1);
      nphi = 2*(nmax + 1);
      if axisymmetry
          ntheta = 2*(nmax+1);
          nphi = 3;
          if ~strcmp(beam.gtype, 'lg')
              nphi = paraxial_order+3-rem(paraxial_order,2);
          end
      end
      
      % If we have an offset, we need to have a high enough resolution
      % to match the phase gradient across the far-field
      if ~isempty(p.Results.offset)
        offset_lambda = vecnorm(p.Results.offset)*beam.k_medium/(2*pi);
        ntheta = max(ntheta, 3*ceil(offset_lambda));
        nphi = max(nphi, 2*3*ceil(offset_lambda));
      end

      [theta,phi] = ott.utils.angulargrid(ntheta,nphi);

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
      rw = 2*(wscaling * w0)^2 * tan(theta).^2 ;
      dr = (wscaling * w0) * (sec(theta)).^2 ;
      
      if strcmpi(p.Results.angular_scaling, 'tantheta')
        % Nothing to do
      elseif strcmpi(p.Results.angular_scaling, 'sintheta')

        wscaling=1/sin(abs(beam_angle_deg/180*pi));

        rw = 2*(wscaling * w0)^2 * sin(theta).^2 ;
        dr = (wscaling * w0) * abs(cos(theta)) ;
        
      else
        error('Unknown angular_scaling parameter value');
      end
      
      beam.angular_scaling = p.Results.angular_scaling;

      % degree and order of all modes
      total_modes = nmax^2 + 2*nmax;
      [nn,mm] = ott.utils.combined_index((1:total_modes)');

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
          
          mode_index_vector=[mode_index_vector; ...
              find(mm==azimuthal_mode+1-max([azimuthal,radial]) ...
              | mm==azimuthal_mode-1+max([azimuthal,radial]))];

      end
      mode_index_vector=unique(mode_index_vector);

      beam_envelope=sum(beam_envelope,2);
      outbeam = theta < pi-beam.truncation_angle;
      beam_envelope(outbeam) = 0;

      if ~isempty(offset)
        rhat = rtpv2xyzv( ones(size(theta)), zeros(size(theta)), ...
            zeros(size(theta)), ones(size(theta)), theta, phi );
        [offset,rhat] = matchsize(offset.',rhat);
        phase_shift = exp( 1i * beam.k_medium * dot(offset,rhat,2) );
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

      e_field = [ Etheta(:); Ephi(:) ];

      if axisymmetry
        nn=nn(mode_index_vector);
        mm=mm(mode_index_vector);

        removeels=find(abs(mm)>paraxial_order+1);
        nn(removeels)=[];
        mm(removeels)=[];
      end

      % Do the point matching and store the result
      [beam.a, beam.b] = beam.bsc_farfield(nn, mm, e_field, theta, phi, ...
        'zero_rejection_level', p.Results.zero_rejection_level);

      % If no Nmax supplied, shrink the beam to the smallest size that
      % preserves the beam power
      if isempty(p.Results.Nmax)
        beam = beam.shrink_Nmax();
      end

      % Normalize the beam power
      if ~isempty(p.Results.power)
        beam.power = p.Results.power;
      end

      ott.warning('external');
    end
    
    function varargout = translateZ(beam, varargin)
      %TRANSLATEZ translate a beam along the z-axis
      %
      % beam = TRANSLATEZ(z) translates by a distance z along the z axis.
      % If the translation_method is NewBeamOffset, a new beam is generated.
      %
      % [beam, A, B] = TRANSLATEZ(z) returns the translation matrices
      % and the translated beam.  See also Bsc.TRANSLATE.
      %
      % [beam, AB] = TRANSLATEZ(z) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATEZ(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.
      
      if strcmpi(beam.translation_method, 'Default')
        
        % Use translation matrix method
        [varargout{1:nargout}] = translateZ@ott.BscPointmatch(beam, varargin{:});
        
      elseif strcmpi(beam.translation_method, 'NewBeamOffset')
        
        p = inputParser;
        p.addOptional('z', []);
        p.addParameter('Nmax', beam.Nmax);
        p.parse(varargin{:});
        
        % Generate the new beam
        varargout{1} = ott.BscPmGauss(beam.gtype, beam.mode, ...
          'offset', beam.offset + [0;0;p.Results.z], ...
          'omega', beam.omega, 'power', beam.power, ...
          'wavelength_medium', beam.wavelength, ...
          'polarisation', beam.polarisation, ...
          'truncation_angle', beam.truncation_angle, ...
          'Nmax', p.Results.Nmax, 'angle', beam.angle);
        varargout{1}.type = beam.type;
        varargout{1}.basis = beam.basis;
      end
    end
    
    function varargout = translateXyz(beam, varargin)
      %TRANSLATEXYZ translate the beam given Cartesian coordinates
      %
      % beam = TRANSLATEXYZ(xyz) translate the beam to locations given by
      % the xyz coordinates, where xyz is a 3xN matrix of coordinates.
      % If the translation_method is NewBeamOffset, a new beam is generated.
      %
      % TRANSLATEXYZ(Az, Bz, D) translate the beam using
      % z-translation and rotation matricies.
      %
      % [beam, Az, Bz, D] = TRANSLATEXYZ(...) returns the z-translation
      % matrices, the rotation matrix D, and the translated beam.
      %
      % [beam, A, B] = TRANSLATEXYZ(...) returns the translation matrices
      % and the translated beam.
      %
      % [beam, AB] = TRANSLATEXYZ(...) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATEXYZ(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.

      if strcmpi(beam.translation_method, 'Default')
        
        % Use translation matrix method
        [varargout{1:nargout}] = translateXyz@ott.BscPointmatch(beam, varargin{:});
        
      elseif strcmpi(beam.translation_method, 'NewBeamOffset')
        
        p = inputParser;
        p.addOptional('opt1', []);    % xyz or Az
        p.addOptional('opt2', []);    % [] or Bz
        p.addOptional('opt3', []);    % [] or D
        p.addParameter('Nmax', beam.Nmax);
        p.parse(varargin{:});
        
        assert(isempty(p.Results.opt2) && isempty(p.Results.opt3), ...
          'Rotation and translation matries not supported with this method');
        
        % Generate the new beam
        varargout{1} = ott.BscPmGauss(beam.gtype, beam.mode, ...
          'offset', beam.offset + p.Results.opt1, ...
          'omega', beam.omega, 'power', beam.power, ...
          'wavelength_medium', beam.wavelength, ...
          'polarisation', beam.polarisation, ...
          'truncation_angle', beam.truncation_angle, ...
          'Nmax', p.Results.Nmax, 'angle', beam.angle);
        varargout{1}.type = beam.type;
        varargout{1}.basis = beam.basis;
      end
    end
    
    function varargout = translateRtp(beam, varargin)
      %TRANSLATERTP translate the beam given spherical coordinates
      %
      % beam = TRANSLATERTP(rtp) translate the beam to locations given by
      % the xyz coordinates, where rtp is a 3xN matrix of coordinates.
      % If the translation_method is NewBeamOffset, a new beam is generated.
      %
      % TRANSLATERTP(Az, Bz, D) translate the beam using
      % z-translation and rotation matricies.
      %
      % [beam, Az, Bz, D] = TRANSLATERTP(...) returns the z-translation
      % matrices, the rotation matrix D, and the translated beam.
      %
      % [beam, A, B] = TRANSLATERTP(...) returns the translation matrices
      % and the translated beam.
      %
      % [beam, AB] = TRANSLATERTP(...) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATERTP(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.
      
      if strcmpi(beam.translation_method, 'Default')
        
        % Use translation matrix method
        [varargout{1:nargout}] = translateRtp@ott.BscPointmatch(beam, varargin{:});
        
      elseif strcmpi(beam.translation_method, 'NewBeamOffset')
        
        p = inputParser;
        p.addOptional('opt1', []);    % xyz or Az
        p.addOptional('opt2', []);    % [] or Bz
        p.addOptional('opt3', []);    % [] or D
        p.addParameter('Nmax', beam.Nmax);
        p.parse(varargin{:});
        
        assert(isempty(p.Results.opt2) && isempty(p.Results.opt3), ...
          'Rotation and translation matries not supported with this method');
        
        % Generate the new beam
        xyz = ott.utils.rtp2xyz(p.Results.opt1);
        varargout{1} = ott.BscPmGauss(beam.gtype, beam.mode, ...
          'offset', beam.offset + xyz(:), ...
          'omega', beam.omega, 'power', beam.power, ...
          'wavelength_medium', beam.wavelength, ...
          'polarisation', beam.polarisation, ...
          'truncation_angle', beam.truncation_angle, ...
          'Nmax', p.Results.Nmax, 'angle', beam.angle);
        varargout{1}.type = beam.type;
        varargout{1}.basis = beam.basis;
      end
    end
  end
end

