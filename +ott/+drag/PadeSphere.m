classdef PadeSphere < ott.drag.StokesSphereWall
% Creeping flow around a sphere in shear flow near to a wall.
% Inherits from :class:`StokesSphereWall`.
%
% This class implements the Pade series approximation
% for spherical particles moving near a planar surface.
% The approximation should work for spacing between the sphere surface
% and plane larger than :math:`10^{-2}\times`radius.
%
% This class uses the coefficients included in
%
%   M. Chaoui and F. Feuillebois,
%   Creeping Flow Around a Sphere in a Shear Flow Close to a Wall.
%   The Quarterly Journal of Mechanics and Applied Mathematics,
%   Volume 56, Issue 3, August 2003, Pages 381--410,
%   https://doi.org/10.1093/qjmam/56.3.381
%
% This class assumes the surface is perpendicular to the z axis, positioned
% bellow the spherical particle.
%
% Properties
%   - radius      -- Radius of the sphere
%   - viscosity   -- Viscosity of the medium
%   - separation  -- Distance between the sphere centre and plane
%   - forward     -- Computed drag tensor
%   - inverse     -- Computed from `forward`.
%
% See :class:`Stokes` for other methods/parameters.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    forwardInternal
  end

  properties (Hidden, Constant)
    % Table 9 from the text
    table9 = [
      0 1.0000000000000000 1.0000000000000000
      1 0.4908044826583015 -0.0716955173416986
      2 -5.9089832134410898 -6.1850607349363846
      3 -3.5131616290239789 -0.0643585714633657
      4 15.8084387577986870 17.6701614844569620
      5 11.3026454259591860 1.5230014975500095
      6 -25.4204600037753040 -31.0142617248087940
      7 -22.1471886954821320 -4.9543202018606971
      8 27.7560683302883220 37.7762831912656620
      9 29.9526704556636040 8.6523092950119675
      10 -22.5314407803093640 -34.4875979228363790
      11 -29.8531751466299120 -9.7299979046697445
      12 15.1220023457727230 25.1038379563044960
      13 22.7366718877015930 7.4296786502390768
      14 -9.3218018365645854 -15.2740433268190130
      15 -13.5313752142074220 -3.8514631848696728
      16 5.3051235132574242 7.9065568461965441
      17 6.4704436057717105 1.3310163345710742
      18 -2.4665892743626681 -3.4254372671567435
      19 -2.6437641575552702 -0.3653047198870711
      20 0.7771314187961326 1.2135904493519249
      21 1.0105989058549130 0.1773693978765500
      22 -0.0982971888072048 -0.3562252726334449
      23 -0.3696484621597229 -0.1215738645617847
      24 -0.0388403557770505 0.0924829508022898
      25 0.1178052153278498 0.0570985338450345
      26 0.0239731345108428 -0.0223176692443156
      27 -0.0296092414452849 -0.0163260360921365
      28 -0.0051932812835911 0.0049336345493305
      29 0.0051959782885691 0.0024740275532692
      30 0.0004467306955688 -0.0007151856250505
      31 -0.0005090795488796 -0.0001540322172991
      32 -0.0002283634635380 -0.0000432110135947
      33 0.0001112835192454 0.0001218740679878
      34 0.0001947778117732 0.0000315676817145
      35 -0.0000777784823648 -0.0000901590239588
      36 -0.0000608432770399 0.0000136958535388
      37 0.0000204062418828 0.0000157506941068
      38 0.0000083269893273 -0.0000020983554747];

    % Table 10 from the text
    table10 = [
      0 0.0000000000000000 1.0000000000000000
      1 0.0000000000000000 0.7530928302472617
      2 0.0000000000000000 -5.6876478972560180
      3 0.0000000000000000 -4.6001625314679311
      4 0.0937500000000000 14.7688208805162220
      5 0.0354462028356808 12.9110364270167750
      6 -0.5794683008061320 -23.2985876342930280
      7 -0.1591806968088001 -22.2529692918405930
      8 1.6842716605855981 25.2949343649485390
      9 0.2590332240370517 26.5642830424431010
      10 -3.0399314551829391 -20.7275046216581170
      11 -0.0645486159987336 -23.5088999441590650
      12 3.7757798299644287 14.1983612865681630
      13 -0.4424917088086144 16.1120087317810350
      14 -3.3696060827666314 -8.9979341047037718
      15 0.8677197790058497 -8.7887426960210036
      16 2.1860298844040456 5.4292449518028008
      17 -0.8761868169616034 3.8751690379208039
      18 -1.0192554133470433 -2.9222312270789970
      19 0.5698017880599795 -1.3986551339586100
      20 0.3294504504024299 1.2921105899769376
      21 -0.2542214172218537 0.4215695555772232
      22 -0.0689937284425019 -0.4481976016305430
      23 0.0806987107131334 -0.1083283065391708
      24 0.0097549766071051 0.1241851919484450
      25 -0.0201417209825254 0.0237060050706818
      26 -0.0030113344253303 -0.0323930735345231
      27 0.0057346579076300 -0.0041438522586214
      28 0.0016237250709576 0.0102787025544772
      29 -0.0024134043883321 0.0000794121908992
      30 -0.0002857520641057 -0.0034959095257686
      31 0.0009078237305454 0.0004233705581285
      32 -0.0001481233143099 0.0009223046179774
      33 -0.0002045995112581 -0.0002180191880766
      34 0.0001062038845836 -0.0001333546201417
      35 0.0000104657146092 0.0000637844613933
      36 -0.0000268057642468 -0.0000096296739955
      37 0.0000075227992668 -0.0000116959858281
      38 -0.0000007108447371 -0.0000022007617119];

    % Table 12 from the text
    table12 = [
      0 1.0000000000000000 1.0000000000000000
      1 -7.5838964357215906 -7.5838964357215906
      2 -4.5207915154522071 -4.5207915154522071
      3 67.1126014324299830 66.8001014324299830
      4 36.4257465110401740 38.7957141472031710
      5 -227.5434101958574600 -226.1306628472786400
      6 -148.8868387007147800 -169.9181203983491500
      7 425.6988546482454000 414.7601777953258900
      8 322.9424771530482900 394.3029642171122800
      9 -511.7610421660615400 -469.1169610571635100
      10 -433.3830394473827700 -568.1911141292858900
      11 435.5861115275547500 347.3528087832515900
      12 393.9363791481090400 559.5050863555197800
      13 -289.7726603465931700 -177.4883257131791400
      14 -257.2538676239515800 -401.8573592578553100
      15 169.8215795297231200 74.3740852975847130
      16 127.1823575386862300 223.4938946349586800
      17 -94.8269310329440460 -37.7348164935060740
      18 -50.9124722094202300 -103.2373672310176600
      19 48.3584924238158750 22.7665029241817970
      20 18.1172702688163230 42.6727488787292800
      21 -20.2524941264342860 -10.3674675169848310
      22 -6.1301321875823511 -16.2280741114435580
      23 6.4841086908589745 2.4385277545819979
      24 1.8629524433630982 5.2925499604112982
      25 -1.6546410482740759 0.0727424296455357
      26 -0.4418054918716946 -1.2955329800216224
      27 0.4474639134552081 -0.1458737659618503
      28 0.0855305535945268 0.2286349421460986
      29 -0.1492902692086656 -0.0234720672648906
      30 -0.0198209692030472 -0.0497699599389199
      31 0.0387758704005577 0.0277057758052075
      32 -0.0004298844398339 0.0127677985914393
      33 -0.0066929563118655 -0.0062955861393704
      34 0.0026047337056495 -0.0007922566399001
      35 0.0008666597716715 0.0007258586021021
      36 -0.0005406029268510 0.0000013255383408
      37 -0.0000112860358089 0.0000755918974562];
  end

  methods
    function drag = PadeSphere(varargin)
      % Construct a new creeping flow sphere-wall drag tensor.
      %
      % Usage:
      %   drag = PadeSphere(radius, separation, viscosity, ...)
      %   radius and separation should be specified in the same units.
      %
      % Parameters:
      %   - radius     -- Radius of sphere (default: 1.0)
      %   - separation -- Separation between sphere centre and surface
      %   - viscosity  -- Viscosity of medium (default: 1.0)
      %
      % Parameters can also be passed as named arguments.
      % Unmatched parameters are passed to :class:`Stokes`.

      % Only need a constructor for help/doc functionality
      drag = drag@ott.drag.StokesSphereWall(varargin{:});
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)

      % Calculate stokes sphere drag
      D = ott.drag.StokesSphere(drag.radius, drag.viscosity).forward;

      % l/a: ratio of the distance of the sphere centre from the
      % wall to the sphere radius
      as = drag.radius ./ drag.separation;

      % Check epsilon range
      % This threshold is determines from the text in Chaoui
      if (drag.separation ./ drag.radius - 1) < 0.1
        warning('ott:drag:PadeSphere:small_epsilon', ...
          ['Apprxomation may be inaccurate for small separation', ...
          newline, 'Consider using Chaoui polynomial approximation']);
      end

      % Properties calculated using Equation 5.1
      fxx = sum(drag.table9(:, 2).*as.^drag.table9(:, 1)) ...
          ./ sum(drag.table9(:, 3).*as.^drag.table9(:, 1));
      cyx = sum(drag.table10(:, 2).*as.^drag.table10(:, 1)) ...
          ./ sum(drag.table10(:, 3).*as.^drag.table10(:, 1));
      cyy = sum(drag.table12(:, 2).*as.^drag.table12(:, 1)) ...
          ./ sum(drag.table12(:, 3).*as.^drag.table12(:, 1));

      % Modify diagonal terms (cyy, fxx)
      D(1, 1) = D(1, 1) .* fxx;
      D(2, 2) = D(2, 2) .* fxx;
      D(4, 4) = D(4, 4) .* cyy;
      D(5, 5) = D(5, 5) .* cyy;

      % Add cross-terms (cyx)
      D(1, 5) = -6*pi*drag.viscosity*drag.radius^2*cyx;
      D(2, 4) = 6*pi*drag.viscosity*drag.radius^2*cyx;
      D(5, 1) = -8*pi*drag.viscosity*drag.radius^2*cyx;
      D(4, 2) = 8*pi*drag.viscosity*drag.radius^2*cyx;

      % Add warning about using Faxen for D(3, 3) and D(6, 6)
      warning('ott:drag:ChaouiSphere:faxen_perp_terms', ...
          'Using Faxen corrections for D(3, 3) and D(6, 6)');
      gammaS = 1./(1 - (9/8).*as + (1/2)*as^3);
      betaP = 1./(1 - (1/8)*as^3);
      D(3, 3) = D(3, 3)*gammaS;
      D(6, 6) = D(6, 6)*betaP;
    end
  end
end

