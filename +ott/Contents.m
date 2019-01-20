% ott optical tweezers toolbox
%
% Scattering (T-Matrix) objects
%   Tmatrix                   - Tmatrix abstract class representing T-matrix of a scattering particle
%   TmatrixEbcm               - TmatrixEbcm constructs a T-matrix using extended boundary conditions method
%   TmatrixMie                - TmatrixMie construct T-matrix from Mie scattering coefficients
%   TmatrixPm                 - TmatrixPm constructs a T-matrix using the point matching method
%
% Beam (Bsc) objects
%   Bsc                       - Bsc abstract class representing beam shape coefficients
%   BscBessel                 - BscBessel representation of a bessel beam and bessel-like beams with OAM
%   BscPlane                  - BscPlane representation of a plane wave in VSWF coefficients
%   BscPmGauss                - BscPmGauss provides HG, LG and IG beams using point matching method
%   BscPmParaxial             - BscPmParaxial calculate representation from farfield/paraxial beam
%   BscPointmatch             - BscPointmatch base class for BSC generated using point matching
%
% Functions operating on beams/particles
%   axial_equilibrium         - find equilibrium position and stiffness along beam axis
%   find_equilibrium          - estimates equilibrium positions from position-force data
%   forcetorque               - calculate force/torque/spin in a 3D orthogonal space
%   trap_stiffness            - calculate the force and torque trap stiffness
%
% Other files
%   change_warnings           - enables or disables move/depreciation warnings
%   warning                   - overload of MATLAB warning function for OTT
