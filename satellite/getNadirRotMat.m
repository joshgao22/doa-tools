function rotMat = getNadirRotMat(satPosEci, satVelEci)
%GETNADIRROTMAT Compute local-to-ECI rotation matrix for nadir pointing.
% Constructs a right-handed local coordinate frame assuming the satellite
% array boresight points toward the Earth center (nadir).
%
%Syntax:
%  rotMat = getNadirRotMat(satPosEci, satVelEci)
%
%Inputs:
%  satPosEci   : 3x1 satellite position vector in ECI coordinates
%  satVelEci   : 3x1 satellite velocity vector in ECI coordinates
%
%Output:
%  rotMat : 3x3 rotation matrix mapping local frame to ECI frame
%
%Local frame definition:
%  zLocal : nadir direction (pointing to Earth center)
%  yLocal : orbit normal direction
%  xLocal : along-track direction
%
%Construction:
%  zLocal = -satPosEci / ||satPosEci||
%  yLocal = (satPosEci × satVelEci) / ||satPosEci × satVelEci||
%  xLocal = yLocal × zLocal
%
%The rotation matrix is
%
%  rotMat = [xLocal, yLocal, zLocal]
%
%which maps vectors from the local frame to the ECI frame.
%
%Example:
%  rotMat = getNadirRotMat(satPosEci, satVelEci);
%
%See also: cross, norm

arguments
  satPosEci (3,1) {mustBeNumeric}
  satVelEci (3,1) {mustBeNumeric}
end

% -------------------------------------------------------------------------
% Nadir direction (z-axis)
% -------------------------------------------------------------------------
zLocal = -satPosEci / norm(satPosEci);

% -------------------------------------------------------------------------
% Orbit normal direction (y-axis)
% -------------------------------------------------------------------------
orbitNormal = cross(satPosEci, satVelEci);
yLocal = orbitNormal / norm(orbitNormal);

% -------------------------------------------------------------------------
% Along-track direction (x-axis)
% -------------------------------------------------------------------------
xLocal = cross(yLocal, zLocal);
xLocal = xLocal / norm(xLocal);

% Re-orthogonalize y-axis for numerical stability
yLocal = cross(zLocal, xLocal);

% -------------------------------------------------------------------------
% Assemble rotation matrix (local → ECI)
% -------------------------------------------------------------------------
rotMat = [xLocal, yLocal, zLocal];

end