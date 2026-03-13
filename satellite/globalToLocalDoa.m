function doasLoc = globalToLocalDoa(posGlob, posLoc, rotLocToGlob)
%GLOBALTOLOCALDOA Convert user positions in global coordinates to local DOAs
%as observed by an array or satellite.
%
%   doasLoc = globalToLocalDoa(posGlob, posLoc, rotLocToGlob)
%
%Inputs:
%   posGlob : 3 × K
%       User positions in the global frame (meters). Each column corresponds 
%       to one user.
%
%   posLoc : 3 × 1
%       Array/satellite position in the global frame (meters).
%
%   rotLocToGlob : 3 × 3
%       Local-to-global rotation matrix of the array.
%       Each column is a unit local axis expressed in the global frame.
%       Thus, for a vector uLocal:    uGlobal = rotLocToGlob * uLocal.
%       This function uses its transpose to obtain local coordinates:
%           uLocal = rotLocToGlob' * uGlobal.
%
%Output:
%   GLOBALTOLOCALDOA : 2 × K
%       Local azimuth–elevation DOAs in radians:
%           az = atan2(uY, uX)
%           el = asin(uZ)
%       where [uX; uY; uZ] is the line-of-sight unit vector expressed
%       in the array's local coordinate frame.
%
%Description:
%   For each user position posGlob(:,k), the function computes the line-of-sight
%   direction from the array at posLoc to the user, normalizes it, rotates the
%   resulting vector into the local array frame using rotLocToGlob', and finally
%   converts it to azimuth and elevation angles. The implementation is fully
%   vectorized over K users.
%
%Notes:
%   - The azimuth convention follows atan2(y, x) measured counterclockwise
%     from the local +x axis toward the local +y axis.
%   - The elevation is measured from the local x–y plane toward +z.
%

arguments
  posGlob (3,:) {mustBeNumeric, mustBeReal}
  posLoc (3,1) {mustBeNumeric, mustBeReal}
  rotLocToGlob (3,3) {mustBeNumeric, mustBeReal}
end

% Global line-of-sight vectors: array -> user (un-normalized)
vGlob = posGlob - posLoc;              % 3 × K

% Normalize column-wise
vNorm = vGlob ./ vecnorm(vGlob, 2, 1); % 3 × K

% Transform to local coordinates: uLocal = qL2G' * vGlobal
uLocal = rotLocToGlob' * vNorm;        % 3 × K

uX = uLocal(1, :);
uY = uLocal(2, :);
uZ = uLocal(3, :);

% Azimuth and elevation
azL = mod(atan2(uY, uX), 2*pi);  % 1 × K
elL = asin(uZ);                  % 1 × K

doasLoc = [azL; elL];            % 2 × K
end