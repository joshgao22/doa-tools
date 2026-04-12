% Regression check for buildUserStateFromLatlon.
% This script focuses on two narrow contracts only:
%   1) with a motion-perturbed sceneSeq, buildUserStateFromLatlon must
%      reproduce the stored usrPosEci / usrVelEci rather than falling back
%      to Earth-fixed nominal motion only;
%   2) datetime, UTC datevec, and datenum inputs must all yield the same
%      nominal single-epoch state.
clear(); close all;

localAddProjectPath();
fixture = localBuildRegressionFixture();
sceneSeq = fixture.sceneSeq;
latlonTrueDeg = fixture.usrLla(1:2, 1);

fprintf('Running regressionUserStateFromLatlon ...\n');

[userPosSeq, userVelSeq, auxSeq] = buildUserStateFromLatlon(latlonTrueDeg, sceneSeq);
userPosTruth = reshape(sceneSeq.usrPosEci(:, 1, :), 3, []);
userVelTruth = reshape(sceneSeq.usrVelEci(:, 1, :), 3, []);
userPosNominalTruth = reshape(sceneSeq.usrPosNominalEci(:, 1, :), 3, []);
userVelNominalTruth = reshape(sceneSeq.usrVelNominalEci(:, 1, :), 3, []);

posSeqErr = max(vecnorm(userPosSeq - userPosTruth, 2, 1));
velSeqErr = max(vecnorm(userVelSeq - userVelTruth, 2, 1));
nominalGap = max(vecnorm(userPosSeq - userPosNominalTruth, 2, 1));
offsetPosErr = max(vecnorm(auxSeq.userPosOffsetEci - (userPosTruth - userPosNominalTruth), 2, 1));
offsetVelErr = max(vecnorm(auxSeq.userVelOffsetEci - (userVelTruth - userVelNominalTruth), 2, 1));

if posSeqErr > 1e-9 || velSeqErr > 1e-9
  error('regressionUserStateFromLatlon:SceneSeqMismatch', ...
    ['Motion-aware sceneSeq reconstruction failed ', ...
     '(max pos err %.3e m, max vel err %.3e m/s).'], ...
    posSeqErr, velSeqErr);
end
if nominalGap < 1e-3
  error('regressionUserStateFromLatlon:MissingMotionOffset', ...
    'The regression fixture did not create a meaningful user-motion offset.');
end
if offsetPosErr > 1e-9 || offsetVelErr > 1e-9
  error('regressionUserStateFromLatlon:OffsetMismatch', ...
    ['Recovered motion offsets do not match sceneSeq truth ', ...
     '(pos %.3e, vel %.3e).'], offsetPosErr, offsetVelErr);
end
if strcmp(string(auxSeq.offsetSource), "none")
  error('regressionUserStateFromLatlon:OffsetSourceMissing', ...
    'A motion-perturbed sceneSeq must report a non-empty offsetSource.');
end

sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
utcScalar = sceneRef.utc;
[userPosDatetime, userVelDatetime] = buildUserStateFromLatlon(latlonTrueDeg, utcScalar);
[userPosDatevec, userVelDatevec] = buildUserStateFromLatlon(latlonTrueDeg, datevec(utcScalar));
[userPosDatenum, userVelDatenum] = buildUserStateFromLatlon(latlonTrueDeg, datenum(utcScalar));

userPosNominalRef = reshape(sceneSeq.usrPosNominalEci(:, 1, sceneSeq.refFrameIdx), 3, 1);
userVelNominalRef = reshape(sceneSeq.usrVelNominalEci(:, 1, sceneSeq.refFrameIdx), 3, 1);

posDatetimeErr = norm(userPosDatetime - userPosNominalRef);
velDatetimeErr = norm(userVelDatetime - userVelNominalRef);
posDatevecErr = norm(userPosDatevec - userPosNominalRef);
velDatevecErr = norm(userVelDatevec - userVelNominalRef);
posDatenumErr = norm(userPosDatenum - userPosNominalRef);
velDatenumErr = norm(userVelDatenum - userVelNominalRef);

if max([posDatetimeErr, velDatetimeErr, posDatevecErr, velDatevecErr, posDatenumErr, velDatenumErr]) > 1e-9
  error('regressionUserStateFromLatlon:UtcCompatibilityMismatch', ...
    'datetime / datevec / datenum inputs do not reproduce the same nominal state.');
end

fprintf('  max sceneSeq pos err (m)    : %.3e\n', posSeqErr);
fprintf('  max sceneSeq vel err (m/s)  : %.3e\n', velSeqErr);
fprintf('  max nominal gap (m)         : %.3e\n', nominalGap);
fprintf('  offset source               : %s\n', string(auxSeq.offsetSource));
fprintf('  UTC compatibility max err   : %.3e\n', ...
  max([posDatetimeErr, velDatetimeErr, posDatevecErr, velDatevecErr, posDatenumErr, velDatenumErr]));
fprintf('PASS: regressionUserStateFromLatlon\n');


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one motion-aware dynamic scene sequence.

rng(253);

numFrame = 7;
refFrameIdx = ceil(numFrame / 2);
frameIntvlSec = 1 / 750;
carrierFreq = 11.7e9;
lightSpeed = 299792458;
wavelen = lightSpeed / carrierFreq;

utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVec = utcRef + seconds(((1:numFrame) - refFrameIdx) * frameIntvlSec);
usrLla = [37.78; 36.59; 0];

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
arrUpa = createUpa([4, 4], wavelen / 2);

motionOpt = struct();
motionOpt.enableUsrMotion = true;
motionOpt.usrVelEciExtra = [55; -35; 12];
motionOpt.usrAccEciExtra = [0.2; -0.1; 0.05];
motionOpt.motionRefFrameIdx = refFrameIdx;
motionOpt.recomputeAccess = true;

sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1, refFrameIdx, motionOpt);

fixture = struct();
fixture.sceneSeq = sceneSeq;
fixture.usrLla = usrLla;
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
