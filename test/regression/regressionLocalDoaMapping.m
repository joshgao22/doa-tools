% Regression check for global-to-local DoA mapping.
% This script focuses on one narrow contract only:
%   1) truth localDoa stored in sceneSeq must match direct recomputation by
%      globalToLocalDoa from the frame user position, satellite position,
%      and satellite local rotation;
%   2) the implied local unit direction must rotate back to the original
%      global line-of-sight direction.
clear(); close all; clc;

localAddProjectPath();
fixture = localBuildRegressionFixture();
sceneSeq = fixture.sceneSeq;
localDoaTruth = localExtractTruthLocalDoa(sceneSeq);

fprintf('Running regressionLocalDoaMapping ...\n');

maxAbsAzErrDeg = 0;
maxAbsElErrDeg = 0;
maxDirClosureErr = 0;

for iFrame = 1:sceneSeq.numFrame
  userPosEci = reshape(sceneSeq.usrPosEci(:, 1, iFrame), 3, 1);
  satPosMat = sceneSeq.satPosEci(:, :, iFrame);

  for iSat = 1:sceneSeq.numSat
    rotMat = sceneSeq.rotMat{iFrame, iSat};
    doaDirect = globalToLocalDoa(userPosEci, satPosMat(:, iSat), rotMat);
    doaScene = localDoaTruth(:, iSat, iFrame);

    azErrDeg = abs(rad2deg(localWrapToPi(doaDirect(1) - doaScene(1))));
    elErrDeg = abs(rad2deg(doaDirect(2) - doaScene(2)));
    maxAbsAzErrDeg = max(maxAbsAzErrDeg, azErrDeg);
    maxAbsElErrDeg = max(maxAbsElErrDeg, elErrDeg);

    losGlobal = userPosEci - satPosMat(:, iSat);
    losGlobal = losGlobal / norm(losGlobal);
    losLocal = doa2dir(doaDirect);
    losGlobalBack = rotMat * losLocal;
    dirClosureErr = norm(losGlobalBack - losGlobal);
    maxDirClosureErr = max(maxDirClosureErr, dirClosureErr);
  end
end

if isfield(sceneSeq, 'localDoaRef') && ~isempty(sceneSeq.localDoaRef)
  localDoaRefTruth = localDoaTruth(:, :, sceneSeq.refFrameIdx);
  refDiff = localDoaRefTruth - reshape(sceneSeq.localDoaRef, 2, sceneSeq.numSat);
  refAzErrDeg = max(abs(rad2deg(localWrapToPi(refDiff(1, :)))));
  refElErrDeg = max(abs(rad2deg(refDiff(2, :))));
else
  refAzErrDeg = NaN;
  refElErrDeg = NaN;
end

tolDoaDeg = 1e-9;
tolDirClosure = 1e-12;
if maxAbsAzErrDeg > tolDoaDeg || maxAbsElErrDeg > tolDoaDeg
  error('regressionLocalDoaMapping:TruthMismatch', ...
    ['Direct globalToLocalDoa recomputation does not match sceneSeq.localDoa ', ...
     '(max az err %.3e deg, max el err %.3e deg).'], ...
    maxAbsAzErrDeg, maxAbsElErrDeg);
end
if maxDirClosureErr > tolDirClosure
  error('regressionLocalDoaMapping:DirectionClosureMismatch', ...
    'The local/global direction closure error %.3e exceeds tol %.3e.', ...
    maxDirClosureErr, tolDirClosure);
end
if isfinite(refAzErrDeg) && (refAzErrDeg > tolDoaDeg || refElErrDeg > tolDoaDeg)
  error('regressionLocalDoaMapping:ReferenceSliceMismatch', ...
    ['sceneSeq.localDoaRef is inconsistent with the ref-frame slice of ', ...
     'sceneSeq.localDoa (az %.3e deg, el %.3e deg).'], ...
    refAzErrDeg, refElErrDeg);
end

fprintf('  max azimuth error (deg)  : %.3e\n', maxAbsAzErrDeg);
fprintf('  max elevation error (deg): %.3e\n', maxAbsElErrDeg);
if isfinite(refAzErrDeg)
  fprintf('  ref-frame localDoaRef err: [%.3e, %.3e] deg\n', refAzErrDeg, refElErrDeg);
end
fprintf('  max direction closure err: %.3e\n', maxDirClosureErr);
fprintf('PASS: regressionLocalDoaMapping\n');


function localDoa = localExtractTruthLocalDoa(sceneSeq)
%LOCALEXTRACTTRUTHLOCALDOA Extract 2 x Ns x Nf local DoA truth.

localDoaRaw = sceneSeq.localDoa;
if ndims(localDoaRaw) == 3
  localDoa = localDoaRaw;
  return;
end
if ndims(localDoaRaw) == 4
  localDoa = localDoaRaw(:, :, 1, :);
  return;
end

error('regressionLocalDoaMapping:InvalidLocalDoaSize', ...
  'sceneSeq.localDoa must have size 2xNsxNf or 2xNsxNuxNf.');
end


function wrappedVal = localWrapToPi(val)
%LOCALWRAPTOPI Wrap angles in radians to [-pi, pi).

wrappedVal = mod(val + pi, 2 * pi) - pi;
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact dynamic scene sequence.

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
sceneSeq = genMultiFrameScene(utcVec, tle, usrLla, [1, 2], [], arrUpa, ...
  15, 55, "satellite", 1, refFrameIdx);

fixture = struct();
fixture.sceneSeq = sceneSeq;
end


function localAddProjectPath()
%LOCALADDPROJECTPATH Add the repository folders to the MATLAB path.

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(projectRoot));
end
