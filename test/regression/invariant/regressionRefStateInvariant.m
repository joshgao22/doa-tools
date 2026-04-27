function regressionRefStateInvariant(varargin)
% Regression check for reference-satellite Doppler state invariants.
% This script focuses on one narrow contract only:
%   1) fdSat(refSat,:) = fdRef
%   2) deltaFd(refSat,:) = 0
%   3) fdSat = fdRef + deltaFd
% and verifies that the same contract still holds after selecting a
% single-satellite subset scene sequence.

opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = localBuildRegressionFixture();
sceneSeq = fixture.sceneSeq;
sceneRef = fixture.sceneRef;

  fprintf('Running regressionRefStateInvariant ...\n');

fullState = buildReferenceDopplerState(sceneRef, ...
  sceneSeq.satPosEci, sceneSeq.satVelEci, ...
  sceneSeq.usrPosEci, sceneSeq.usrVelEci, ...
  fixture.wavelen, sceneSeq.timeOffsetSec);
localCheckReferenceInvariant(fullState, 'fullScene');

sceneSeqRefOnly = selectSatSceneSeq(sceneSeq, fullState.refSatIdxLocal);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
refOnlyState = buildReferenceDopplerState(sceneRefOnly, ...
  sceneSeqRefOnly.satPosEci, sceneSeqRefOnly.satVelEci, ...
  sceneSeqRefOnly.usrPosEci, sceneSeqRefOnly.usrVelEci, ...
  fixture.wavelen, sceneSeqRefOnly.timeOffsetSec);
localCheckReferenceInvariant(refOnlyState, 'refOnlySubset');

if refOnlyState.refSatIdxLocal ~= 1
  error('regressionRefStateInvariant:InvalidSubsetReference', ...
    'The selected single-satellite subset must keep refSatIdxLocal = 1.');
end

  fprintf('  fullScene refSatIdxLocal      : %d\n', fullState.refSatIdxLocal);
  fprintf('  fullScene refSatIdxGlobal     : %d\n', fullState.refSatIdxGlobal);
  fprintf('  subset refSatIdxLocal         : %d\n', refOnlyState.refSatIdxLocal);
  fprintf('  subset fdRef at ref frame (Hz): %.6f\n', refOnlyState.fdRefRefFrame);
  fprintf('PASS: regressionRefStateInvariant\n');


end

function localCheckReferenceInvariant(dopplerState, caseTag)
%LOCALCHECKREFERENCEINVARIANT Check the shared reference-state identities.

refSatIdx = dopplerState.refSatIdxLocal;
fdRefMat = repmat(dopplerState.fdRef, size(dopplerState.fdSat, 1), 1);
tol = 1e-9 * max(1, max(abs(dopplerState.fdSat(:))));

deltaRefMax = max(abs(dopplerState.deltaFd(refSatIdx, :)));
fdRefMatchMax = max(abs(dopplerState.fdSat(refSatIdx, :) - dopplerState.fdRef));
fdComposeMax = max(abs(dopplerState.fdSat(:) - (fdRefMat(:) + dopplerState.deltaFd(:))));

if deltaRefMax > tol
  error('regressionRefStateInvariant:InvalidDeltaFd', ...
    '[%s] max abs deltaFd(refSat,:) = %.3e exceeds tol %.3e.', ...
    caseTag, deltaRefMax, tol);
end
if fdRefMatchMax > tol
  error('regressionRefStateInvariant:InvalidFdRef', ...
    '[%s] max abs fdSat(refSat,:) - fdRef = %.3e exceeds tol %.3e.', ...
    caseTag, fdRefMatchMax, tol);
end
if fdComposeMax > tol
  error('regressionRefStateInvariant:InvalidFdComposition', ...
    '[%s] max abs(fdSat - fdRef - deltaFd) = %.3e exceeds tol %.3e.', ...
    caseTag, fdComposeMax, tol);
end
end


function fixture = localBuildRegressionFixture()
%LOCALBUILDREGRESSIONFIXTURE Build one compact dynamic two-satellite scene.

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
fixture.sceneRef = sceneSeq.sceneCell{sceneSeq.refFrameIdx};
fixture.wavelen = wavelen;
end
