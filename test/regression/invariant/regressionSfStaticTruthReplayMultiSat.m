function regressionSfStaticTruthReplayMultiSat(varargin)
% Regression check for SF static multi-satellite truth replay.
% This script focuses on one narrow contract only:
%   1) the multi-sat static truth point must replay both satellites with a
%      consistent per-satellite objective decomposition;
%   2) the static reference-Doppler semantics must satisfy
%         fdSat = fdRef + deltaFd
%      at the truth point.

opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fixture = buildSfStaticRegressionFixture();

  fprintf('Running regressionSfStaticTruthReplayMultiSat ...\n');

truthOptVar = [fixture.truth.latlonTrueDeg(:); fixture.truth.fdRefTrueHz];
baseOpt = struct();
baseOpt.useLogObjective = false;

modelAllOpt = baseOpt;
modelAllOpt.satWeight = [1; 1];
modelRefOpt = baseOpt;
modelRefOpt.satWeight = [1; 0];
modelOtherOpt = baseOpt;
modelOtherOpt.satWeight = [0; 1];

[modelAll, ~, ~] = buildDoaDopplerSfModel( ...
  fixture.scene, fixture.rxSig, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, 1, modelAllOpt);
[modelRef, ~, ~] = buildDoaDopplerSfModel( ...
  fixture.scene, fixture.rxSig, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, 1, modelRefOpt);
[modelOther, ~, ~] = buildDoaDopplerSfModel( ...
  fixture.scene, fixture.rxSig, fixture.pilotWave, fixture.carrierFreq, ...
  fixture.sampleRate, fixture.viewMs.doaGrid, fixture.fdRange, 1, modelOtherOpt);

[objAll, ~, ~, auxAll] = evalDoaDopplerSfProfileLike(modelAll, truthOptVar);
[objRef, ~, ~, auxRef] = evalDoaDopplerSfProfileLike(modelRef, truthOptVar);
[objOther, ~, ~, auxOther] = evalDoaDopplerSfProfileLike(modelOther, truthOptVar);

if ~all(isfinite([objAll, objRef, objOther]))
  error('regressionSfStaticTruthReplayMultiSat:NonFiniteObjective', ...
    'The SF static truth replay objective must stay finite.');
end
if any(~isfinite(auxAll.fd(:))) || any(~isfinite(auxAll.deltaFd(:)))
  error('regressionSfStaticTruthReplayMultiSat:NonFiniteDopplerState', ...
    'The SF static truth replay Doppler state must stay finite.');
end

objAddErr = abs(objAll - (objRef + objOther));
objSatRefErr = abs(auxAll.objectiveSat(1) - objRef);
objSatOtherErr = abs(auxAll.objectiveSat(2) - objOther);
fdComposeErr = max(abs(auxAll.fd(:) - (fixture.truth.fdRefTrueHz + auxAll.deltaFd(:))));
fdTruthErr = max(abs(auxAll.fd(:) - fixture.truth.fdSatTrueHz(:)));

if objAddErr > 1e-8 * max(abs(objAll), 1)
  error('regressionSfStaticTruthReplayMultiSat:ObjectiveAdditivityMismatch', ...
    'The multi-sat truth objective must equal the sum of its per-satellite slices.');
end
if objSatRefErr > 1e-8 * max(abs(objRef), 1)
  error('regressionSfStaticTruthReplayMultiSat:RefSliceMismatch', ...
    'The ref-satellite truth objective slice is inconsistent.');
end
if objSatOtherErr > 1e-8 * max(abs(objOther), 1)
  error('regressionSfStaticTruthReplayMultiSat:OtherSliceMismatch', ...
    'The non-reference-satellite truth objective slice is inconsistent.');
end
if fdComposeErr > 1e-9
  error('regressionSfStaticTruthReplayMultiSat:FdCompositionMismatch', ...
    'The SF static truth replay must satisfy fdSat = fdRef + deltaFd.');
end
if fdTruthErr > 1e-3
  error('regressionSfStaticTruthReplayMultiSat:TruthFdMismatch', ...
    'The SF static truth replay must recover the exact geometric per-satellite fd values.');
end
if abs(auxAll.deltaFd(fixture.truth.refSatIdxLocal)) > 1e-9
  error('regressionSfStaticTruthReplayMultiSat:ReferenceDeltaMismatch', ...
    'The reference satellite must satisfy deltaFdRef(refSat)=0 at the truth point.');
end

  fprintf('  total truth objective    : %.6f\n', objAll);
  fprintf('  per-sat objective slices : [%.6f, %.6f]\n', auxAll.objectiveSat(1), auxAll.objectiveSat(2));
  fprintf('  additivity error         : %.6g\n', objAddErr);
  fprintf('  fd composition error (Hz): %.6g\n', fdComposeErr);
  fprintf('  fd truth error max (Hz)  : %.6g\n', fdTruthErr);
  fprintf('PASS: regressionSfStaticTruthReplayMultiSat\n');

end
