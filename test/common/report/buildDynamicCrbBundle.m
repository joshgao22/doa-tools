function crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar)
%BUILDDYNAMICCRBBUNDLE Build one compact CRB bundle for a periodic fixture.
% The dev/perf scripts share the same CRB anchor construction so the entry
% points can stay focused on orchestration and display.

arguments
  periodicFixture (1,1) struct
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  noiseVar (1,1) double
end

truth = periodicFixture.truth;
sceneRef = periodicFixture.sceneRef;
sceneRefOnly = periodicFixture.sceneSeqRefOnly.sceneCell{periodicFixture.sceneSeqRefOnly.refFrameIdx};

crbSfOpt = struct();
crbSfOpt.doaType = 'latlon';
[crbSfRef, auxCrbSfRef] = localTryBuildStaticCrb( ...
  sceneRefOnly, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbSfOpt);
[crbSfMs, auxCrbSfMs] = localTryBuildStaticCrb( ...
  sceneRef, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbSfOpt);

pathGainRef = ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
noiseVarRef = noiseVar * ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
pathGainMs = ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);
noiseVarMs = noiseVar * ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);

crbMfKnownOpt = struct();
crbMfKnownOpt.doaType = 'latlon';
crbMfKnownOpt.phaseMode = 'continuous';
crbMfKnownOpt.fdRateMode = 'known';
crbMfKnownOpt.steeringMode = 'framewise';
crbMfKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
[crbMfRefKnown, auxCrbMfRefKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, crbMfKnownOpt);
[crbMfMsKnown, auxCrbMfMsKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, crbMfKnownOpt);

crbMfUnknownOpt = crbMfKnownOpt;
crbMfUnknownOpt.fdRateMode = 'unknown';
[crbMfRefUnknown, auxCrbMfRefUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, crbMfUnknownOpt);
[crbMfMsUnknown, auxCrbMfMsUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, crbMfUnknownOpt);

crbBundle = struct();
crbBundle.truth = truth;
crbBundle.crbSfRef = crbSfRef;
crbBundle.auxCrbSfRef = auxCrbSfRef;
crbBundle.crbSfMs = crbSfMs;
crbBundle.auxCrbSfMs = auxCrbSfMs;
crbBundle.crbMfRefKnown = crbMfRefKnown;
crbBundle.auxCrbMfRefKnown = auxCrbMfRefKnown;
crbBundle.crbMfMsKnown = crbMfMsKnown;
crbBundle.auxCrbMfMsKnown = auxCrbMfMsKnown;
crbBundle.crbMfRefUnknown = crbMfRefUnknown;
crbBundle.auxCrbMfRefUnknown = auxCrbMfRefUnknown;
crbBundle.crbMfMsUnknown = crbMfMsUnknown;
crbBundle.auxCrbMfMsUnknown = auxCrbMfMsUnknown;
end

function [crb, aux] = localTryBuildStaticCrb(scene, pilotWave, carrierFreq, sampleRate, doaParam, fdRefHz, numSource, noiseVar, crbOpt)
try
  [crb, aux] = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, numSource, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = "zero";
  aux.phaseMode = "single-frame";
  aux.errorMessage = string(ME.message);
end
end

function [crb, aux] = localTryBuildDynamicCrb(sceneSeq, pilotWave, carrierFreq, sampleRate, doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt)
try
  [crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = getDoaDopplerFieldOrDefault(crbOpt, 'fdRateMode', "unknown");
  aux.phaseMode = getDoaDopplerFieldOrDefault(crbOpt, 'phaseMode', "continuous");
  aux.errorMessage = string(ME.message);
end
end
