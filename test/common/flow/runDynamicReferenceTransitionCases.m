function [caseDynRefKnown, caseDynRefUnknown, refRunMeta] = runDynamicReferenceTransitionCases(periodicFixture, ...
  pilotWave, carrierFreq, sampleRate, optVerbose, flowOpt, initParamStaticRef, caseStaticRefOnly)
%RUNDYNAMICREFERENCETRANSITIONCASES Run SS-MF CP-K/CP-U transition cases.

refRunMeta = struct();
caseTimer = tic;
caseDynRefKnown = localRunRefKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, initParamStaticRef);
refRunMeta.refKnownSec = toc(caseTimer);
caseTimer = tic;
caseDynRefUnknown = localRunRefUnknownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynRefKnown, caseStaticRefOnly);
refRunMeta.refUnknownSec = toc(caseTimer);
end

function caseDynRefKnown = localRunRefKnownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, initParamStaticRef)
%LOCALRUNREFKNOWNCASE Run one SS-MF-CP-K case.

dynRefKnownOpt = flowOpt.dynBaseOpt;
dynRefKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynRefKnownOpt.initDoaHalfWidth = flowOpt.refKnownDoaHalfWidth;
dynRefKnownOpt.initDoaParam = reshape(initParamStaticRef(1:2), [], 1);

caseDynRefKnown = runDynamicDoaDopplerCase("SS-MF-CP-K", "single", ...
  periodicFixture.viewRefOnly, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynRefKnownOpt, true, ...
  periodicFixture.debugTruthRef, initParamStaticRef);
end


function caseDynRefUnknown = localRunRefUnknownCase(periodicFixture, pilotWave, carrierFreq, sampleRate, ...
  optVerbose, flowOpt, caseDynRefKnown, caseStaticRefOnly)
%LOCALRUNREFUNKNOWNCASE Run one SS-MF-CP-U case.

dynRefUnknownOpt = flowOpt.dynBaseOpt;
dynRefUnknownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
dynRefUnknownOpt.initDoaParam = caseStaticRefOnly.estResult.doaParamEst(:);
dynRefUnknownOpt.initDoaHalfWidth = flowOpt.refUnknownDoaHalfWidth;

initParamRefUnknownCpK = buildDynamicInitParamFromCase(caseDynRefKnown, false, caseDynRefKnown.estResult.fdRateEst);
initParamRefUnknownStatic = buildDynamicInitParamFromCase(caseStaticRefOnly, false, periodicFixture.truth.fdRateFit);
refUnknownCand = buildUnknownInitCandidateSet( ...
  caseDynRefKnown, caseStaticRefOnly, initParamRefUnknownCpK, initParamRefUnknownStatic, ...
  dynRefUnknownOpt.initDoaHalfWidth, flowOpt.refKnownDoaHalfWidth);
caseDynRefUnknown = runDynamicDoaDopplerCase("SS-MF-CP-U", "single", ...
  periodicFixture.viewRefOnly, periodicFixture.truth, pilotWave, carrierFreq, sampleRate, ...
  periodicFixture.fdRange, periodicFixture.fdRateRange, optVerbose, dynRefUnknownOpt, false, ...
  periodicFixture.debugTruthRef, refUnknownCand);
end


