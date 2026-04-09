function caseInfo = runDoaOnlyCase(displayName, satMode, view, wavelen, ...
  pilotWave, verbose, modelOpt)
%RUNDOAONLYCASE Run one DoA-only estimator case using a unified wrapper.

[estRaw, ~, ~] = estimatorDoaMlePilotOpt( ...
  view.sceneRef.array, wavelen, view.rxSigSf, pilotWave, ...
  view.doaGrid, [], verbose, modelOpt);

estResult = localWrapDoaOnlyResult(estRaw, view);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa", "none", estResult);
end

function estOut = localWrapDoaOnlyResult(estIn, view)
%LOCALWRAPDOAONLYRESULT Convert DoA-only output to the common summary format.

estOut = estIn;
estOut.modelType = 'doa-only';
estOut.fdRefEst = NaN;
estOut.fdRateEst = NaN;
if ~isfield(estOut, 'aux') || ~isstruct(estOut.aux)
  estOut.aux = struct();
end
if isfield(estIn, 'latlonEst') && ~isfield(estOut.aux, 'latlonEst')
  estOut.aux.latlonEst = estIn.latlonEst;
end

timeOffsetSec = 0;
if isfield(view, 'sceneSeq') && ~isempty(view.sceneSeq)
  timeOffsetSec = view.sceneSeq.timeOffsetSec(view.sceneSeq.refFrameIdx);
end
estOut.timeOffsetSec = timeOffsetSec;
end
