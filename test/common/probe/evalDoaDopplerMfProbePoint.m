function [probeEval, debugAux, probePoint] = evalDoaDopplerMfProbePoint(model, pointCur, probeOpt)
%EVALDOADOPPLERMFPROBEPOINT Evaluate one MF debug probe point.
% This helper centralizes the one-point probe path used by the development
% landscape / comb-scan scripts. The point is injected through
% model.debugTruth and evaluated by buildDoaDopplerMfDebug so that all probe
% scripts share the same compact evaluator glue.

if nargin < 3 || isempty(probeOpt)
  probeOpt = struct();
end

probeTemplate = localGetStructField(probeOpt, 'probeTemplate', struct());
debugEnable = logical(localGetStructField(probeOpt, 'debugEnable', false));

probePoint = localBuildProbePoint(model, pointCur, probeTemplate);

modelScan = model;
modelScan.debugEnable = debugEnable;
modelScan.debugTruth = struct();
modelScan.debugTruth.doaParam = probePoint.doaParam(:);
modelScan.debugTruth.fdRef = probePoint.fdRef;
modelScan.debugTruth.fdRate = probePoint.fdRate;

debugAux = buildDoaDopplerMfDebug(modelScan, struct(), probePoint, probePoint, struct());
probeEval = localGetStructField(debugAux, 'probe', struct());
probeEval = localGetStructField(probeEval, 'truth', struct());
end

function probePoint = localBuildProbePoint(model, pointCur, probeTemplate)
%LOCALBUILDPROBEPOINT Build one standardized probe-point struct.

probePoint = probeTemplate;
probePoint.doaParam = localResolveDoaParam(model, pointCur);
probePoint.fdRef = localResolveScalarField(pointCur, probeTemplate, 'fdRef', NaN);
probePoint.fdRate = localResolveFdRate(model, pointCur, probeTemplate);
probePoint.obj = localResolveScalarField(probeTemplate, struct(), 'obj', NaN);
probePoint.residualNorm = localResolveScalarField(probeTemplate, struct(), 'residualNorm', NaN);
end

function doaParam = localResolveDoaParam(model, pointCur)
%LOCALRESOLVEDOAPARAM Resolve one DoA parameter vector for the probe point.

if isfield(pointCur, 'doaParam') && ~isempty(pointCur.doaParam)
  doaParam = pointCur.doaParam(:);
  return;
end

switch string(localGetStructField(model, 'doaType', 'latlon'))
  case "latlon"
    if isfield(pointCur, 'latlon') && ~isempty(pointCur.latlon)
      doaParam = pointCur.latlon(:);
      return;
    end
  otherwise
    if isfield(pointCur, 'angle') && ~isempty(pointCur.angle)
      doaParam = pointCur.angle(:);
      return;
    end
end

error('evalDoaDopplerMfProbePoint:MissingDoaParam', ...
  'pointCur must provide doaParam or a mode-compatible field (latlon/angle).');
end

function fdRate = localResolveFdRate(model, pointCur, probeTemplate)
%LOCALRESOLVEFDRATE Resolve one fdRate value consistent with the mode.

switch string(localGetStructField(model, 'fdRateMode', 'unknown'))
  case "known"
    fdRate = localGetStructField(model, 'fdRateKnown', NaN);
  case "zero"
    fdRate = 0;
  otherwise
    fdRate = localResolveScalarField(pointCur, probeTemplate, 'fdRate', NaN);
end
end

function val = localResolveScalarField(primaryStruct, fallbackStruct, fieldName, defaultVal)
%LOCALRESOLVESCALARFIELD Resolve one scalar field from two candidate structs.

val = localGetStructField(primaryStruct, fieldName, []);
if isempty(val)
  val = localGetStructField(fallbackStruct, fieldName, defaultVal);
end
if isempty(val)
  val = defaultVal;
end
end

function val = localGetStructField(s, fieldName, defaultVal)
%LOCALGETSTRUCTFIELD Get one struct field with a default fallback.

if nargin < 3
  defaultVal = [];
end

val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
