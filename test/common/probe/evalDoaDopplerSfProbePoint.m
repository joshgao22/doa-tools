function probe = evalDoaDopplerSfProbePoint(view, pilotWave, carrierFreq, sampleRate, ...
  fdRange, modelOpt, probePoint)
%EVALDOADOPPLERSFPROBEPOINT Evaluate one fixed SF static probe point.
% This helper centralizes the fixed-point single-frame static evaluator path
% used by the MS-vs-SS regression scripts. It keeps the regressions focused
% on contracts instead of re-implementing model build + profile evaluation
% glue in multiple files.

arguments
  view (1, 1) struct
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  fdRange (1, 2) double
  modelOpt (1, 1) struct = struct()
  probePoint (1, 1) struct = struct()
end

numSource = 1;
[model, ~, ~, ~] = buildDoaDopplerSfModel( ...
  view.sceneRef, view.rxSigSf, pilotWave, carrierFreq, sampleRate, ...
  view.doaGrid, fdRange, numSource, modelOpt);
optVar = localBuildOptVar(model, probePoint);
[obj, pathGain, noiseVar, aux] = evalDoaDopplerSfProfileLike(model, optVar);

probe = struct();
probe.model = model;
probe.optVar = optVar;
probe.obj = obj;
probe.pathGain = pathGain;
probe.noiseVar = noiseVar;
probe.aux = aux;
end


function optVar = localBuildOptVar(model, probePoint)
%LOCALBUILDOPTVAR Build one [DoA; fdRef] vector for the SF evaluator.

doaParam = localResolveDoaParam(model, probePoint);
fdRef = localResolveFdRef(probePoint);
optVar = [doaParam(:); fdRef(:)];
end

function doaParam = localResolveDoaParam(model, probePoint)
%LOCALRESOLVEDOAPARAM Resolve one DoA vector compatible with the model.

if isfield(probePoint, 'doaParam') && ~isempty(probePoint.doaParam)
  doaParam = reshape(probePoint.doaParam, [], 1);
  return;
end

switch string(localGetStructField(model, 'doaType', 'latlon'))
  case "latlon"
    if isfield(probePoint, 'latlon') && ~isempty(probePoint.latlon)
      doaParam = reshape(probePoint.latlon, [], 1);
      return;
    end
  otherwise
    if isfield(probePoint, 'angle') && ~isempty(probePoint.angle)
      doaParam = reshape(probePoint.angle, [], 1);
      return;
    end
end

error('evalDoaDopplerSfProbePoint:MissingDoaParam', ...
  'probePoint must provide doaParam or a mode-compatible field (latlon/angle).');
end

function fdRef = localResolveFdRef(probePoint)
%LOCALRESOLVEFDREF Resolve one fdRef scalar.

fdRef = localGetStructField(probePoint, 'fdRef', []);
if isempty(fdRef)
  error('evalDoaDopplerSfProbePoint:MissingFdRef', ...
    'probePoint.fdRef must be provided for SF static probing.');
end
fdRef = reshape(fdRef, [], 1);
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
