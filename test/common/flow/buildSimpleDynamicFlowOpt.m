function flowOpt = buildSimpleDynamicFlowOpt(varargin)
%BUILDSIMPLEDYNAMICFLOWOPT Build the simplified dynamic-flow option bundle.
% The simplified flow keeps only the pieces needed for
%   static seed -> subset free search -> periodic narrow refine
% where the periodic replay first compares narrow frozen same-tooth
% candidates and only then, for clearly trusted same-tooth hard cases,
% runs one very-small DoA polish.
% The flow intentionally excludes wide periodic fallback, known-rate
% anchors, and blanket same-tooth rescue branches.

opt = localParseFlowOpt(varargin{:});

flowOpt = struct();
flowOpt.doaOnlyOpt = localMergeStruct(struct('useLogObjective', true), ...
  localGetFieldOrDefault(opt, 'doaOnlyOpt', struct()));
flowOpt.staticBaseOpt = localMergeStruct(struct('useLogObjective', true), ...
  localGetFieldOrDefault(opt, 'staticBaseOpt', struct()));
flowOpt.staticMsHalfWidth = reshape(localGetFieldOrDefault(opt, 'staticMsHalfWidth', [0.002; 0.002]), [], 1);

flowOpt.dynBaseOpt = localMergeStruct(struct( ...
  'useLogObjective', true, ...
  'initFdCount', 81, ...
  'useAccessMask', false, ...
  'phaseMode', 'continuous', ...
  'steeringMode', 'framewise', ...
  'continuousPhaseConsistencyWeight', 0.05, ...
  'continuousPhaseCollapsePenaltyWeight', 0.10, ...
  'continuousPhaseNegativeProjectionPenaltyWeight', 0.10, ...
  'unknownWarmAnchorUseScaledSolve', true, ...
  'unknownWarmAnchorFallbackSqp', true, ...
  'debugEnable', true, ...
  'debugStoreEvalTrace', false, ...
  'debugMaxEvalTrace', 120), ...
  localGetFieldOrDefault(opt, 'dynBaseOpt', struct()));

flowOpt.singleDoaHalfWidth = reshape(localGetFieldOrDefault(opt, 'singleDoaHalfWidth', [0.003; 0.003]), [], 1);
flowOpt.multiDoaHalfWidth = reshape(localGetFieldOrDefault(opt, 'multiDoaHalfWidth', [0.002; 0.002]), [], 1);
flowOpt.subsetDoaHalfWidthDeg = reshape(localGetFieldOrDefault(opt, 'subsetDoaHalfWidthDeg', [0.01; 0.01]), [], 1);
flowOpt.periodicRefineFdHalfWidthHz = localGetFieldOrDefault(opt, 'periodicRefineFdHalfWidthHz', 50);
flowOpt.periodicRefineFdRateHalfWidthHzPerSec = localGetFieldOrDefault(opt, 'periodicRefineFdRateHalfWidthHzPerSec', 100);
flowOpt.periodicRefineDoaHalfWidthDeg = reshape(localGetFieldOrDefault(opt, 'periodicRefineDoaHalfWidthDeg', [1e-8; 1e-8]), [], 1);
flowOpt.periodicRefineFreezeDoa = logical(localGetFieldOrDefault(opt, 'periodicRefineFreezeDoa', true));
flowOpt.periodicRefineDoaSeedMode = string(localGetFieldOrDefault(opt, 'periodicRefineDoaSeedMode', "dualWhenMulti"));
flowOpt.periodicRefineMaxSubsetDoaDriftDeg = localGetFieldOrDefault(opt, 'periodicRefineMaxSubsetDoaDriftDeg', 0.003);
flowOpt.periodicRefinePreferStaticWhenSubsetDriftLarge = logical(localGetFieldOrDefault(opt, 'periodicRefinePreferStaticWhenSubsetDriftLarge', true));
flowOpt.periodicRefineSubsetTrustMinRelativeMargin = localGetFieldOrDefault(opt, 'periodicRefineSubsetTrustMinRelativeMargin', 5e-4);
flowOpt.periodicRefineMaxFrozenDoaDisagreementDeg = localGetFieldOrDefault(opt, 'periodicRefineMaxFrozenDoaDisagreementDeg', 0.003);
flowOpt.periodicRefineMaxFrozenRelativeObjGapForTie = localGetFieldOrDefault(opt, 'periodicRefineMaxFrozenRelativeObjGapForTie', 5e-4);
flowOpt.periodicRefineEnableVerySmallDoaPolish = logical(localGetFieldOrDefault(opt, 'periodicRefineEnableVerySmallDoaPolish', true));
flowOpt.periodicRefinePolishDoaHalfWidthDeg = reshape(localGetFieldOrDefault(opt, 'periodicRefinePolishDoaHalfWidthDeg', [0.004; 0.004]), [], 1);
flowOpt.periodicRefinePolishTriggerDoaDriftDeg = localGetFieldOrDefault(opt, 'periodicRefinePolishTriggerDoaDriftDeg', 0.0015);
flowOpt.periodicRefinePolishMaxFrozenRelativeObjGap = localGetFieldOrDefault(opt, 'periodicRefinePolishMaxFrozenRelativeObjGap', 5e-4);
flowOpt.parallelOpt = localMergeStruct(struct('enableSubsetEvalParfor', true, 'minSubsetEvalParfor', 4), ...
  localGetFieldOrDefault(opt, 'parallelOpt', struct()));
end

function opt = localParseFlowOpt(varargin)
opt = struct();
if nargin == 1 && isstruct(varargin{1})
  opt = varargin{1};
elseif mod(nargin, 2) == 0
  for iArg = 1:2:nargin
    name = varargin{iArg};
    if ~(ischar(name) || (isstring(name) && isscalar(name)))
      error('buildSimpleDynamicFlowOpt:InvalidNameValue', ...
        'Name-value inputs must use character or string names.');
    end
    opt.(char(name)) = varargin{iArg + 1};
  end
elseif nargin ~= 0
  error('buildSimpleDynamicFlowOpt:InvalidInput', ...
    'Use either one override struct or name-value pairs.');
end
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end

function outStruct = localMergeStruct(baseStruct, overrideStruct)
outStruct = baseStruct;
if nargin < 2 || isempty(overrideStruct)
  return;
end
fieldList = fieldnames(overrideStruct);
for iField = 1:numel(fieldList)
  outStruct.(fieldList{iField}) = overrideStruct.(fieldList{iField});
end
end
