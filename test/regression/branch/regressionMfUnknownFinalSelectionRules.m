function regressionMfUnknownFinalSelectionRules(varargin)
% Regression checks for final unknown-stage candidate selection rules.
%
% This script keeps the small synthetic final-selection contracts in one
% table-driven guard.  Each case builds only summary structs and calls
% selectFinalDynamicUnknownSummary; no truth-aware selection signal is used.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

  fprintf('Running regressionMfUnknownFinalSelectionRules ...\n');

caseList = { ...
  localObjectiveOnlyCase(); ...
  localSubsetAnchorGuardCase(); ...
  localFdWideWinsCase(); ...
  localWideRefineFromFdWideCase(); ...
  localSameToothWideWinsCase(); ...
  localCloserToZeroToothWinsCase() ...
  };

for iCase = 1:numel(caseList)
  localRunCase(caseList{iCase}, verbose);
end

  fprintf('PASS: regressionMfUnknownFinalSelectionRules\n');



end

function localRunCase(caseSpec, verbose)
if isfield(caseSpec, 'selectOpt') && ~isempty(caseSpec.selectOpt)
  [selectedIdx, selectedTag, reasonText, scoreMat] = selectFinalDynamicUnknownSummary( ...
    caseSpec.summaryList, caseSpec.tagList, caseSpec.selectOpt);
else
  [selectedIdx, selectedTag, reasonText, scoreMat] = selectFinalDynamicUnknownSummary( ...
    caseSpec.summaryList, caseSpec.tagList);
end

if selectedIdx ~= caseSpec.expectedIdx || selectedTag ~= caseSpec.expectedTag
  error('regressionMfUnknownFinalSelectionRules:UnexpectedSelection', ...
    '[%s] expected %s at index %d, but got %s at index %d.', ...
    caseSpec.name, caseSpec.expectedTag, caseSpec.expectedIdx, selectedTag, selectedIdx);
end

if strlength(caseSpec.reasonContains) > 0 && ~contains(reasonText, caseSpec.reasonContains)
  error('regressionMfUnknownFinalSelectionRules:UnexpectedReason', ...
    '[%s] expected selection reason to mention "%s", but got %s.', ...
    caseSpec.name, caseSpec.reasonContains, reasonText);
end

if strlength(caseSpec.reasonNotContains) > 0 && contains(reasonText, caseSpec.reasonNotContains)
  error('regressionMfUnknownFinalSelectionRules:UnexpectedGuardReason', ...
    '[%s] selection reason should not mention "%s": %s.', ...
    caseSpec.name, caseSpec.reasonNotContains, reasonText);
end

  fprintf('  %-34s -> idx=%d, tag=%s, reason=%s\n', ...
    caseSpec.name, selectedIdx, selectedTag, reasonText);
  for iRow = 1:size(scoreMat, 1)
    fprintf('    score %-24s : %s\n', caseSpec.tagList(iRow), localFormatNumericRow(scoreMat(iRow, :)));
  end
end


function caseSpec = localObjectiveOnlyCase()
caseSpec = localEmptyCase("objective-only final ranking");
caseSpec.summaryList = { ...
  localBuildTruthErrorSummary(true,  -1.5082e6, 6.5616e4, 1.0e-3,   1.0e-2,  10.0, 0, 12.0); ...
  localBuildTruthErrorSummary(true,  -1.5086e6, 6.5576e4, 5.0e-2, 500.0, 5000.0, 0, 18.0) ...
  };
caseSpec.tagList = ["periodic-wide"; "periodic-in-tooth"];
caseSpec.expectedIdx = 2;
caseSpec.expectedTag = "periodic-in-tooth";
caseSpec.reasonContains = "final objective";
end


function caseSpec = localSubsetAnchorGuardCase()
caseSpec = localEmptyCase("subset-anchor guard");
caseSpec.summaryList = { ...
  localBuildStateSummary(true, -1.5086e6, 6.5576e4, [37.7900, 36.6000], -29190.0, -3826.0, 1, 12.0); ...
  localBuildStateSummary(true, -1.50855e6, 6.5580e4, [37.7800, 36.5900], -28940.0, -3833.5, 0, 11.0); ...
  localBuildStateSummary(true, -1.5084e6, 6.5590e4, [37.7810, 36.5910], -28939.8, -3832.0, 0, 18.0) ...
  };
caseSpec.tagList = ["periodic-wide"; "periodic-fd-anchor"; "periodic-in-tooth"];
caseSpec.selectOpt = localBuildAnchorGuardOpt(2);
caseSpec.expectedIdx = 2;
caseSpec.expectedTag = "periodic-fd-anchor";
caseSpec.reasonContains = "subset-anchor guard";
end


function caseSpec = localFdWideWinsCase()
caseSpec = localEmptyCase("fd-wide wins");
caseSpec.summaryList = { ...
  localBuildStateSummary(true, -1.5079e6, 6.5660e4, [37.7810, 36.5920], -28939.9, -3833.5, 0, 12.0); ...
  localBuildStateSummary(true, -1.5084e6, 6.5600e4, [37.7795, 36.5905], -28940.0, -3833.6, 0, 11.0); ...
  localBuildStateSummary(true, -1.5081e6, 6.5635e4, [37.7890, 36.5990], -28190.0, -3834.0, 1, 10.0); ...
  localBuildStateSummary(true, -1.5080e6, 6.5640e4, [37.7808, 36.5917], -28939.8, -3834.1, 0, 18.0) ...
  };
caseSpec.tagList = ["periodic-wide"; "periodic-fd-wide"; "periodic-fd-anchor"; "periodic-in-tooth"];
caseSpec.selectOpt = localBuildAnchorGuardOpt(3);
caseSpec.expectedIdx = 2;
caseSpec.expectedTag = "periodic-fd-wide";
caseSpec.reasonNotContains = "subset-anchor guard";
end


function caseSpec = localWideRefineFromFdWideCase()
caseSpec = localEmptyCase("wide refine from fd-wide");
caseSpec.summaryList = { ...
  localBuildStateSummary(true, -1.50845e6, 6.5595e4, [37.7798, 36.5902], -28940.0, -3833.5, 0, 13.0); ...
  localBuildStateSummary(true, -1.50840e6, 6.5600e4, [37.7795, 36.5905], -28940.1, -3833.6, 0, 11.0); ...
  localBuildStateSummary(true, -1.50810e6, 6.5630e4, [37.7890, 36.5990], -28190.0, -3834.0, 1, 10.0); ...
  localBuildStateSummary(true, -1.50835e6, 6.5605e4, [37.7797, 36.5904], -28939.9, -3833.7, 0, 18.0) ...
  };
caseSpec.tagList = ["periodic-wide"; "periodic-fd-wide"; "periodic-fd-anchor"; "periodic-in-tooth"];
caseSpec.selectOpt = localBuildAnchorGuardOpt(3);
caseSpec.expectedIdx = 1;
caseSpec.expectedTag = "periodic-wide";
caseSpec.reasonNotContains = "subset-anchor guard";
end


function caseSpec = localSameToothWideWinsCase()
caseSpec = localEmptyCase("same-tooth wide wins");
caseSpec.summaryList = { ...
  localBuildStateSummary(true, -1.5083e6, 6.5607e4, [37.7790, 36.5910], -28940.01, -3833.5, 0, 12.0); ...
  localBuildStateSummary(true, -1.5080e6, 6.5643e4, [37.7890, 36.5990], -28940.01, -3837.2, 0, 11.0); ...
  localBuildStateSummary(true, -1.5080e6, 6.5640e4, [37.7890, 36.5990], -28940.03, -3823.3, 0, 18.0) ...
  };
caseSpec.tagList = ["periodic-wide"; "periodic-fd-anchor"; "periodic-in-tooth"];
caseSpec.selectOpt = localBuildAnchorGuardOpt(2);
caseSpec.expectedIdx = 1;
caseSpec.expectedTag = "periodic-wide";
caseSpec.reasonNotContains = "subset-anchor guard";
end


function caseSpec = localCloserToZeroToothWinsCase()
caseSpec = localEmptyCase("closer-to-zero tooth wins");
caseSpec.summaryList = { ...
  localBuildStateSummary(true, -1.5081e6, 6.5632e4, [37.7800, 36.5900], -28940.01, -3844.6, 0, 12.0); ...
  localBuildStateSummary(true, -1.5052e6, 6.5920e4, [37.7490, 36.5870], -30440.0, -3883.5, -2, 11.0); ...
  localBuildStateSummary(true, -1.5052e6, 6.5922e4, [37.7492, 36.5872], -30440.0, -3883.6, -2, 18.0) ...
  };
caseSpec.tagList = ["periodic-wide"; "periodic-fd-anchor"; "periodic-in-tooth"];
caseSpec.selectOpt = localBuildAnchorGuardOpt(2);
caseSpec.expectedIdx = 1;
caseSpec.expectedTag = "periodic-wide";
caseSpec.reasonNotContains = "subset-anchor guard";
end


function caseSpec = localEmptyCase(name)
caseSpec = struct();
caseSpec.name = name;
caseSpec.summaryList = {};
caseSpec.tagList = strings(0, 1);
caseSpec.selectOpt = [];
caseSpec.expectedIdx = NaN;
caseSpec.expectedTag = "";
caseSpec.reasonContains = "";
caseSpec.reasonNotContains = "";
end


function opt = localBuildAnchorGuardOpt(subsetAnchorIdx)
opt = struct( ...
  'enableSubsetAnchorGuard', true, ...
  'subsetAnchorIdx', subsetAnchorIdx, ...
  'subsetAnchorTag', "periodic-fd-anchor", ...
  'maxFdRefDriftHz', 100, ...
  'maxFdRateDriftHzPerSec', 100, ...
  'maxDoaDriftDeg', 0.01, ...
  'minObjGainToLeaveAnchor', 1e4, ...
  'minResidualGainToLeaveAnchor', 1e3);
end


function summary = localBuildTruthErrorSummary(isResolved, finalObj, finalResidualNorm, angleErrDeg, fdRefErrHz, fdRateErrHzPerSec, toothIdx, runTimeMs)
summary = struct();
summary.isResolved = isResolved;
summary.finalObj = finalObj;
summary.finalResidualNorm = finalResidualNorm;
summary.angleErrDeg = angleErrDeg;
summary.fdRefErrHz = fdRefErrHz;
summary.fdRateErrHzPerSec = fdRateErrHzPerSec;
summary.toothIdx = toothIdx;
summary.toothResidualHz = fdRefErrHz;
summary.runTimeMs = runTimeMs;
end


function summary = localBuildStateSummary(isResolved, finalObj, finalResidualNorm, doaParamEst, fdRefEst, fdRateEst, toothIdx, runTimeMs)
summary = struct();
summary.isResolved = isResolved;
summary.finalObj = finalObj;
summary.finalResidualNorm = finalResidualNorm;
summary.doaParamEst = reshape(doaParamEst, 1, []);
summary.fdRefEst = fdRefEst;
summary.fdRateEst = fdRateEst;
summary.angleErrDeg = 0;
summary.fdRefErrHz = 0;
summary.fdRateErrHzPerSec = 0;
summary.toothIdx = toothIdx;
summary.toothResidualHz = 0;
summary.runTimeMs = runTimeMs;
end


function textOut = localFormatNumericRow(valueVec)
valueVec = reshape(valueVec, 1, []);
cellText = arrayfun(@(x) sprintf('%.6e', x), valueVec, 'UniformOutput', false);
textOut = ['[', strjoin(cellText, ', '), ']'];
end
