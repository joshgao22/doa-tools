function result = runRegressionDynamicFlow(varargin)
%RUNREGRESSIONDYNAMICFLOW Run the fast dynamic-flow regression guardrail suite.
%
% This suite sits between quick invariants and heavier pipeline/dev checks.
% It is intentionally kept fast, and only covers the current simplified
% dynamic-flow contracts that should be rerun after routine flow edits.
%
% Heavy representative/pipeline diagnostics such as
% probeMfSubsetToothSelectThenPeriodicRefine and
% regressionMfDoaProfileWithHealthyFd are intentionally excluded here.
% Those belong in the slower pipeline/dev layer.

caseRelPathList = [ ...
  "test/regression/branch/regressionMfSubsetSelectNoTruthLeak.m"; ...
  "test/regression/branch/regressionMfUnknownFinalSelectionRules.m"; ...
  "test/regression/branch/regressionSimpleFlowPolishTrigger.m"; ...
  "test/regression/branch/regressionSimplePeriodicReplaySeedSelection.m" ...
  ];

result = runRegressionSuite("dynamic_flow", caseRelPathList, varargin{:});
end
