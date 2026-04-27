function result = runRegressionPerfSmoke(varargin)
%RUNREGRESSIONPERFSMOKE Run the perf-oriented regression smoke suite.
%
% Syntax:
%   result = runRegressionPerfSmoke()
%   result = runRegressionPerfSmoke(runOpt)
%   result = runRegressionPerfSmoke(Name=Value)
%
% Notes:
%   - This runner is intentionally separate from quick/dynamic_flow/pipeline.
%   - Use it when touching checkpoint/resume, summary, snapshot, or cleanup
%     orchestration around long-task perf flows.

caseRelPathList = [ ...
  "test/regression/perf/regressionDynamicPerfSmoke.m" ...
  ];

result = runRegressionSuite("perf_smoke", caseRelPathList, varargin{:});
end
