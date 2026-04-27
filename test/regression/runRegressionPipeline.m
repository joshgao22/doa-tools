function result = runRegressionPipeline(varargin)
%RUNREGRESSIONPIPELINE Run the slower pipeline-level regression suite.
%
% Syntax:
%   result = runRegressionPipeline()
%   result = runRegressionPipeline(runOpt)
%   result = runRegressionPipeline(Name=Value)
%
% Inputs:
%   runOpt             - optional options structure, or Name=Value pairs
%     .projectRoot     - repository root. Default: auto-resolve
%     .stopOnFailure   - true to stop on first failure. Default: true
%     .verbose         - true to print progress. Default: true
%
% Output:
%   result             - suite result structure returned by
%                        runRegressionSuite
%
% Notes:
%   - This suite is reserved for active end-to-end behavior contracts that
%     should remain automatic but are too expensive for quick/dynamic_flow.
%   - Strategy reports, plots, and checkpoint/resume smoke tests are kept
%     outside this runner: use test/dev/strategy, test/dev/probe, and
%     test/regression/perf/runRegressionPerfSmoke for those paths.

caseRelPathList = [ ...
  "test/regression/pipeline/regressionMfDoaProfileWithHealthyFd.m"; ...
  "test/regression/pipeline/regressionMfCpVsIpRepresentative.m" ...
  ];

result = runRegressionSuite("pipeline", caseRelPathList, varargin{:});
end
