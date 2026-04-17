function result = runRegressionPipeline(varargin)
%RUNREGRESSIONPIPELINE Run the slower pipeline-level regression suite.
%
%Syntax:
%  result = runRegressionPipeline()
%  result = runRegressionPipeline(runOpt)
%  result = runRegressionPipeline(Name=Value)
%
%Inputs:
%  runOpt             - optional options structure, or Name=Value pairs
%    .projectRoot     - repository root. Default: auto-resolve
%    .stopOnFailure   - true to stop on first failure. Default: true
%    .verbose         - true to print progress. Default: true
%
%Output:
%  result             - suite result structure returned by
%                       runRegressionSuite
%
%Notes:
%  - This suite is reserved for the higher-cost pipeline regression that
%    exercises the current subset-select -> periodic in-tooth refine flow.
%  - Keep this list short and representative. It should guard the active
%    mainline, not every exploratory scan script.

scriptRelPathList = [ ...
  "test/regression/pipeline/regressionMfSubsetToothSelectThenPeriodicRefine.m" ...
  ];

result = runRegressionSuite("pipeline", scriptRelPathList, varargin{:});
end
