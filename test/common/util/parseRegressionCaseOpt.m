function opt = parseRegressionCaseOpt(varargin)
%PARSEREGRESSIONCASEOPT Parse common regression case options.
%
% Syntax:
%   opt = parseRegressionCaseOpt()
%   opt = parseRegressionCaseOpt(runOpt)
%   opt = parseRegressionCaseOpt(Name, Value)
%
% Supported options:
%   verbose - enable detailed case diagnostics. Default: inherited runner
%             verbose, or false when the case is executed directly.
%

opt = struct();
opt.verbose = localGetInheritedVerbose();

if nargin == 0
  return;
end

if nargin == 1 && isstruct(varargin{1})
  optIn = varargin{1};
  if isfield(optIn, 'verbose') && ~isempty(optIn.verbose)
    opt.verbose = logical(optIn.verbose);
  end
  return;
end

if mod(nargin, 2) ~= 0
  error('parseRegressionCaseOpt:InvalidInput', ...
    'Optional inputs must be one struct or Name/Value pairs.');
end

for iArg = 1:2:nargin
  name = lower(char(string(varargin{iArg})));
  value = varargin{iArg + 1};
  switch name
    case 'verbose'
      opt.verbose = logical(value);
    otherwise
      error('parseRegressionCaseOpt:UnknownOption', ...
        'Unknown option name: %s', char(string(varargin{iArg})));
  end
end
end


function verbose = localGetInheritedVerbose()
%LOCALGETINHERITEDVERBOSE Resolve runner-published verbose flag.

verbose = false;

try
  if evalin('base', "exist(''doaToolsVerbose'', ''var'')") == 1
    verbose = logical(evalin('base', 'doaToolsVerbose'));
    return;
  end
catch
end

try
  if isappdata(0, 'doaToolsVerbose')
    verbose = logical(getappdata(0, 'doaToolsVerbose'));
  end
catch
end
end
