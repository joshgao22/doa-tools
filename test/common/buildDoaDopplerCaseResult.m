function caseInfo = buildDoaDopplerCaseResult(displayName, satMode, frameMode, ...
  paramMode, dynamicMode, estResult)
%BUILDDOADOPPLERCASERESULT Pack one estimator result with display metadata.

caseInfo = struct();
caseInfo.displayName = string(displayName);
caseInfo.satMode = string(satMode);
caseInfo.frameMode = string(frameMode);
caseInfo.paramMode = string(paramMode);
caseInfo.dynamicMode = string(dynamicMode);
caseInfo.estResult = estResult;
end
