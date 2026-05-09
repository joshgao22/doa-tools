function optVar = makeMfProbeOptVar(method, doaParam, fdRefHz, fdRateHzPerSec)
%MAKEMFPROBEOPTVAR Build one MF probe optimization vector.

doaParam = takeFirstTwoValues(doaParam);
if method.isKnownRate
  optVar = [doaParam(:); fdRefHz];
else
  optVar = [doaParam(:); fdRefHz; fdRateHzPerSec];
end
end
