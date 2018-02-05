function out=nle(dataIn)
% assumes time in dimension 2
out=dataIn.^2 - circshift(dataIn,1,2).*circshift(dataIn,-1,2);