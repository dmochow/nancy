function [cleanData,gsvs] = gsvdDenoise(noisyData,refData,Kmin)
%denoise the time-space matrix origData by removing all sources that are at
%least 1/Kmin more powerful in origData than in refData
%  
%
if nargin<3, Kmin=0.5; end % default is to remove sources twice as strong in noisyData
if size(noisyData,1)<size(noisyData,2), 
    noisyData=noisyData.'; 
    warning('Permuting origData'); 
    permuteFlag=1;
end
if size(refData,1)<size(refData,2), refData=refData.'; warning('Permuting refData'); end

[U,V,XX,C,S] = gsvd(refData,noisyData,0);
gsvs=diag(C)./diag(S);
K=sum(gsvs<Kmin);
fprintf('%d dimensions removed \n',K);
cleanData=V*S(:,K+1:end)*XX(:,K+1:end).';   
if permuteFlag, cleanData=cleanData.'; end
%figure; plot(gsvs,'*k'); % use this to determine K
    
   
end

