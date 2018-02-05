function Idown = reduceMontage(E,I,R,K,Imax)

A=R.'; s=I; 
sN = [s; -sum(s)];  % minus sign on the reference electrode
A = [A, zeros(size(A(:,1)))];
N = size(A,2);
[tmp,indx] = sort(abs(sN));

T = eye(N);   % this will be the tranform (including the old ref!)
T(indx(end-K+1),:) = -1; % the Kth smallest one becomes the new reference
invT = inv(T);

indx = indx(end:-1:end-K+1); % keep the top K
sreref = invT(indx,:)*sN;     % compute the rereferenced currents for the top K
sreref = sreref-mean(sreref); % make sure they sum to zero
sK=zeros(N,1); sK(indx) = sreref;

Areref = A*T(:,indx(1:K-1)); % we only need the top K-1,

[Itmp,optval] = solve_lsq_l1(E,Areref.',Imax)
sreref(1:K-1)=Itmp; sreref(K)=-sum(sreref(1:K-1));
Idown = T(:,indx(1:K-1))*sreref(1:K-1); % now figure out what that is in the original space
