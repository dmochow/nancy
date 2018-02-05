function Il1 = reciprocateL1(V,head,Imax)
% wrapper to solve_l1 that tries to overcome failed cvx runs where no
% solution was found
%optval=NaN;
%maxIters=100;
%iter=1;
%while isnan(optval) && iter<maxIters

cvx_solver SeDuMi
[I,optval] = solve_l1(V,head.R,Imax);
if sum(isnan(I))
    %cvx_solver SeDuMi
    cvx_solver SDPT3 
    [I,optval] = solve_l1(V,head.R,Imax);
end
Il1=I;
cvx_solver SDPT3  % go back to default solver

%    iter=iter+1;
%end
Il1=I;