function [I,optval] = solve_l1(V,R,Imax)
%I=(k*(R*R'))\V;
n=size(V,1); % number of electrodes
cvx_begin
    variable I(n);
    minimize( I'*(R*R')*I - 2*V'*I );
     subject to
         sum(abs([I; -sum(I)])) <= 2*Imax;
cvx_end
optval=cvx_optval;