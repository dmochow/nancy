function [I,optval] = solve_lsq_l1(E,R,Imax)
n=size(R,1); % number of electrodes
cvx_begin
    variable I(n);
    minimize( (R.'*I-E).'*(R.'*I-E) );
     subject to
         sum(abs([I; -sum(I)])) <= 2*Imax;
cvx_end
optval=cvx_optval;