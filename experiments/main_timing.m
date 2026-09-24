function main_timing(c)
% MAIN_TIMING  dedicated timing pass (D5): re-measure one representative cell per
% scenario under a UNIFORM condition -- one process per physical P-core, nothing
% else running -- and record the runtimes with run_mode = 'uniform6'.
%
%   main_timing(1)   s1    @ 9 dB,  seed 5
%   main_timing(2)   s2    @ 50%,   seed 5
%   main_timing(3)   s3    @ 9 dB,  seed 5
%   main_timing(4)   s4sp  @ q=10,  seed 5
%
% The four calls are launched concurrently, each pinned to
% a different physical P-core, so every method of every representative cell sees
% the same uncontended single-core condition.  Rows land in
% <pkg>/results/timing/cells/<scenario>/<setting>_seed5.csv (run_mode=uniform6)
% and are the numbers to quote in the paper.

if nargin < 1 || isempty(c), c = 1; end
switch c
    case 1, sc = 's1';   st = 9;    sd = 5;
    case 2, sc = 's2';   st = 0.5;  sd = 5;
    case 3, sc = 's3';   st = 9;    sd = 5;
    case 4, sc = 's4sp'; st = 10;   sd = 5;
    otherwise, error('cell index 1..4');
end
fprintf('[timing] cell %d : %s setting=%g seed=%d\n', c, sc, st, sd);
main_run(sc, st, sd, {}, 'uniform6', 'timing');
end
