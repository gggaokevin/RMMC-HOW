function main_worker(w)
% MAIN_WORKER  run every cell assigned to worker w (see schedule.csv), sequentially.
% Called as  matlab -batch "cd('.../code/experiments'); main_worker(w)"  with the process
% pinned to worker w's logical CPU (masks in schedule.csv / main_plan.m).

% ---- package layout ---------------------------------------------------------
%   <pkg>/experiments/main_*.m   the S1-S4 experiment programs (this file)
%   <pkg>/                       RMMC.m, MMC*.m, utils/, baselines/, third_party/
%   <pkg>/data/orl/orl_face/     ORL face database (S4 input)
%   <pkg>/results/               every result these programs write
here     = fileparts(mfilename('fullpath'));          % <pkg>/experiments
codeRoot = fileparts(here);                           % <pkg>
% ----------------------------------------------------------------------------
outdir = fullfile(codeRoot,'results');
sch = readtable(fullfile(outdir,'schedule.csv'));
sel = find(sch.worker == w);
if isempty(sel)
    fprintf('[w%d] no cells assigned\n', w); return;
end
fprintf('[w%d] %d cells, estimated %.1f min\n', w, numel(sel), sum(sch.est_sec(sel))/60);
t0 = tic;
for i = sel(:)'
    ms = strsplit(sch.methods{i}, ',');
    fprintf('[w%d] %d/%d  %s setting=%g seed=%d  {%s}\n', w, find(sel==i), numel(sel), ...
            sch.scenario{i}, sch.setting(i), sch.seed(i), sch.methods{i});
    try
        main_run(sch.scenario{i}, sch.setting(i), sch.seed(i), ms);
    catch ME
        fprintf('[w%d] CELL FAILED %s %g %d : %s\n', w, sch.scenario{i}, sch.setting(i), sch.seed(i), ME.message);
    end
    fprintf('[w%d] elapsed %.1f min\n', w, toc(t0)/60);
end
fprintf('[w%d] DONE in %.1f min\n', w, toc(t0)/60);
end
