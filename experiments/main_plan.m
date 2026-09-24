function main_plan(nW)
% MAIN_PLAN(nW)  build the missing-work schedule for nW workers.
%   nW = 6  : one worker per PHYSICAL P-core (masks 0x3,0xC,0x30,0xC0,0x300,0xC00),
%             the uniform 6-way-independent standard;
%   nW = 16 (default): the 12 P-threads plus the 4 E-cores (E workers get cheap
%             cells only).
% MAIN_PLAN  build the cell list, the cost estimates and the 16-worker schedule.
%
% Writes <pkg>/results/plan_est.csv and <pkg>/results/schedule.csv, honouring an
% existing <pkg>/results/_gate/import_list.csv (rows taken from the frozen
% 3-trial runs; those (cell,method) pairs are not scheduled).
%
% 16 workers: logical CPUs 0-11 = the six P-cores (two threads each), 12-15 = the
% four E-cores.  E-core workers only receive cheap cells (est <= 150 s), because an
% E-core is measured ~2.3-2.7x slower than a P-core.  Long cells always go to
% P-thread workers.  Assignment is greedy LPT on the estimated cost.

% ---- package layout ---------------------------------------------------------
%   <pkg>/experiments/main_*.m   the S1-S4 experiment programs (this file)
%   <pkg>/                       RMMC.m, MMC*.m, utils/, baselines/, third_party/
%   <pkg>/data/orl/orl_face/     ORL face database (S4 input)
%   <pkg>/results/               every result these programs write
here     = fileparts(mfilename('fullpath'));          % <pkg>/experiments
codeRoot = fileparts(here);                           % <pkg>
% ----------------------------------------------------------------------------
outdir = fullfile(codeRoot,'results');
if ~exist(outdir,'dir'), mkdir(outdir); end
gate = fullfile(outdir,'_gate'); if ~exist(gate,'dir'), mkdir(gate); end

ALL6 = {'MMC','RMMC','HOW','L0IQR','TriTD','RTCUR','RGNMR'};
ALL5 = {'MMC','RMMC','HOW','L0IQR','TriTD','RTCUR'};
scen = { 's1',   [0 3 6 9 12 15],        ALL6
         's2',   0:0.1:0.9,              ALL6
         's3',   [0 3 6 9 12 15],        ALL6
         's4sp', [5 10 15],              ALL5 };
seeds = 5:9;

% ---- optional import list from the gate ------------------------------------
imp = containers.Map('KeyType','char','ValueType','logical');
fimp = fullfile(gate,'import_list.csv');
fro  = fullfile(fileparts(codeRoot),'data','archive','v5_main_frozen');  % author's archive only
if exist(fimp,'file') && exist(fro,'dir')
    T = readcell(fimp);
    for i = 2:size(T,1)
        imp(sprintf('%s|%g|%d|%s', T{i,1}, T{i,2}, T{i,3}, T{i,4})) = true;
    end
    fprintf('import list: %d rows\n', imp.Count);
    materialise_imports(fro, outdir, T);
elseif exist(fimp,'file')
    fprintf(['import list present but the frozen rows are missing (%s):\n' ...
             '  ignored - every cell is scheduled and run from scratch\n'], fro);
end

% ---- what is still missing (idempotent: re-run after any partial pass) -----
have = containers.Map('KeyType','char','ValueType','logical');
for f = dir(fullfile(outdir,'cells','*','*.csv'))'
    C = readcell(fullfile(f.folder,f.name));
    if size(C,1) < 2, continue; end
    hdr = C(1,:);
    im = find(strcmp(hdr,'method'),1); is = find(strcmp(hdr,'stop_reason'),1);
    ic = find(strcmp(hdr,'scenario'),1);
    for j = 2:size(C,1)
        if ~ischar(C{j,im}), continue; end
        sr = C{j,is}; if ~ischar(sr), sr = ''; end
        if strcmp(sr,'error'), continue; end                 % failed -> retry
        key = sprintf('%s|%g|%d|%s', ch(C{j,ic}), C{j,3}, C{j,4}, C{j,im});   % scenario|setting|seed|method
        have(key) = true;
    end
end
for k = imp.keys, have(k{1}) = true; end

% ---- plan ------------------------------------------------------------------
rows = {};
for si = 1:size(scen,1)
    sc = scen{si,1}; axis = scen{si,2}; meth = scen{si,3};
    for a = axis
        for sd = seeds
            keep = {};
            for mi = 1:numel(meth)
                if isKey(have, sprintf('%s|%g|%d|%s', sc, a, sd, meth{mi})), continue; end
                keep{end+1} = meth{mi}; %#ok<AGROW>
            end
            if isempty(keep), continue; end
            est = 0;
            for mi = 1:numel(keep), est = est + est_cost(sc, a, keep{mi}); end
            rows(end+1,:) = {sc, a, sd, strjoin(keep,','), est, numel(keep)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'scenario','setting','seed','methods','est_sec','nmethod'});
writetable(T, fullfile(outdir,'plan_est.csv'));
fprintf('cells to run: %d ; estimated total %.0f s (%.2f h single-thread P-core)\n', ...
        height(T), sum(T.est_sec), sum(T.est_sec)/3600);

% ---- assignment: greedy LPT, E-core workers cheap only ---------------------
if nargin < 1 || isempty(nW), nW = 16; end
if nW == 6
    masks = [3 12 48 192 768 3072];            % 0x3 .. 0xC00: one per physical P-core
    isE   = false(1,6);
else
    nW = 16;
    masks  = [1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768];  % logical 0..15
    isE    = false(1,nW); isE(13:16) = true;
end
EMAX   = 150;                                   % E-core workers: cells up to this est
load   = zeros(1,nW);
assign = zeros(height(T),1);
[~, ord] = sort(T.est_sec, 'descend');
for ii = 1:numel(ord)
    i = ord(ii);
    cand = 1:nW;
    if T.est_sec(i) > EMAX, cand = cand(~isE(cand)); end
    [~, j] = min(load(cand)); w = cand(j);
    assign(i) = w; load(w) = load(w) + T.est_sec(i);
end
T.worker = assign; T.mask = masks(assign)';
writetable(T, fullfile(outdir,'schedule.csv'));
fprintf('per-worker estimated load (s):\n');
for w = 1:nW
    tag = 'P'; if isE(w), tag = 'E'; end
    fprintf('  w%-2d %s mask=0x%X  est %6.0f s (%5.2f h)  cells=%d\n', ...
            w, tag, masks(w), load(w), load(w)/3600, sum(assign==w));
end
fprintf('LPT makespan estimate: %.2f h\n', max(load)/3600);
end

% ---------------------------------------------------------------------------
function materialise_imports(fro, outdir, T)
% Write the imported (cell,method) rows into the same cells/*.csv files, reading
% the metrics from the frozen 3-trial CSVs.  Already-present methods are skipped,
% so re-running the plan never duplicates a row.
F = {};
for f = dir(fullfile(fro,'s*.csv'))'
    F{end+1} = readcell(fullfile(fro,f.name)); %#ok<AGROW>
end
n = 0;
for i = 2:size(T,1)
    sc = T{i,1}; av = T{i,2}; sd = T{i,3}; m = T{i,4};
    d = fullfile(outdir,'cells',sc); if ~exist(d,'dir'), mkdir(d); end
    cf = fullfile(d, sprintf('%g_seed%d.csv', av, sd));
    if exist(cf,'file')
        C = readcell(cf);
        have = false;
        for j = 2:size(C,1)
            if ischar(C{j,5}) && strcmp(C{j,5}, m), have = true; end
        end
        if have, continue; end
    else
        fid = fopen(cf,'w');
        fprintf(fid,'cell_id,scenario,setting,seed,method,codehash,cap,iters,iters_stage2,stop_reason,run_mode,source,seconds,apsnr_unit,apsnr_matlab,assim,relerr,rmse,relerr_on,psnr_min_frame,band_psnr,band_bias,rmiss_psnr,rmiss_bias,note\n');
        fclose(fid);
    end
    r = [];
    for q = 1:numel(F)
        TT = F{q};
        for j = 2:size(TT,1)
            if strcmp(TT{j,1},sc) && abs(TT{j,2}-av) < 1e-9 && TT{j,3}==sd && strcmp(TT{j,4},m)
                r = TT(j,:); break;
            end
        end
        if ~isempty(r), break; end
    end
    if isempty(r), fprintf('  !! frozen row not found for %s %g %d %s\n', sc, av, sd, m); continue; end
    rel = r{6}; ps = r{8}; aps = r{9}; assim = r{10}; rmse = r{11}; it = r{12}; cap = r{5};
    fid = fopen(cf,'a');
    fprintf(fid,'%s_%g_%d,%s,%g,%d,%s,%s,%d,%d,NA,%s,v5exp_batched,v5exp_import,NA,%.6f,NA,%.6f,%.6f,%.6f,NA,NA,NA,NA,NA,imported_from_frozen\n', ...
        sc, av, sd, sc, av, sd, m, 'v5exp', cap, it, 'import', aps, assim, rel, rmse);
    fclose(fid);
    n = n + 1;
end
fprintf('materialised %d imported rows into cells/*.csv\n', n);
end

% ---------------------------------------------------------------------------
% ---------------------------------------------------------------------------
function s = ch(x)
if ischar(x), s = x; elseif isstring(x), s = char(x); else, s = ''; end
end
function s = est_cost(scenario, setting, method)
% seconds on one P-core, from the measured isolated/batched times
switch scenario
    case 's1'
        base = struct('MMC',20,'RMMC',30,'HOW',4,'L0IQR',5,'TriTD',66,'RTCUR',6,'RGNMR',125);
    case 's3'
        base = struct('MMC',20,'RMMC',30,'HOW',5,'L0IQR',5,'TriTD',66,'RTCUR',32,'RGNMR',110);
    case 's2'
        base = struct('MMC',20,'RMMC',30,'HOW',5,'L0IQR',5,'TriTD',65,'RTCUR',40,'RGNMR',300);
        if setting >= 0.65 && setting <= 0.75, base.RGNMR = 2600; end     % the 70% blow-up
        if setting >= 0.80, base.RGNMR = 25; end                          % early flag exit
        if setting >= 0.80, base.RMMC  = 2;  end                          % d collapses: abort
        if setting <= 0.2,  base.RGNMR = 300; end
    case 's4sp'
        base = struct('MMC',170,'RMMC',330,'HOW',170,'L0IQR',210,'TriTD',426,'RTCUR',1800,'RGNMR',0);
    otherwise
        error('scenario');
end
s = base.(method);
end
