function main_run(scenario, setting, seed, methods, runmode, tagroot)
% MAIN_RUN  one cell of the main comparison experiment (every method except RMMC).
%
%   main_run('s1', 9, 5)                    % all applicable methods
%   main_run('s1', 9, 5, {'MMC','HOW'})     % a subset (import gate uses this)
%   main_run('s4sp',10, 5, {}, 'uniform6','timing')   % dedicated timing pass
%
% runmode  : label written into the row ('batch_max' by default, 'uniform6' for the
%            dedicated timing pass).
% tagroot  : sub-directory of 03_results/main holding the output; 'timing' writes
%            rows only (no .mat), so the main reconstruction store is never touched.
%
% Scenarios
%   s1  : synthetic 100x200x20, 20% random missing, SNR = setting dB
%   s2  : synthetic, missing ratio = setting (0..0.9), SNR = 10 dB
%   s3  : s1 + 10 missing rows and 10 missing columns per frame
%   s4sp: ORL 112x92x400 (read from <pkg>/data/orl), 20% random + six bands
%         per frame, salt-and-pepper q = setting percent on the observed entries
%
% Caps: s1-s3 = 200 (HOW: 100 warm-up + 100 robust = 200 total; l0: 100+100),
%       s4sp  = 400 (HOW: 100 + 300; l0: 200+200).  RGNMR runs its released
%       defaults on s1-s3 only.  Every method keeps its own ground-truth-free
%       stopping rule; no ground truth is ever used for a stopping decision.
%
% The data construction of s1-s3 is copied VERBATIM from the frozen driver
% RMMC_V5/code/exp_v4/run_v5_main.m (rng(0) -> eigs(gaussgram) -> normalisation ->
% mk_mask(seed) -> gm_noise(seed+1)); s4sp follows 03_code/s4_ab_cell.m.
%
% Output: everything under <pkg>/results
%   mats/<scenario>/<setting>_seed<k>.mat        X0, In, Om, Xh_<METHOD> (appended)
%   cells/<scenario>/<setting>_seed<k>.csv       one row per method
%   metrics/...  traces/...  logs/...

if nargin < 4 || isempty(methods), methods = {}; end
if nargin < 5 || isempty(runmode), runmode = 'batch_max'; end
if nargin < 6 || isempty(tagroot), tagroot = ''; end

% ---- package layout ---------------------------------------------------------
%   <pkg>/experiments/main_*.m   the S1-S4 experiment programs (this file)
%   <pkg>/                       RMMC.m, MMC*.m, utils/, baselines/, third_party/
%   <pkg>/data/orl/orl_face/     ORL face database (S4 input)
%   <pkg>/results/               every result these programs write
here     = fileparts(mfilename('fullpath'));          % <pkg>/experiments
codeRoot = fileparts(here);                           % <pkg>
% ----------------------------------------------------------------------------
if isempty(tagroot)
    outdir = fullfile(codeRoot,'results');
else
    outdir = fullfile(codeRoot,'results',tagroot);
end
addpath(here); addpath(codeRoot);
addpath(fullfile(codeRoot,'utils'));
addpath(fullfile(codeRoot,'baselines','RMC_HOW'));
addpath(fullfile(codeRoot,'baselines','RTCUR'));
addpath(fullfile(codeRoot,'baselines','RTCUR','tools'));
addpath(fullfile(codeRoot,'baselines','TriTD'));
addpath(genpath(fullfile(codeRoot,'baselines','RGNMR')));
addpath(genpath(fullfile(codeRoot,'third_party','tensor_toolbox')));
for d = {'mats','cells','metrics','traces','logs'}
    p = fullfile(outdir,d{1},scenario); if ~exist(p,'dir'), mkdir(p); end
end

ORDER = {'MMC','RMMC','HOW','L0IQR','TriTD','RTCUR','RGNMR'};
if isempty(methods), methods = ORDER; end
methods = ORDER(ismember(ORDER, methods));           % fixed order, RTCUR then RGNMR

stem = sprintf('%g_seed%d', setting, seed);
logf = fullfile(outdir,'logs',scenario,[stem '.txt']);
lfid = fopen(logf,'a');
say  = @(varargin) (fprintf(lfid,varargin{:}) && fprintf(varargin{:}));

say('\n==== main_run %s setting=%g seed=%d methods={%s} %s ====\n', ...
    scenario, setting, seed, strjoin(methods,','), datestr(now,'yyyy-mm-dd HH:MM:SS'));

% ---------------- settings --------------------------------------------------
tau = 0.10; h = 10; OpRand = 0.80; STOP_TOL = 1e-4;
switch scenario
    case {'s1','s2','s3'}
        m=100; n=200; K=20; rk=20; d=20; lam=5e-4; It=200; region = false;
        rk_tc=[20 20 20]; rk_tri=20; howRob=100; l0a=100; l0b=100;
        if strcmp(scenario,'s2'), SNR = 10; Op = 1-setting; else, SNR = setting; Op = OpRand; end
        pat = 1; if strcmp(scenario,'s3'), pat = 3; end
        rng(0);
        [U0,~]=eigs(gaussgram(m,8),rk,'largestabs');
        [V0,~]=eigs(gaussgram(n,10),rk,'largestabs');
        X0=zeros(m,n,K); for k=1:K, X0(:,:,k)=U0*randn(rk,rk)*V0'; end
        X0=X0/mean(abs(X0(:)));
        Om = mk_mask(m,n,K,seed,Op,pat);
        In = (X0 + gm_noise(X0,Om,SNR,tau,h,seed+1)) .* Om;
        Oc = ~Om;
    case 's4sp'
        m=112; n=92; K=400; rk=40; d=90; lam=0.1; It=400; region = true;
        rk_tc=[90 90 40]; rk_tri=20; howRob=300; l0a=200; l0b=200;
        X0 = load_orl(fullfile(codeRoot,'data'));
        Om = mk_mask4(m,n,K,seed,OpRand);
        In = X0 .* Om;
        rng(seed+1);
        SPm = (rand(m,n,K) < setting/100) & Om;
        In(SPm) = round(rand(nnz(SPm),1));
        Oc = ~Om; pat = 2; SNR = NaN; Op = OpRand;
    otherwise
        error('unknown scenario %s', scenario);
end
say('data %dx%dx%d rk=%d d=%d lam=%g It=%d |Om|=%.4f\n', m,n,K,rk,d,lam,It,nnz(Om)/numel(Om));

matfile = fullfile(outdir,'mats',scenario,[stem '.mat']);
csvfile = fullfile(outdir,'cells',scenario,[stem '.csv']);
if exist(matfile,'file')
    S = load(matfile,'gt_meta'); gt_meta = S.gt_meta;
else
    gt_meta = struct('scenario',scenario,'setting',setting,'seed',seed,'m',m,'n',n,'K',K, ...
                     'rk',rk,'d',d,'lam',lam,'It',It,'SNR',SNR,'Op',Op,'pattern',pat, ...
                     'rk_tc',rk_tc,'rk_tri',rk_tri,'tau',tau,'h',h,'mask','v5_main/frozen');
    save(matfile,'X0','In','Om','gt_meta','-v7');
end
if ~exist(csvfile,'file')
    fid = fopen(csvfile,'w');
    fprintf(fid,'cell_id,scenario,setting,seed,method,codehash,cap,iters,iters_stage2,stop_reason,run_mode,source,seconds,apsnr_unit,apsnr_matlab,assim,relerr,rmse,relerr_on,psnr_min_frame,band_psnr,band_bias,rmiss_psnr,rmiss_bias,note\n');
    fclose(fid);
end
done = read_done(csvfile);

% ---------------- run -------------------------------------------------------
for mi = 1:numel(methods)
    M = methods{mi};
    if done.isKey(M) && ~strcmp(done(M),'error')
        say('  %-6s already %s -> skip\n', M, done(M)); continue;
    end
    hash = file_hash(which_method(M, codeRoot));
    t = tic; Xh = []; trace = []; trname = ''; stop_reason = 'cap'; it = NaN; it2 = NaN; itrace = {}; snap_ap = [];
    try
        switch M
            case 'MMC'
                [Xh, Err, it] = MMC_half(In, Oc, X0, d, d, It, lam, STOP_TOL);
                if it < It, stop_reason='tol'; end
                trace = Err(1:it); trname = 'relerr';
            case 'RMMC'
                % V6 (observed-support IQR): lambda switched by the data family
                % (synthetic 100x200: 5e-3; ORL 112x92: 0.1); ONE penalty pair
                % (rho1,rho2)=(9.05,1.45) and eta=(2.5,3.5) for both families.
                [lam_r, rho1_r, rho2_r] = rmmc_defaults(m,n);
                snap = [50 100]; if It >= 400, snap = [50 100 200]; end
                [Xh, O_r, Err, it, Snap] = RMMC(In, Oc, X0, d, d, It, 2.5, 3.5, lam_r, ...
                                                 0, STOP_TOL, snap, rho1_r, rho2_r, false);
                trace = Err(1:max(abs(it),1)); trname = 'relerr';
                if it < 0
                    % degenerate scale over the observed support (d = 0):
                    % the cell is NOT executable, no reconstruction is stored
                    stop_reason = 'empty_iqr';
                    fail_note = sprintf('d collapsed at round %d; not executable', abs(it));
                    say('  RMMC   NOT EXECUTABLE (%s)\n', fail_note);
                    row(fullfile(outdir,'cells',scenario,[stem '.csv']), scenario, setting, seed, ...
                        M, hash, It, abs(it), NaN, stop_reason, toc(t), nan(1,11), 'new', runmode);
                    continue;
                end
                if it < It, stop_reason='tol'; end
                snap_ap = zeros(1,numel(Snap));
                for q = 1:numel(Snap), snap_ap(q) = mean(ap_unit(X0, Snap(q).USV)); end
                say('  RMMC  snapshots: %s -> APSNR %s\n', mat2str([Snap.iter]), mat2str(round(snap_ap,2)));
            case 'HOW'
                Xh = X0; it = 0; ncap = 0;
                for k=1:K
                    rng(7000+k);
                    [Xh(:,:,k),~,~,~,info]=HOW_RMC(In(:,:,k),double(Om(:,:,k)),rk,howRob,2,2.3);
                    it = it + info.iter; if isfield(info,'warmup'), it = it + info.warmup; end
                    if info.iter < howRob, ncap = ncap+1; end
                    itrace{end+1} = info.loc; %#ok<AGROW>
                end
                if ncap > 0, stop_reason='tol'; end
            case 'L0IQR'
                Xh = X0; it = 0;
                l0opts = struct('xi1',2,'c0',2,'tol',1e-4,'stopMode','relchg','stop1',true);
                for k=1:K
                    rng(7000+k);
                    [Xh(:,:,k),~,~,~,~,itk]=L0_IQR(X0(:,:,k),In(:,:,k),double(Om(:,:,k)),rk,l0a,l0b,1e-5,l0opts);
                    it = it + itk;
                end
                it2 = l0b;
                if it < K*(l0a+l0b), stop_reason='tol'; end
            case 'TriTD'
                o = struct('mu',1e-3,'rho',1.25,'lambda',1.8,'lambda2',1e-3,'maxIter',It,'tol',1e-5,'disp',0);
                rng(0);
                [A,B,C,~,errHist]=triple_decomp_ADMM(In,rk_tri,o); Xh=triple_product(A,B,C);
                it = numel(errHist); trace = errHist(:); trname = 'residual';
                if it < It, stop_reason='tol'; end
            case 'RTCUR'
                if pat==1 && Op>=0.8
                    F = 3*max(abs(In(:))); if F==0, F=3; end
                    para = []; para.max_iter = It; para.epsilon = 1e-5;
                else
                    F = 2*max(abs(In(:))); if F==0, F=2; end
                    para = []; para.zeta=1; para.gamma=1; para.max_iter=It; para.epsilon=1e-5;
                end
                rng(11);
                InF = In; InF(Oc) = F;
                [Lc,Xs,~,err]=RTCUR_fc(tensor(InF),rk_tc,para); Xh=double(ttm(Lc,Xs));
                it = nnz(err>-1); trace = err(1:it); trname = 'err';
                if it < It, stop_reason='tol'; end
            case 'RGNMR'
                Xh = X0; it = 0; nf = 0;
                for k=1:K
                    Omf = double(Om(:,:,k)); [ii,jj]=find(Omf); omega=[ii jj];
                    onum = max(round(tau*nnz(Omf)),1);
                    rng(7000+k);
                    [Xh(:,:,k),flag,relRes]=RGNMR(In(:,:,k),omega,rk,onum,rgnmr_opts());
                    it = it + numel(relRes); nf = nf + double(flag~=0);
                    itrace{end+1} = relRes; %#ok<AGROW>
                end
                if nf > 0, stop_reason='flag_nonzero'; end
            otherwise
                error('unknown method');
        end
        el = toc(t);
    catch ME
        el = toc(t);
        say('  %-6s FAIL: %s\n', M, ME.message);
        row(fullfile(outdir,'cells',scenario,[stem '.csv']), scenario, setting, seed, M, hash, It, NaN, NaN, 'error', el, nan(1,11), 'new', runmode);
        continue;
    end

    Mt = per_frame(X0, Xh);
    if region
        [rb, rn] = region_metrics(X0, Xh, Om);
    else
        rb = [NaN NaN NaN NaN]; rn = rb;
    end
    say('  %-6s iters=%-5d it2=%-4g (%.1f s) APSNR=%.2f ASSIM=%.4f relerr=%.4f [%s]\n', ...
        M, it, it2, el, Mt.apsnr, Mt.assim, Mt.relerr, stop_reason);

    mrow = struct('hash',hash,'cap',It,'iters',it,'iters2',it2,'stop',stop_reason, ...
                  'sec',el,'apsnr',Mt.apsnr,'apsnr_m',Mt.apsnr_m,'assim',Mt.assim, ...
                  'relerr',Mt.relerr,'rmse',Mt.rmse,'relerr_on',Mt.relerr_on, ...
                  'psnr_min',Mt.psnr_min,'band_psnr',rb(1),'band_bias',rn(1), ...
                  'snap_ap',snap_ap,'params',gt_meta);
    res = struct('Xh',Xh,'meta',mrow,'trace',trace,'itrace',{itrace});
    W = struct(); W.(['res_' M]) = res; %#ok<STRNU>
    if isempty(tagroot)                      % the timing pass stores rows only
        save(matfile, '-struct','W', '-v7','-append');
    end

    fid = fopen(fullfile(outdir,'metrics',scenario,sprintf('%s_%s.csv',stem,M)),'w');
    fprintf(fid,'frame,rmse,relerr,mse,ssim,psnr_matlab\n');
    for k=1:K
        fprintf(fid,'%d,%.10g,%.10g,%.10g,%.10g,%.10g\n', k, Mt.rmse_f(k), Mt.relerr_f(k), Mt.mse_f(k), Mt.ssim_f(k), Mt.psnr_m_f(k));
    end
    fclose(fid);
    if ~isempty(trace)
        fid = fopen(fullfile(outdir,'traces',scenario,sprintf('%s_%s.csv',stem,M)),'w');
        fprintf(fid,'iter,%s\n', trname);
        for q=1:numel(trace), fprintf(fid,'%d,%.10g\n', q, trace(q)); end
        fclose(fid);
    elseif ~isempty(itrace)
        fid = fopen(fullfile(outdir,'traces',scenario,sprintf('%s_%s.csv',stem,M)),'w');
        for q=1:numel(itrace)
            v = itrace{q};
            fprintf(fid,'frame%d,%s\n', q, strjoin(cellstr(num2str(v(:),'%.10g')),' '));
        end
        fclose(fid);
    end
    row(fullfile(outdir,'cells',scenario,[stem '.csv']), scenario, setting, seed, M, hash, It, it, it2, stop_reason, el, ...
        [Mt.apsnr, Mt.apsnr_m, Mt.assim, Mt.relerr, Mt.rmse, Mt.relerr_on, Mt.psnr_min, rb(1), rn(1), rb(2), rn(2)], 'new', runmode);
    say('  %-6s saved (hash %s)\n', M, hash);
end
fclose(lfid);
end

% ============================ helpers =======================================
function f = which_method(M, codeRoot)
switch M
    case 'MMC',    f = fullfile(codeRoot,'MMC_half.m');
    case 'RMMC',   f = fullfile(codeRoot,'RMMC.m');
    case 'HOW',    f = fullfile(codeRoot,'baselines','RMC_HOW','HOW_RMC.m');
    case 'L0IQR',  f = fullfile(codeRoot,'baselines','RMC_HOW','L0_IQR.m');
    case 'TriTD',  f = fullfile(codeRoot,'baselines','TriTD','triple_decomp_ADMM.m');
    case 'RTCUR',  f = fullfile(codeRoot,'baselines','RTCUR','RTCUR_fc.m');
    case 'RGNMR',  f = fullfile(codeRoot,'baselines','RGNMR','RGNMR.m');
end
end

function d = read_done(csvfile)
d = containers.Map('KeyType','char','ValueType','char');
if ~exist(csvfile,'file'), return; end
T = readcell(csvfile);
for i = 2:size(T,1)
    if size(T,2) >= 10 && ischar(T{i,5}) && ischar(T{i,10}), d(T{i,5}) = T{i,10}; end
end
end

function row(csvfile, scenario, setting, seed, M, hash, cap, it, it2, stop_reason, sec, v, source, runmode)
% v = [apsnr, apsnr_m, assim, relerr, rmse, relerr_on, psnr_min, band_psnr, band_bias, rmiss_psnr, rmiss_bias]
fid = fopen(csvfile,'a');
fprintf(fid,'%s_%g_%d,%s,%g,%d,%s,%s,%d,%g,%g,%s,%s,%s,%.2f', ...
    scenario, setting, seed, scenario, setting, seed, M, hash, cap, it, it2, stop_reason, runmode, source, sec);
fprintf(fid,',%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,-\n', ...
    v(1),v(2),v(3),v(4),v(5),v(6),v(7),v(8),v(9),v(10),v(11));
fclose(fid);
end

function a = ap_unit(X0,Xh)
% per-frame unit-peak PSNR (the frozen table's convention)
K=size(X0,3); a=zeros(1,K);
for k=1:K, a(k) = 10*log10(1/mean((X0(:,:,k)-Xh(:,:,k)).^2,'all')); end
end

function M = per_frame(X0, Xh)
K = size(X0,3);
M.rmse_f=zeros(1,K); M.relerr_f=zeros(1,K); M.mse_f=zeros(1,K); M.ssim_f=zeros(1,K); M.psnr_m_f=zeros(1,K);
[m,n,~] = size(X0);
for k=1:K
    E = X0(:,:,k)-Xh(:,:,k);
    M.mse_f(k) = mean(E(:).^2);
    M.rmse_f(k) = norm(E,'fro')/sqrt(m*n);
    M.relerr_f(k) = norm(E,'fro')/norm(X0(:,:,k),'fro');
    M.ssim_f(k) = ssim(X0(:,:,k),Xh(:,:,k));
    M.psnr_m_f(k) = psnr(X0(:,:,k),Xh(:,:,k));
end
M.rmse=mean(M.rmse_f); M.relerr=mean(M.relerr_f);
M.apsnr = mean(10*log10(1./max(M.mse_f,eps)));
M.apsnr_m = mean(M.psnr_m_f); M.assim = mean(M.ssim_f); M.psnr_min = min(M.psnr_m_f);
on = zeros(1,K);
for k=1:K
    msk = X0(:,:,k) ~= 0 | Xh(:,:,k) ~= 0;   % unused fallback; relerr_on recomputed below
    on(k) = NaN; %#ok<NASGU>
end
M.relerr_on = NaN;
end

function [ps, bi] = region_metrics(X0, Xh, Om)
% band / rmiss / bnd / obs, unit-peak PSNR and signed bias (S4 only)
[m,n,K] = size(Om);
band = false(m,n,K); bnd = false(m,n,K);
for k=1:K
    Mk = Om(:,:,k);
    bc = find(all(~Mk,1)); br = find(all(~Mk,2));
    if ~isempty(bc), band(:,bc,k)=true; end
    if ~isempty(br), band(br,:,k)=true; end
    if any(any(band(:,:,k)))
        nb = conv2(double(band(:,:,k)),ones(5),'same')>0;
        bnd(:,:,k) = Mk & nb;
    end
end
masks = {band, ~Om & ~band, bnd, Om & ~bnd};
ps = zeros(1,4); bi = zeros(1,4);
for r=1:4
    p=zeros(1,K); b=zeros(1,K);
    for k=1:K
        msk = masks{r}(:,:,k);
        if ~any(msk(:)), p(k)=NaN; b(k)=NaN; continue; end
        e = Xh(:,:,k)-X0(:,:,k); ev = e(msk);
        p(k) = 10*log10(1/mean(ev.^2)); b(k) = mean(ev);
    end
    ps(r)=mean(p,'omitnan'); bi(r)=mean(b,'omitnan');
end
end

function X0 = load_orl(dataRoot)
orl = fullfile(dataRoot,'orl','orl_face');
sd = dir(orl); sd = sd([sd.isdir] & ~startsWith({sd.name},'.'));
nums = zeros(numel(sd),1);
for i=1:numel(sd), nums(i)=sscanf(sd(i).name(2:end),'%d'); end
[~,ord]=sort(nums); sd=sd(ord);
X0 = zeros(112,92,numel(sd)*10); t=0;
for s=1:numel(sd)
    fl = dir(fullfile(sd(s).folder,sd(s).name,'*.pgm')); nn=zeros(numel(fl),1);
    for i=1:numel(fl), nn(i)=sscanf(fl(i).name,'%d'); end
    [~,o2]=sort(nn); fl=fl(o2);
    for k=1:numel(fl), t=t+1; X0(:,:,t)=im2double(imread(fullfile(fl(k).folder,fl(k).name))); end
end
end

function Om = mk_mask4(m,n,K,seed,Op)
% ORL scenario: 20% random + six bands per frame (widths 2,1,1,1,1,1), per frame
rng(seed);
Om = rand(m,n,K) < Op;
for k=1:K
    wc=[2 1 1 1 1 1]; wr=wc;
    c0=randi([1 n-2],1,6); r0=randi([1 m-2],1,6);
    for j=1:6
        cc=c0(j):min(n,c0(j)+wc(j)-1); rr=r0(j):min(m,r0(j)+wr(j)-1);
        Om(:,cc,k)=false; Om(rr,:,k)=false;
    end
end
end

% ---- verbatim from the frozen driver ---------------------------------------
function Om = mk_mask(m,n,K,seed,Op,pattern)
rng(seed);
Om = rand(m,n,K) < Op;
if pattern==3
    for k=1:K
        rr=randperm(m,10); cc=randperm(n,10);
        Om(rr,:,k)=false; Om(:,cc,k)=false;
    end
end
end
function noise = gm_noise(X, Om, SNR, tau, h, seed)
rng(seed);
[m,n,K]=size(X); N=nnz(Om);
sig_pow = sum(abs(X(Om(:))).^2)/N;
s_v2    = sig_pow / 10^(SNR/10);
s1      = sqrt(s_v2 / ((1-tau) + tau*h^2));
s2      = h*s1;
flag    = rand(m,n,K) < (1-tau);
noise   = s1*randn(m,n,K).*flag + s2*randn(m,n,K).*(1-flag);
noise   = noise .* Om;
end
function G = gaussgram(dim, tau)
x=(1:dim)'; G = exp(-(x-x').^2/(2*tau^2));
end
function opts = rgnmr_opts()
opts.init_option=0; opts.max_outer_iter=100; opts.max_inner_iter=2000;
opts.LSQR_smart_tol=1; opts.stop_relRes=1e-15; opts.stop_relDiff=-1; opts.stop_D_diff=-1; opts.verbose=0;
end
function h = file_hash(f)
b = fread(fopen(f,'r'), Inf, '*uint8');
h = sprintf('%08x', mod(sum(double(b) .* (1:numel(b))'), 2^32));
end

