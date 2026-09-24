function main_report()
% MAIN_REPORT  aggregate <pkg>/results/cells/**/*.csv into tidy tables.
%
% Tolerant by design: columns are located by NAME, missing values (NaN, 'NA',
% imported rows without a runtime) are skipped, and RMMC is simply absent until
% the algorithm is frozen -- the same script then picks it up unchanged.
%
% Output: <pkg>/results/tables/all_rows.csv
%         <pkg>/results/tables/by_method_<metric>_<scenario>.csv
%         <pkg>/results/tables/summary.txt

% ---- package layout ---------------------------------------------------------
%   <pkg>/experiments/main_*.m   the S1-S4 experiment programs (this file)
%   <pkg>/                       RMMC.m, MMC*.m, utils/, baselines/, third_party/
%   <pkg>/data/orl/orl_face/     ORL face database (S4 input)
%   <pkg>/results/               every result these programs write
here     = fileparts(mfilename('fullpath'));          % <pkg>/experiments
codeRoot = fileparts(here);                           % <pkg>
% ----------------------------------------------------------------------------
outdir = fullfile(codeRoot,'results');
tdir = fullfile(outdir,'tables'); if ~exist(tdir,'dir'), mkdir(tdir); end

files = dir(fullfile(outdir,'cells','*','*.csv'));
sc = {}; st = []; sd = []; me = {}; src = {}; cap = []; it = []; stop = {}; sec = []; aps = []; asm = []; rel = []; rms = []; hsh = {};
for i = 1:numel(files)
    T = readcell(fullfile(files(i).folder, files(i).name));
    if size(T,1) < 2, continue; end
    hdr = T(1,:);
    ix = @(nm) find(strcmp(hdr, nm), 1);
    for j = 2:size(T,1)
        r = @(nm) getv(T{j, ix(nm)});
        sc{end+1,1} = ch(T{j,ix('scenario')}); %#ok<AGROW>
        st(end+1,1) = r('setting');            %#ok<AGROW>
        sd(end+1,1) = r('seed');               %#ok<AGROW>
        me{end+1,1} = ch(T{j,ix('method')});   %#ok<AGROW>
        src{end+1,1} = ch(T{j,ix('source')});  %#ok<AGROW>
        cap(end+1,1) = r('cap');               %#ok<AGROW>
        it(end+1,1) = r('iters');              %#ok<AGROW>
        stop{end+1,1} = ch(T{j,ix('stop_reason')}); %#ok<AGROW>
        sec(end+1,1) = r('seconds');           %#ok<AGROW>
        aps(end+1,1) = r('apsnr_unit');        %#ok<AGROW>
        asm(end+1,1) = r('assim');             %#ok<AGROW>
        rel(end+1,1) = r('relerr');            %#ok<AGROW>
        rms(end+1,1) = r('rmse');              %#ok<AGROW>
        hsh{end+1,1} = ch(T{j,ix('codehash')});%#ok<AGROW>
    end
end
if isempty(sc), fprintf('no rows yet\n'); return; end
Tab = table(string(sc), st, sd, string(me), string(src), cap, it, string(stop), sec, aps, asm, rel, rms, string(hsh), ...
    'VariableNames', {'scenario','setting','seed','method','source','cap','iters','stop_reason','seconds', ...
                      'apsnr_unit','assim','relerr','rmse','codehash'});
writetable(Tab, fullfile(tdir,'all_rows.csv'));
fprintf('rows: %d ; methods: %s\n', height(Tab), strjoin(unique(Tab.method)', ', '));

METH = {'MMC','HOW','L0IQR','TriTD','RTCUR','RGNMR','RMMC'};
SCEN = {'s1','s2','s3','s4sp'};
fid = fopen(fullfile(tdir,'summary.txt'),'w');
for metric = {'apsnr_unit','assim','relerr','seconds'}
    mt = metric{1};
    for c = 1:numel(SCEN)
        s = SCEN{c};
        stv = unique(Tab.setting(Tab.scenario==s))';
        if isempty(stv), continue; end
        vals = nan(numel(METH), numel(stv));
        for m = 1:numel(METH)
            for q = 1:numel(stv)
                sel = Tab.scenario==s & Tab.method==METH{m} & Tab.setting==stv(q);
                if ~any(sel), continue; end
                x = Tab.(mt)(sel);
                if ~isnumeric(x), x = cellfun(@getv, x); end
                vals(m,q) = mean(x, 'omitnan');
            end
        end
        if all(all(isnan(vals))), continue; end
        fprintf(fid,'\n== %s : %s ==\nsetting: %s\n', s, mt, sprintf('%8g ', stv));
        for m = 1:numel(METH)
            if all(isnan(vals(m,:))), continue; end
            fprintf(fid,'%-6s %s\n', METH{m}, sprintf('%9.4f ', vals(m,:)));
        end
        T = array2table(vals, 'VariableNames', matlab.lang.makeValidName(compose('%g', stv)), 'RowNames', METH);
        writetable(T, fullfile(tdir, sprintf('by_method_%s_%s.csv', mt, s)), 'WriteRowNames', true);
    end
end
fprintf(fid,'\n== completeness ==\n');
keys = unique(strcat(string(Tab.scenario),'|',string(Tab.setting),'|',string(Tab.seed)));
miss = 0;
for k = keys'
    sel = strcat(string(Tab.scenario),'|',string(Tab.setting),'|',string(Tab.seed))==k;
    ms = unique(Tab.method(sel));
    want = ["MMC","RMMC","HOW","L0IQR","TriTD","RTCUR"];
    if ~startsWith(k,'s4sp'), want(end+1) = "RGNMR"; end
    lack = setdiff(want, ms);
    if ~isempty(lack)
        fprintf(fid,'MISSING %-18s -> %s\n', k, strjoin(lack, ','));
        miss = miss + 1;
    end
end
fprintf(fid,'cells incomplete: %d / %d  (RMMC excluded by design)\n', miss, numel(keys));
fclose(fid);
fprintf('wrote %s\n', fullfile(tdir,'summary.txt'));
end

% ---------------------------------------------------------------------------
function v = getv(x)
if isnumeric(x)
    v = x;
elseif ischar(x) || isstring(x)
    s = char(x);
    if isempty(s) || strcmpi(s,'NA'), v = NaN; else, v = str2double(s); end
else
    v = NaN;
end
end
function s = ch(x)
if ischar(x), s = x; elseif isstring(x), s = char(x); else, s = ''; end
end
