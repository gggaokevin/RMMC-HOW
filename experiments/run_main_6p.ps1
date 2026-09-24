# run_main_6p.ps1 -- run the missing S1-S4 cells with SIX workers, one per PHYSICAL P-core.
#
# Every worker is pinned to both threads of one P-core (masks 0x3,0xC,0x30,0xC0,
# 0x300,0xC00) and runs single-threaded with its own prefdir and TMP, so all cells
# see the same uncontended-core condition.  Rounds repeat until nothing is missing
# (start-up failures of MATLAB are retried), then main_report aggregates.
#
#   powershell -NoProfile -ExecutionPolicy Bypass -File experiments\run_main_6p.ps1
#
#   -MatlabExe   full path to matlab.exe (default: whatever 'matlab' resolves to)
#   -Rounds      maximum planning/run rounds (default 5)
#
# Layout it assumes (see code/README.md):
#   <pkg>/experiments   these programs          <- $PSScriptRoot
#   <pkg>/results       schedule.csv, cells/, logs/, tables/
param(
    [string]$MatlabExe = 'matlab',
    [int]$Rounds     = 5,
    [int]$StaggerSec = 8
)

$ErrorActionPreference = 'Continue'
$pkg     = Split-Path -Parent $PSScriptRoot                       # <pkg>
$code    = Join-Path $pkg 'experiments'
$main    = Join-Path $pkg 'results'
$log     = Join-Path $main 'logs'
$scratch = Join-Path $env:TEMP 'rmmc_how'
New-Item -ItemType Directory -Force -Path $log, $scratch | Out-Null
$masks = @(3, 12, 48, 192, 768, 3072)      # 0x3 .. 0xC00

for ($r = 1; $r -le $Rounds; $r++) {
    Write-Host "===== round $r : planning (6 P-cores) ====="
    $env:MATLAB_PREFDIR = Join-Path $scratch 'prefs_plan'
    $env:TMP = $env:TEMP = Join-Path $scratch 'tmp_plan'
    New-Item -ItemType Directory -Force -Path $env:MATLAB_PREFDIR, $env:TMP | Out-Null
    & $MatlabExe -singleCompThread -batch "cd('$code'); main_plan(6)" *>> (Join-Path $log "plan6_r$r.txt")
    if ($LASTEXITCODE -ne 0) { Write-Host 'plan failed'; break }

    $sch = Import-Csv (Join-Path $main 'schedule.csv')
    Write-Host "round $r : $($sch.Count) cells remaining"
    if ($sch.Count -eq 0) { Write-Host 'nothing left'; break }

    $procs = @()
    for ($w = 1; $w -le 6; $w++) {
        $mine = $sch | Where-Object { [int]$_.worker -eq $w }
        if (-not $mine) { continue }
        $env:MATLAB_PREFDIR = Join-Path $scratch "prefs_p$w"
        $env:TMP = $env:TEMP = Join-Path $scratch "tmp_p$w"
        New-Item -ItemType Directory -Force -Path $env:MATLAB_PREFDIR, $env:TMP | Out-Null
        $cmd = "cd('$code'); main_worker($w)"
        $p = Start-Process -FilePath $MatlabExe -ArgumentList @('-singleCompThread','-batch',"`"$cmd`"") -PassThru -NoNewWindow `
             -RedirectStandardOutput (Join-Path $log "p$w.out") -RedirectStandardError (Join-Path $log "p$w.err")
        try { $p.ProcessorAffinity = [IntPtr]$masks[$w-1] } catch { }
        $procs += $p
        Start-Sleep -Seconds $StaggerSec
    }
    Write-Host "round $r : launched $($procs.Count) workers"
    $procs | Wait-Process
    Write-Host "round $r : workers exited"
}

Write-Host '===== aggregation ====='
$env:MATLAB_PREFDIR = Join-Path $scratch 'prefs_main'
$env:TMP = $env:TEMP = Join-Path $scratch 'tmp_main'
& $MatlabExe -singleCompThread -batch "cd('$code'); main_report" *>> (Join-Path $log 'report6.txt')
Write-Host 'done'
