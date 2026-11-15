# Restart guard (J. Picione approved 2026-11-15, after an unexplained stop of every Git Bash lane at 05:36).
# Runs under PowerShell only (no Git Bash), started through WMI like the lanes. Every 10 s it checks that each lane
# still has a process. Conditions set by J. Picione:
#  1. a lane that stops twice within one hour is not restarted again: docs/GUARD_HOLD is written and the guard waits;
#  2. never a second copy: a lane is restarted only when none of its processes runs and no job of it is active in its
#     folder (no held entry in the core ledger with a live process, no orphaned ORCA/MPI process);
#  3. restarts go through the core queue: the restarted lanes start every ORCA job through pipeline/cores.py as before;
#  4. every stop is logged in docs/guard_log.md: time, processes running before and after, Git Bash runtime version;
#  5. it never changes a setting, never deletes anything, never reports by itself. Killed sides (incomplete, not
#     results) are moved to 02_irc/attempts/killed_<time>/ with WHY.txt, because a side starts only in a clean folder.
# Usage (PowerShell): & ([scriptblock]::Create((Get-Content -Raw pipeline\restart_guard.ps1))) [-DryRun] [-Root <dir>]
param([switch]$DryRun, [string]$Root = 'C:\Users\giova\Mechanism Animations', [switch]$SelfTest)
$ErrorActionPreference = 'Continue'
$Log = "$Root\docs\guard_log.md"; $Prod = "$Root\docs\production_true_irc.log"; $Hold = "$Root\docs\GUARD_HOLD"
$Until = Get-Date '2026-11-20 12:00'; $Poll = 10
$Bash = 'C:\Program Files\Git\usr\bin\bash.exe'
$Pre = "export PATH=/usr/bin:/bin:`$PATH; cd '/c/Users/giova/Mechanism Animations' && "
$Chain = @('hshift_2_3_dimethyl_2_butyl_m062x', 'hshift_3_methyl_2_butyl_m062x', 'e2_ethoxide_2_bromopropane_m062x')
$Lanes = @(
  @{ key = 'chain'; name = 'production chain'; match = 'run_true_irc_prod\.sh|resume_true_irc_prod\.sh|pipeline\.true_irc_prod' },
  @{ key = 'watchdog'; name = 'stall watchdog'; match = 'pipeline\.stall_watchdog' })
# Self-test (-SelfTest -Root <scratch folder>): one dummy lane (a Git Bash 'sleep'), 2 s polls, 3 minutes
if ($SelfTest) { $Lanes = @(@{ key = 'test'; name = 'self-test lane'; match = 'guard_selftest_marker' }); $Poll = 2; $Until = (Get-Date).AddMinutes(3) }

function Now { (Get-Date).ToString('yyyy-MM-dd HH:mm:ss') }
function W([string]$t) { if ($DryRun) { $t } else { Add-Content -Path $Log -Value $t -Encoding UTF8 } }
function Snap { $h = @{}; Get-CimInstance Win32_Process | ForEach-Object { $h["$($_.ProcessId)|$($_.CreationDate.Ticks)"] = $_ }; $h }
function Line($p) { $c = "$($p.CommandLine)" -replace '\s+', ' '; '  - {0} pid {1} (parent {2}), started {3}: {4}' -f $p.Name, $p.ProcessId, $p.ParentProcessId, $p.CreationDate.ToString('MM-dd HH:mm:ss'), $c.Substring(0, [Math]::Min(150, $c.Length)) }
function LaneProcs($snap, $lane) { @($snap.Values | Where-Object { $_.Name -match '^(bash|python)\.exe$' -and "$($_.CommandLine)" -match $lane.match -and "$($_.CommandLine)" -notmatch 'shell-snapshots' }) }
function Remaining { $t = if (Test-Path $Prod) { Get-Content $Prod -Raw } else { '' }; @($Chain | Where-Object { $t -notmatch [regex]::Escape("production true IRC ${_}: join exit") }) }
function Runtime { $v = (Get-Item "$(Split-Path $Bash)\msys-2.0.dll").VersionInfo.FileVersion; $g = (& 'C:\Program Files\Git\cmd\git.exe' --version) 2>$null; "msys-2.0.dll $v; $g" }
function Busy($snap, [string]$rx) {
  $why = @(); $pids = @{}; $snap.Values | ForEach-Object { $pids[[int]$_.ProcessId] = 1 }
  if ($rx) { $snap.Values | Where-Object { $_.Name -match '^(orca|mpiexec|smpd)' -and -not $pids.ContainsKey([int]$_.ParentProcessId) } | ForEach-Object { $why += "orphaned $($_.Name) pid $($_.ProcessId)" } }
  try { $q = Get-Content "$Root\docs\level_test_queue\allocation.json" -Raw | ConvertFrom-Json
        foreach ($h in @($q.held)) { if ($rx -and "$($h.folder)" -like "*\reactions\$rx\*" -and $pids.ContainsKey([int]$h.pid)) { $why += "core-queue job $($h.what) (pid $($h.pid))" } } } catch { $why += 'core ledger unreadable' }
  $why }
function Archive([string]$rx, [string]$stamp) {
  $Irc = "$Root\reactions\$rx\02_irc"; $A = "$Irc\attempts\killed_$stamp"; $moved = @()
  foreach ($side in 'backward', 'forward') {
    if ((Test-Path "$Irc\true_irc\$side") -and -not (Test-Path "$Irc\true_irc\$side\path.json")) {
      $n = @(Get-ChildItem "$Irc\true_irc\$side" -Directory -Filter 'g*').Count
      if (-not $DryRun) { New-Item -ItemType Directory -Force $A | Out-Null; Move-Item "$Irc\true_irc\$side" "$A\$side"
        if (Test-Path "$Irc\console_true_irc_$side.txt") { Move-Item "$Irc\console_true_irc_$side.txt" "$A\" } }
      $moved += "$side ($n gradient jobs)" } }
  if ($moved.Count -and -not $DryRun) { Set-Content "$A\WHY.txt" -Encoding UTF8 -Value ("NOT results. Side(s) " + ($moved -join ', ') + " stopped before finishing; found by the restart guard at $(Now) (details in docs/guard_log.md). Rerun from the TS in a clean folder with identical settings.") }
  $moved }
function Start-Lane($lane) {
  if ($lane.key -eq 'chain') { $c = "`"$Bash`" -c `"$Pre bash pipeline/resume_true_irc_prod.sh $((Remaining) -join ' ') >> docs/console_production_chain_guard.txt 2>&1`"" }
  elseif ($lane.key -eq 'test') { $c = "`"$Bash`" -c `"export PATH=/usr/bin:/bin:`$PATH; sleep 600; echo guard_selftest_marker`"" }
  else { $c = "`"$Bash`" -c `"$Pre .venv/Scripts/python.exe -X utf8 -m pipeline.stall_watchdog >> docs/console_stall_watchdog.txt 2>&1`"" }
  if ($DryRun) { return "would start: $c" }
  $r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $c; CurrentDirectory = $Root }
  "return code $($r.ReturnValue), pid $($r.ProcessId)" }

$prev = Snap; $prevT = Now; $alive = @{}; $stops = @{}; $pend = @{}; $held = @{}; $done = @{}; $waitNote = @{}
W "`n## Guard started $(Now) (PowerShell pid $PID; until $($Until.ToString('yyyy-MM-dd HH:mm')); poll every $Poll s)"
W "- Git Bash runtime: $(Runtime)"
foreach ($l in $Lanes) { $n = @(LaneProcs $prev $l).Count; $alive[$l.key] = $n -gt 0; $stops[$l.key] = @(); W "- $($l.name): $n process(es) running" }
W "- production chain, reactions not yet joined: $((Remaining) -join ', ')"
if ($DryRun) { foreach ($rx in Remaining) { "busy check ${rx}: $((Busy $prev $rx) -join '; ')"; "archive check ${rx}: $((Archive $rx 'dryrun') -join ', ')" }; foreach ($l in $Lanes) { Start-Lane $l }; return }

while ((Get-Date) -lt $Until) {
  Start-Sleep -Seconds $Poll
  $cur = Snap; $t = Now
  foreach ($l in $Lanes) {
    $k = $l.key
    if ($held[$k] -or $done[$k]) { continue }
    if ($k -eq 'chain' -and @(Remaining).Count -eq 0) { $done[$k] = 1; W "`n### $t production chain finished: every reaction has a join line. Guard no longer watches it."; continue }
    if (@(LaneProcs $cur $l).Count -gt 0) { $alive[$k] = $true; continue }
    if ($alive[$k]) {
      $alive[$k] = $false; $stops[$k] = @($stops[$k] | Where-Object { ((Get-Date) - $_).TotalMinutes -lt 60 }) + (Get-Date)
      $gone = @($prev.Keys | Where-Object { -not $cur.ContainsKey($_) } | ForEach-Object { $prev[$_] })
      $new = @($cur.Keys | Where-Object { -not $prev.ContainsKey($_) } | ForEach-Object { $cur[$_] })
      W "`n### STOP $($l.name): found at $t; last seen running at $prevT"
      W "- Git Bash runtime: $(Runtime)"
      W "- lane processes at $prevT (now gone):"; $prev.Values | Where-Object { $_.Name -match '^(bash|python)\.exe$' -and "$($_.CommandLine)" -match $l.match } | ForEach-Object { W (Line $_) }
      W "- every process that ended between $prevT and ${t}:"; $gone | Sort-Object Name | ForEach-Object { W (Line $_) }
      W "- every process that started between $prevT and ${t}:"; $new | Sort-Object Name | ForEach-Object { W (Line $_) }
      W "- bash, python, ORCA, MPI, PowerShell, Claude and WMI host processes running at ${t}:"; $cur.Values | Where-Object { $_.Name -match '^(bash|python|orca|mpiexec|smpd|powershell|claude|WmiPrvSE)' } | Sort-Object Name | ForEach-Object { W (Line $_) }
      if ($stops[$k].Count -ge 2) {
        $held[$k] = 1; $m = "$t $($l.name) stopped twice within one hour; not restarted. Waiting for J. Picione."
        Add-Content $Hold $m -Encoding UTF8; W "- HOLD: $m"; Add-Content $Prod "$t restart guard: $($l.name) stopped twice within one hour; not restarted, waiting for J. Picione" -Encoding UTF8; continue }
      $pend[$k] = Get-Date }
    if ($pend[$k]) {
      $rx = if ($k -eq 'chain') { @(Remaining)[0] } else { '' }
      $b = @(Busy $cur $rx)
      if ($b.Count) {
        if (((Get-Date) - $pend[$k]).TotalMinutes -gt 15) { $held[$k] = 1; $m = "$t $($l.name) not restarted: jobs of the stopped lane still active after 15 min ($($b -join '; ')). Waiting for J. Picione."; Add-Content $Hold $m -Encoding UTF8; W "- HOLD: $m" }
        elseif (-not $waitNote[$k]) { $waitNote[$k] = 1; W "- $t waiting before restart: $($b -join '; ')" }
        continue }
      $arch = @(); if ($k -eq 'chain') { $arch = @(Archive $rx ((Get-Date).ToString('yyyyMMdd_HHmm'))) }
      $res = Start-Lane $l
      W "- $t restarted $($l.name) through WMI: $res$(if ($arch.Count) { "; archived as not results in $rx/02_irc/attempts: $($arch -join ', ')" })"
      Add-Content $Prod "$t restart guard: $($l.name) stopped (last seen $prevT) and restarted$(if ($arch.Count) { "; $rx sides archived ($($arch -join ', ')) and rerun from the TS" })" -Encoding UTF8
      if ($res -notmatch 'return code 0') { $held[$k] = 1; Add-Content $Hold "$t $($l.name) restart failed ($res). Waiting for J. Picione." -Encoding UTF8 }
      $pend[$k] = $null; $waitNote[$k] = $null; $alive[$k] = $true }
  }
  $prev = $cur; $prevT = $t
}
W "`n## Guard ended $(Now) (end date reached)"
