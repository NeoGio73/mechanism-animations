# Start a long-running command outside the Claude app's process tree (2026-11-15: an automatic update of the
# app at 00:45 ended every process it had started, including the production paths). The process is created by
# Windows' WMI service, so it is not a child of the app and survives an app update or restart. Nothing is
# installed and no setting is changed.
# Usage: powershell -NoProfile -File pipeline/start_detached.ps1 -Command "<command line>" [-Dir "<working folder>"]
param([Parameter(Mandatory = $true)][string]$Command,
      [string]$Dir = "C:\Users\giova\Mechanism Animations")
$r = Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{ CommandLine = $Command; CurrentDirectory = $Dir }
if ($r.ReturnValue -ne 0) { throw "Win32_Process.Create failed with code $($r.ReturnValue)" }
"started pid $($r.ProcessId): $Command"
