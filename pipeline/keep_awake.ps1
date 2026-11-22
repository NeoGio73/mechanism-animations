# Keep-awake for the calculations (re-created 2026-11-22 after the power outage; the 2026-11-17 version was started as
# an encoded command only). While this process runs, Windows does not go to sleep on idle (ES_CONTINUOUS |
# ES_SYSTEM_REQUIRED, a request held by this process only); the display may still turn off. It changes no power,
# system or McAfee setting: the request ends when the process ends. Stops by itself at the same end time as the
# restart guard (Wed 2026-12-02 12:00).
# Started hidden through WMI: powershell.exe -NoProfile -File pipeline\keep_awake.ps1
$sig = '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
$k = Add-Type -MemberDefinition $sig -Name Pwr -Namespace KeepAwake -PassThru
$stopAt = Get-Date '2026-12-02 12:00'
while ((Get-Date) -lt $stopAt) {
  [void]$k::SetThreadExecutionState([uint32]2147483649)   # 0x80000001
  Start-Sleep -Seconds 60
}
[void]$k::SetThreadExecutionState([uint32]2147483648)       # ES_CONTINUOUS alone: release the request
