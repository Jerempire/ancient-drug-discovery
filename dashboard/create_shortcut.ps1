$ws = New-Object -ComObject WScript.Shell
$sc = $ws.CreateShortcut("$env:USERPROFILE\Desktop\Binder Explorer.lnk")
$sc.TargetPath = "$env:USERPROFILE\Projects\medical\ancient-drug-discovery\dashboard\index.html"
$sc.Description = "Ancient Drug Discovery - Binder Explorer Dashboard"
$sc.Save()
Write-Host "Shortcut created on Desktop."
