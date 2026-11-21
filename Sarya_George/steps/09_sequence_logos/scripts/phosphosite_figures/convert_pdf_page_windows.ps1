# Format conversion only. Some Windows Inkscape builds keep running after export;
# accept only a complete, parsed SVG, then close this helper's own hidden process.
param(
  [Parameter(Mandatory=$true)][string]$Inkscape,
  [Parameter(Mandatory=$true)][string]$InputPdf,
  [Parameter(Mandatory=$true)][string]$OutputSvg,
  [Parameter(Mandatory=$true)][int]$Page,
  [Parameter(Mandatory=$true)][string]$LogPath
)
$ErrorActionPreference = 'Stop'
$Inkscape = [IO.Path]::GetFullPath($Inkscape)
$InputPdf = [IO.Path]::GetFullPath($InputPdf)
$OutputSvg = [IO.Path]::GetFullPath($OutputSvg)
$LogPath = [IO.Path]::GetFullPath($LogPath)
$temporarySvg = Join-Path ([IO.Path]::GetDirectoryName($LogPath)) ([guid]::NewGuid().ToString() + '.svg')
$arguments = @(('"'+$InputPdf+'"'), ('--pages='+$Page), '--pdf-poppler',
  '--export-type=svg', '--export-plain-svg', '--export-text-to-path',
  ('--export-filename="'+$temporarySvg+'"'))
$process = Start-Process -FilePath $Inkscape -ArgumentList $arguments -WindowStyle Hidden -PassThru `
  -RedirectStandardOutput ($LogPath + '.stdout') -RedirectStandardError ($LogPath + '.stderr')
try {
  $watch = [Diagnostics.Stopwatch]::StartNew()
  $complete = $false
  while ($watch.Elapsed.TotalSeconds -lt 45) {
    if (Test-Path -LiteralPath $temporarySvg) {
      try {
        $document = [xml][IO.File]::ReadAllText($temporarySvg)
        $complete = $document.DocumentElement.LocalName -eq 'svg' -and
          $document.DocumentElement.HasAttribute('viewBox')
      } catch { $complete = $false }
      if ($complete) { break }
    }
    if ($process.HasExited) {
      # The SVG writer may release its file lock just as the process exits.
      Start-Sleep -Milliseconds 250
      try {
        $document = [xml][IO.File]::ReadAllText($temporarySvg)
        $complete = $document.DocumentElement.LocalName -eq 'svg' -and
          $document.DocumentElement.HasAttribute('viewBox')
      } catch { $complete = $false }
      if ($complete) { break }
      throw 'Inkscape exited before producing a complete SVG.'
    }
    Start-Sleep -Milliseconds 150
  }
  if (-not $complete) { throw 'Inkscape did not produce a complete SVG within 45 seconds.' }
  if (-not $process.WaitForExit(1000)) {
    $process.Kill()
    $process.WaitForExit()
  }
  Copy-Item -LiteralPath $temporarySvg -Destination $OutputSvg -Force
  [IO.File]::WriteAllText($LogPath, 'Export verified as complete XML; source PDF read only. Process cleanup limited to the helper launched for this page.')
} catch {
  [IO.File]::WriteAllText($LogPath + '.failed', $_.ToString())
  throw
} finally {
  if (-not $process.HasExited) { $process.Kill() }
  $process.Dispose()
}
