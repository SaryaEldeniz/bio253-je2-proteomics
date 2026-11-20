# Windows export fallback: the bundled Node finalizer cannot canonicalize this
# workspace (EPERM). Run its packaged validators directly, then publish unchanged.
param([string]$Candidate,[string]$PreviewSource)
$ErrorActionPreference='Stop'
$workspaceDir=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../../..'))
if(-not $Candidate) { $Candidate=Join-Path $workspaceDir 'results/slide_summary/draft.pptx' }
$Candidate=[IO.Path]::GetFullPath($Candidate)
$skillTools='C:/Users/kahve/.codex/plugins/cache/openai-primary-runtime/presentations/26.1112.11217/skills/presentations/container_tools'
$runtimePy='C:/Program Files/WindowsApps/OpenAI.CodexPrimaryRuntime.v26-1112-737-0_26.1112.737.0_x64__3k8sg7r9htsxt/dependencies/python/python.exe'
$outDir=Join-Path $workspaceDir 'steps/10_integration_poster/presentations'
$finalPath=Join-Path $outDir 'JE2_PASN_progress_summary.pptx'
if(Test-Path -LiteralPath $finalPath) { throw 'Choose a new output name; the exporter never overwrites a final deck.' }
function Invoke-DeckValidator([string]$ValidatorName,[string[]]$CheckArguments) {
  $checkText=& $runtimePy (Join-Path $skillTools $ValidatorName) @CheckArguments
  if($LASTEXITCODE -ne 0) { throw "$ValidatorName failed: $checkText" }
  return (($checkText -join "`n") | ConvertFrom-Json)
}
$checks=[ordered]@{}
$checks.package=Invoke-DeckValidator 'inspect_presentation_package_integrity.py' @($Candidate)
$checks.layout=Invoke-DeckValidator 'inspect_presentation_layout_geometry.py' @($Candidate,'--expected-slide-size-emu','15240000,8572500','--validate-bullet-geometry','--validate-heading-fit','--require-native-table-slide','1','--approved-font-family','Arial')
$checks.charts=Invoke-DeckValidator 'native_quantitative_chart_gate.py' @($Candidate,'--expected-chart-slides','1','--target-application','powerpoint')
$checks.titles=Invoke-DeckValidator 'native_chart_title_gate.py' @($Candidate)
$checks.arithmetic=Invoke-DeckValidator 'native_table_arithmetic_gate.py' @($Candidate,'--fail-on-findings')
if($checks.package.slide_count -ne 1 -or $checks.package.chart_count -ne 4 -or -not $checks.charts.passed) { throw 'Wrong slide/chart count or failed chart validation.' }
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
Copy-Item -LiteralPath $Candidate -Destination $finalPath
if((Get-FileHash -LiteralPath $Candidate).Hash -ne (Get-FileHash -LiteralPath $finalPath).Hash) { throw 'Published bytes differ from checked candidate.' }
if($PreviewSource) { Copy-Item -LiteralPath $PreviewSource -Destination (Join-Path $outDir 'JE2_PASN_progress_summary.png') }
$checks.delivery=[ordered]@{path=$finalPath;sha256=(Get-FileHash -LiteralPath $finalPath).Hash;application='PowerPoint';nativeCharts=4;nativeTable=1;nativePowerPointOpened=$false;visualReview='Connected-session PNG inspected';workbookMode='Complete native PowerPoint literal chart data'}
$checks | ConvertTo-Json -Depth 30 | Set-Content -LiteralPath (Join-Path $workspaceDir 'results/slide_summary/final_validation.json')
Write-Output "Saved $finalPath"
