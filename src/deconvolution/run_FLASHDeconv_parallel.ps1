<#
  run_FLASHDeconv_parallel.ps1
  Deconvolve every mzML in a folder, ONE thread per file, several files at once.

  Why: single-threaded per file (-threads 1) avoids the intra-file data race
  that produced the heap-corruption crashes (exit 0xC0000374); running many
  files concurrently keeps all cores busy. Each file is its own OS process.

  Idempotent / resumable: any file whose *_ms1.tsv already exists (non-empty)
  is skipped, so you can re-run this to pick up only the failed/missing files.
  Per-file stdout/stderr go to <OutDir>\logs so parallel output isn't interleaved.
#>

param(
    [string]$FLASHDeconv = 'C:\Program Files\OpenMS-3.1.0-pre-FDRdeploy\bin\FLASHDeconv.exe',
    [string]$mzmlDir     = 'F:\SP3_TDMS\mzML',
    [string]$OutDir      = 'F:\SP3_TDMS\flashdeconv',
    [int]   $MaxParallel = [Environment]::ProcessorCount   # files run at once; lower for RAM/CPU headroom
)

$logDir = Join-Path $OutDir 'logs'
New-Item -ItemType Directory -Force -Path $OutDir, $logDir | Out-Null

# Everything except -in / -out / -out_spec1. Note -threads 1.
$commonArgs = @(
    '-write_detail',
    '-FD:report_FDR',
    '-FD:allowed_isotope_error', '0',
    '-SD:tol', '10',
    '-SD:min_mass', '2000', '-SD:max_mass', '70000',
    '-SD:min_charge', '6', '-SD:max_charge', '60',
    '-threads', '1'
)

function Get-Running { param($list) @($list | Where-Object { -not $_.Proc.HasExited }) }

$files   = Get-ChildItem -Path $mzmlDir -Filter *.mzML | Sort-Object Name
$launched = @()
$skipped  = 0

Write-Host "Launching up to $MaxParallel file(s) at a time, 1 thread each.`n"

foreach ($file in $files) {
    $base = $file.BaseName
    $feat = Join-Path $OutDir "$base.tsv"
    $ms1  = Join-Path $OutDir "${base}_ms1.tsv"

    if ((Test-Path $ms1) -and (Get-Item $ms1).Length -gt 0) {
        Write-Host "SKIP  $base (already done)"
        $skipped++
        continue
    }

    # Throttle: wait until a slot frees up.
    while ((Get-Running $launched).Count -ge $MaxParallel) {
        Start-Sleep -Milliseconds 500
    }

    $fdArgs = @('-in', "`"$($file.FullName)`"", '-out', "`"$feat`"", '-out_spec1', "`"$ms1`"") + $commonArgs
    $outLog = Join-Path $logDir "$base.out.log"
    $errLog = Join-Path $logDir "$base.err.log"

    $proc = Start-Process -FilePath $FLASHDeconv -ArgumentList $fdArgs -NoNewWindow -PassThru `
                          -RedirectStandardOutput $outLog -RedirectStandardError $errLog
    $launched += [pscustomobject]@{ Name = $base; Proc = $proc; Ms1 = $ms1 }
    Write-Host ("START {0}  ({1} running)" -f $base, (Get-Running $launched).Count)
}

# Wait for the last batch to finish.
while ((Get-Running $launched).Count -gt 0) {
    Start-Sleep -Milliseconds 500
}

# Summary.
Write-Host "`n===== Summary ====="
$failed = @()
foreach ($job in $launched) {
    $code = $job.Proc.ExitCode
    $ok   = ($code -eq 0) -and (Test-Path $job.Ms1) -and ((Get-Item $job.Ms1).Length -gt 0)
    if ($ok) {
        Write-Host ("OK    {0}" -f $job.Name)
    } else {
        $hex  = ('0x{0:X8}' -f ($code -band 0xFFFFFFFF))
        $note = if ($hex -eq '0xC0000374') { ' = heap corruption' } else { '' }
        Write-Warning ("FAIL  {0}  exit {1} ({2}{3})" -f $job.Name, $code, $hex, $note)
        $failed += $job.Name
    }
}

Write-Host ("`nSkipped {0} (already done), ran {1}, failed {2}." -f $skipped, $launched.Count, $failed.Count)
if ($failed.Count) {
    Write-Host "Still failing: $($failed -join ', ')"
    Write-Host "Re-run this script to retry just those (finished files are skipped). Logs: $logDir"
} else {
    Write-Host "All files completed."
}