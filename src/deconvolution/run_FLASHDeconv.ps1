# --- paths ---
$FLASHDeconv = 'C:\Program Files\OpenMS-3.1.0-pre-FDRdeploy\bin\FLASHDeconv.exe'
$mzmlDir     = 'F:\SP3_TDMS\mzML'
$OutDir      = 'F:\SP3_TDMS\flashdeconv'

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

# --- stage 2: deconvolve every mzML with the FDR module ---
Get-ChildItem -Path $mzmlDir -Filter *.mzML | ForEach-Object {
    $base = $_.BaseName
    $feat = Join-Path $OutDir ($base + '.tsv')
    $ms1  = Join-Path $OutDir ($base + '_ms1.tsv')
    Write-Host "Deconvolving $base ..."
    & $FLASHDeconv -in $_.FullName -out $feat -out_spec1 $ms1 `
        -write_detail -FD:report_FDR -FD:allowed_isotope_error 0 `
        -SD:tol 10 -SD:min_mass 2000 -SD:max_mass 70000 `
        -SD:min_charge 6 -SD:max_charge 60 -threads 10
    if ($LASTEXITCODE -ne 0) { Write-Warning "FLASHDeconv failed on $base (exit $LASTEXITCODE)" }
}