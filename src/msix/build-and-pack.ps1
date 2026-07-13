# build-and-pack.ps1 — Build the CLI (Release) and create an MSIX package.
#
# Usage: .\build-and-pack.ps1
#
# Requires: CMake with Visual Studio 2026, WinApp CLI (winapp).
# Run from: src/msix/

$ErrorActionPreference = "Stop"

$ScriptDir  = Split-Path -Parent $MyInvocation.MyCommand.Path
$SrcDir     = Resolve-Path "$ScriptDir\.."
$BuildDir   = "$SrcDir\b\vs-rel"
$MsixDir    = $ScriptDir

# 1. Configure (skip if already configured) -----------------------------------
if (-not (Test-Path "$BuildDir\CMakeCache.txt")) {
    Write-Host "[1/4] Configuring CMake (vs-rel)..." -ForegroundColor Cyan
    cmake --preset vs-rel -S "$SrcDir"
    if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }
}
else {
    Write-Host "[1/4] CMake already configured, skipping." -ForegroundColor DarkGray
}

# 2. Build --------------------------------------------------------------------
Write-Host "[2/4] Building Release..." -ForegroundColor Cyan
cmake --build "$BuildDir" --config Release
if ($LASTEXITCODE -ne 0) { throw "Build failed" }

# 3. Copy cli.exe into the msix folder ----------------------------------------
$CliSource = "$BuildDir\cli\Release\cli.exe"
if (-not (Test-Path $CliSource)) {
    Write-Error "cli.exe not found at $CliSource"
    exit 1
}
Write-Host "[3/4] Copying cli.exe -> $MsixDir\" -ForegroundColor Cyan
Copy-Item $CliSource "$MsixDir\cli.exe" -Force

# 4. Pack MSIX ----------------------------------------------------------------
Write-Host "[4/4] Packing MSIX..." -ForegroundColor Cyan
Push-Location $MsixDir
try {
    winapp pack .
    if ($LASTEXITCODE -ne 0) { throw "winapp pack failed" }
}
finally {
    Pop-Location
}

Write-Host "Done! MSIX package created in $MsixDir" -ForegroundColor Green
