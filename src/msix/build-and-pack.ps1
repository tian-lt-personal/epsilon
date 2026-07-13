# build-and-pack.ps1 — Build the CLI (Release) and create an MSIX package.
#
# Usage: .\build-and-pack.ps1 [-Version x.y.z] [-Revision N]
#
# Requires: CMake with Visual Studio 2026, WinApp CLI (winapp).
# Run from: src/msix/

param(
    [string]$Version,
    [int]$Revision = 0
)

$ErrorActionPreference = "Stop"

$ScriptDir  = Split-Path -Parent $MyInvocation.MyCommand.Path
$SrcDir     = Resolve-Path "$ScriptDir\.."
$BuildDir   = "$SrcDir\b\vs-rel"
$MsixDir    = $ScriptDir
$VersionFile = "$SrcDir\VERSION"

# 0. Resolve version -----------------------------------------------------------
if (-not $Version) {
    if (-not (Test-Path $VersionFile)) {
        Write-Error "VERSION file not found at $VersionFile"
        exit 1
    }
    $Version = (Get-Content $VersionFile).Trim()
}

$MsixVersion = "$Version.$Revision"
Write-Host "Version: $MsixVersion" -ForegroundColor Cyan

# 1. Configure (skip if already configured) -----------------------------------
if (-not (Test-Path "$BuildDir\CMakeCache.txt")) {
    Write-Host "[1/5] Configuring CMake (vs-rel)..." -ForegroundColor Cyan
    cmake --preset vs-rel -S "$SrcDir"
    if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }
}
else {
    Write-Host "[1/5] CMake already configured, skipping." -ForegroundColor DarkGray
}

# 2. Build --------------------------------------------------------------------
Write-Host "[2/5] Building Release..." -ForegroundColor Cyan
cmake --build "$BuildDir" --config Release
if ($LASTEXITCODE -ne 0) { throw "Build failed" }

# 3. Copy cli.exe into the msix folder ----------------------------------------
$CliSource = "$BuildDir\cli\Release\cli.exe"
if (-not (Test-Path $CliSource)) {
    Write-Error "cli.exe not found at $CliSource"
    exit 1
}
Write-Host "[3/5] Copying cli.exe -> $MsixDir\" -ForegroundColor Cyan
Copy-Item $CliSource "$MsixDir\cli.exe" -Force

# 4. Inject version into manifest and pack ------------------------------------
Write-Host "[4/5] Injecting version into manifest..." -ForegroundColor Cyan
$ManifestPath = "$MsixDir\Package.appxmanifest"
$ManifestBackup = "$MsixDir\Package.appxmanifest.bak"

Copy-Item $ManifestPath $ManifestBackup
try {
    [xml]$manifest = Get-Content $ManifestPath
    $manifest.Package.Identity.Version = $MsixVersion
    $manifest.Save($ManifestPath)

    Write-Host "[5/5] Packing MSIX..." -ForegroundColor Cyan
    Push-Location $MsixDir
    try {
        winapp pack .
        if ($LASTEXITCODE -ne 0) { throw "winapp pack failed" }
    }
    finally {
        Pop-Location
    }
}
finally {
    Copy-Item $ManifestBackup $ManifestPath -Force
    Remove-Item $ManifestBackup -Force
}

Write-Host "Done! MSIX package created in $MsixDir" -ForegroundColor Green
