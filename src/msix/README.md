# MSIX Packaging

The Epsilon CLI is distributed as an [MSIX](https://learn.microsoft.com/en-us/windows/msix/) package for Windows (x64).

## Versioning

The package version is derived from [`../VERSION`](../VERSION) — the single source of truth for the project version. The 3-part semver (`0.0.1`) is extended to a 4-part MSIX version (`0.0.1.0`) by appending a revision number.

To override the version at build time, pass `-Version` to the build script.

## Build Locally

**Prerequisites:**
- Windows with Visual Studio 2026 (or later) and the MSVC toolchain
- CMake 3.20+
- [WinApp CLI](https://learn.microsoft.com/en-us/windows/msix/desktop/desktop-to-uwp-manual-conversion) (`winapp`)

**Steps:**

```powershell
# From src/msix/
cd src/msix

# Build with the version from ../VERSION
.\build-and-pack.ps1

# Or override the version
.\build-and-pack.ps1 -Version "0.2.0"

# Or override both version and revision (→ 0.2.0.1)
.\build-and-pack.ps1 -Version "0.2.0" -Revision 1
```

The script:

1. Reads the version from `../VERSION` (or the `-Version` parameter)
2. Configures CMake (`vs-rel` preset)
3. Builds `cli.exe` in Release
4. Copies `cli.exe` into the MSIX directory
5. Injects the version into `Package.appxmanifest`
6. Runs `winapp pack .` to create the `.msix`
7. Restores `Package.appxmanifest` to its placeholder state

The resulting package is written to the current directory as `OrchidApps.Calcli_<version>_x64.msix`.

**Install locally for testing:**

```powershell
# Install (requires the package to be signed or installed with dev trust)
Add-AppxPackage -Path .\OrchidApps.Calcli_0.0.1.0_x64.msix
```

## Publish via Git Tags

Releases are automated through GitHub Actions. Pushing a `v*` tag triggers the release workflow.

**Steps:**

```bash
# 1. Bump the version
echo -n "0.1.0" > src/VERSION

# 2. Commit
git add src/VERSION
git commit -m "Bump version to 0.1.0"

# 3. Tag and push
git tag v0.1.0
git push origin main
git push origin v0.1.0
```

The CI workflow (`.github/workflows/epsilon-ci.yml`) will:

1. Validate that the tag (`v0.1.0`) matches the VERSION file (`0.1.0`)
2. Build the CLI in Release
3. Build the MSIX package
4. Upload the MSIX as a workflow artifact
5. Create a GitHub Release with the MSIX attached

The release job only runs after the main CI build (GCC, Clang, MSVC) and tests pass.

## Manifest

`Package.appxmanifest` is checked in with `Version="0.0.0.0"` as a placeholder. The actual version is injected by `build-and-pack.ps1` at build time and restored afterwards — the placeholder never leaves the working tree.
