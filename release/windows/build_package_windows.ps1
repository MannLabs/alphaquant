# Build the install package for Windows.
# This script must be run from the root of the repository after running build_installer_windows.ps1

# Download AlphaMap data files (FASTA and CSV files)
Write-Host "Starting downloads of FASTA and CSV files..."

# Create directory for AlphaMap data
$alphaMapDataDir = ".\dist_pyinstaller\alphaquant_gui\_internal\alphamap\data"
New-Item -ItemType Directory -Force -Path $alphaMapDataDir | Out-Null

# Get file list from GitHub
try {
    $downloadList = Invoke-RestMethod -Uri "https://api.github.com/repos/MannLabs/alphamap/contents/alphamap/data?ref=main"
} catch {
    Write-Host "Error: Failed to fetch file list from GitHub API"
    exit 1
}

# Download FASTA and CSV files
$fileCount = 0
foreach ($file in $downloadList) {
    if ($file.name -match "\.(fasta|csv)$") {
        $url = $file.download_url
        $filename = $file.name
        Write-Host "Downloading $filename..."

        try {
            Invoke-WebRequest -Uri $url -OutFile "$alphaMapDataDir\$filename"
            Write-Host "Successfully downloaded $filename"
            $fileCount++
        } catch {
            Write-Host "Error: Failed to download $filename"
            exit 1
        }
    }
}

# Verify downloads
Write-Host "Downloaded $fileCount files"
if ($fileCount -eq 0) {
    Write-Host "Error: No files were downloaded"
    exit 1
}

# Wrapping the pyinstaller folder in a .exe package
& "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" .\release\windows\alphaquant_innoinstaller.iss
