# Build the install package for Windows.
# This script must be run from the root of the repository after running build_installer_windows.ps1

# A directory for handling the AlphaMap (no typo!) data, analogous to the Mac installer, is added via the .inno file

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

# Download and extract AlphaQuant resources from datashare
Write-Host "Creating resources directory..."
$resourcesDir = ".\dist_pyinstaller\alphaquant_gui\_internal\alphaquant\resources"
New-Item -ItemType Directory -Force -Path $resourcesDir | Out-Null

# Download and extract the first zip file
Write-Host "Downloading and extracting first resource from datashare..."
$tempZip1 = [System.IO.Path]::GetTempFileName()
try {
    Invoke-WebRequest -Uri "https://datashare.biochem.mpg.de/s/ezPzeqStEgDD8gg/download" -OutFile $tempZip1
    Write-Host "Successfully downloaded first resource, now extracting..."
    Expand-Archive -Path $tempZip1 -DestinationPath $resourcesDir -Force
    Remove-Item -Path $tempZip1
} catch {
    Write-Host "Error: Failed to download or extract first resource from datashare"
    exit 1
}

# Download and extract the second zip file
Write-Host "Downloading and extracting second resource from datashare..."
$tempZip2 = [System.IO.Path]::GetTempFileName()
try {
    Invoke-WebRequest -Uri "https://datashare.biochem.mpg.de/s/stH9pmNe6O9CRHG/download" -OutFile $tempZip2
    Write-Host "Successfully downloaded second resource, now extracting..."
    Expand-Archive -Path $tempZip2 -DestinationPath $resourcesDir -Force
    Remove-Item -Path $tempZip2
} catch {
    Write-Host "Error: Failed to download or extract second resource from datashare"
    exit 1
}

# Verify the downloads and extraction
Write-Host "Verifying extracted resources..."
if (-not (Test-Path -Path $resourcesDir)) {
    Write-Host "Error: Resources directory not found after extraction"
    exit 1
}
Write-Host "Resources successfully downloaded and extracted"

# Wrapping the pyinstaller folder in a .exe package
&  "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" .\release\windows\alphaquant_innoinstaller.iss
