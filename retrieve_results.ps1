# =============================================================================
# retrieve_results.ps1
# Copies all results from HPC back to local machine
# Usage: .\retrieve_results.ps1
# =============================================================================

$HPC_USER = "sequinlan"
$HPC_HOST = "sonic.ucd.ie"
$HPC_PATH = "/home/people/sequinlan/GeoMx_Analysis/results"
$LOCAL_PATH = "$PSScriptRoot\results"

# Create local results directory if it doesn't exist
New-Item -ItemType Directory -Force -Path $LOCAL_PATH | Out-Null

Write-Host "Copying results from HPC..." -ForegroundColor Cyan

scp -r "${HPC_USER}@${HPC_HOST}:${HPC_PATH}/" "$LOCAL_PATH"

if ($LASTEXITCODE -eq 0) {
    Write-Host "`nDone! Files copied to: $LOCAL_PATH" -ForegroundColor Green
    Write-Host "Launch dashboard with: source('launch_dashboard.R')" -ForegroundColor Green
    Get-ChildItem $LOCAL_PATH | Format-Table Name, @{N="Size_MB";E={[math]::Round($_.Length/1MB,2)}} -AutoSize
} else {
    Write-Host "scp failed. Check your HPC connection." -ForegroundColor Red
}
