$ErrorActionPreference = "Stop"

$ProjectRoot = (Resolve-Path $PSScriptRoot).Path
$ImageName = "plasmodium-docking"

docker build -t $ImageName $ProjectRoot
docker run --rm -v "${ProjectRoot}:/work" $ImageName bash scripts/run_workflow.sh
