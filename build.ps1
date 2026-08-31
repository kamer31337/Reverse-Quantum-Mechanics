Write-Host "================================================================================" -ForegroundColor Cyan
Write-Host " Building Reverse Quantum Mechanics (RQM) Terminal Program" -ForegroundColor Green
Write-Host "================================================================================" -ForegroundColor Cyan

$sources = @(
    "src/main.c",
    "src/rqm_math.c",
    "src/rqm_classical.c",
    "src/rqm_hilbert_pathology.c",
    "src/rqm_ensembles.c",
    "src/rqm_born_entropy.c",
    "src/rqm_unitary.c",
    "src/rqm_equilibration.c",
    "src/rqm_classical_limit.c",
    "src/rqm_triple_equilibria.c",
    "src/rqm_ui.c",
    "src/rqm_tests.c"
)

if (Get-Command clang -ErrorAction SilentlyContinue) {
    Write-Host "[INFO] Compiling with clang..." -ForegroundColor Yellow
    clang -std=c11 -O3 -Wall -Wextra -Iinclude $sources -o rqm.exe
    if ($LASTEXITCODE -eq 0) {
        Write-Host "[SUCCESS] Built rqm.exe successfully!" -ForegroundColor Green
        exit 0
    }
}

if (Get-Command gcc -ErrorAction SilentlyContinue) {
    Write-Host "[INFO] Compiling with gcc..." -ForegroundColor Yellow
    gcc -std=c11 -O3 -Wall -Wextra -Iinclude $sources -lm -o rqm.exe
    if ($LASTEXITCODE -eq 0) {
        Write-Host "[SUCCESS] Built rqm.exe successfully!" -ForegroundColor Green
        exit 0
    }
}

Write-Error "No compatible C compiler found (clang or gcc)."
exit 1
