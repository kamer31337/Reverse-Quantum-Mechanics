@echo off
setlocal enabledelayedexpansion

echo ================================================================================
echo Building Reverse Quantum Mechanics (RQM) Terminal Program
echo ================================================================================

gcc --version >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo [INFO] Found gcc. Compiling with gcc -std=c11...
    gcc -std=c11 -O3 -Wall -Wextra -Iinclude src/main.c src/rqm_math.c src/rqm_classical.c src/rqm_hilbert_pathology.c src/rqm_ensembles.c src/rqm_born_entropy.c src/rqm_unitary.c src/rqm_equilibration.c src/rqm_classical_limit.c src/rqm_triple_equilibria.c src/rqm_ui.c src/rqm_tests.c -lm -o rqm.exe
    if %ERRORLEVEL% equ 0 (
        echo [SUCCESS] Build completed successfully: rqm.exe
        echo Run 'rqm.exe' for interactive UI or 'rqm.exe --test' for theorem verification.
        exit /b 0
    ) else (
        echo [ERROR] gcc compilation failed.
        exit /b 1
    )
)

clang --version >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo [INFO] Found clang. Compiling with clang -std=c11...
    clang -std=c11 -O3 -Wall -Wextra -Iinclude src/main.c src/rqm_math.c src/rqm_classical.c src/rqm_hilbert_pathology.c src/rqm_ensembles.c src/rqm_born_entropy.c src/rqm_unitary.c src/rqm_equilibration.c src/rqm_classical_limit.c src/rqm_triple_equilibria.c src/rqm_ui.c src/rqm_tests.c -o rqm.exe
    if %ERRORLEVEL% equ 0 (
        echo [SUCCESS] Build completed successfully: rqm.exe
        echo Run 'rqm.exe' for interactive UI or 'rqm.exe --test' for theorem verification.
        exit /b 0
    ) else (
        echo [ERROR] clang compilation failed.
        exit /b 1
    )
)

cl >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo [INFO] Found MSVC cl.exe. Compiling...
    cl /std:c11 /O2 /Iinclude src\main.c src\rqm_math.c src\rqm_classical.c src\rqm_hilbert_pathology.c src\rqm_ensembles.c src\rqm_born_entropy.c src\rqm_unitary.c src\rqm_equilibration.c src\rqm_classical_limit.c src\rqm_triple_equilibria.c src\rqm_ui.c src\rqm_tests.c /Fe:rqm.exe
    if %ERRORLEVEL% equ 0 (
        echo [SUCCESS] Build completed successfully: rqm.exe
        echo Run 'rqm.exe' for interactive UI or 'rqm.exe --test' for theorem verification.
        exit /b 0
    ) else (
        echo [ERROR] cl compilation failed.
        exit /b 1
    )
)

echo [ERROR] No supported C compiler found (gcc, clang, or cl). Please install GCC or Clang.
exit /b 1
