# FEM Electrostatics

A minimal finite-element electrostatics solver written in C++17.

---

## 🧱 Prerequisites
- **CMake ≥ 3.16**
- **A C++17 compiler**, e.g.  
  - MSVC 2022 (`cl.exe`)   → for NMake / Visual Studio / Ninja  
  - MinGW-w64 (g++)        → for Ninja  
  - Clang / AppleClang     → for macOS / Linux  
  - Intel oneAPI (icx / icpx)

---

## ⚙️ Build on Windows (using MSVC + Ninja)
```bat
:: open Visual Studio Developer Command Prompt
cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

Executable → `build\bin\Debug\fem_solver.exe`

---

## ⚙️ Build on Windows (using MSVC + NMake)
```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat" -arch=amd64
cmake -G "NMake Makefiles" -S . -B build-nmake -DCMAKE_BUILD_TYPE=Debug
cmake --build build-nmake
```

Executable → `build-nmake\bin\Debug\fem_solver.exe`

---

## ⚙️ Build on Linux / macOS
```bash
cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/bin/Release/fem_solver
```

---

## 🧭 Run & Debug (in Cursor or VS Code)
1. Choose **Run and Debug → Debug (MSVC, Ninja/NMake build)** from the drop-down.  
2. Set breakpoints and press **F5**.  
3. The debugger uses `${workspaceFolder}/build/bin/<Config>/fem_solver.exe`.

---

## 🧰 Other Targets
| Target | Purpose |
|---------|----------|
| `fem_solver` | main executable |
| `build_instructions` | prints quick build/run help |
| `clean` | remove intermediates |

---

## 📄 License
MIT / internal research use — see LICENSE if provided.
