# /run-tests - SLIDE Test Suite

Builds and runs the SLIDE test suite using CMake and CTest.

## Usage

```
/run-tests [--debug|--release] [--filter=PATTERN]
```

## Steps

1. **Configure build** (if not already configured):
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Debug
   ```

2. **Build test targets**:
   ```bash
   cmake --build build --target unit_tests -j
   ```

3. **Run tests**:
   ```bash
   ctest --test-dir build --output-on-failure -j2
   ```

4. **Report results**:
   - List any failing tests with file:line references
   - Summarize pass/fail counts

## Options

- `--debug`: Use Debug build (default, includes symbols)
- `--release`: Use Release build (faster execution)
- `--filter=PATTERN`: Run only tests matching pattern (e.g., `--filter=Cell_SPM`)

## Test Locations

- Unit tests: `tests/unit/`
- Integration tests: `tests/integration/`

## Common Test Files

| File | Coverage |
|------|----------|
| Cell_SPM_test.cpp | Single Particle Model cell |
| Cell_ECM_test.cpp | Equivalent Circuit Model cell |
| Module_s_test.cpp | Series module |
| Module_p_test.cpp | Parallel module |
| Battery_test.cpp | Complete battery system |
| Cycler_test.cpp | Cycling procedures |
| Procedure_test.cpp | High-level aging procedures |

## Troubleshooting

- **Build fails**: Check CMake version (requires 3.31+) and C++20 compiler
- **Tests timeout**: Increase timeout with `ctest --timeout 300`
- **Specific test fails**: Run with verbose: `ctest -V -R test_name`
