# /benchmark - SLIDE Performance Testing

Executes performance benchmarks and records results for comparison.

## Usage

```
/benchmark [--target=NAME] [--record]
```

## Available Benchmarks

| Target | Description | Expected Time |
|--------|-------------|---------------|
| benchmark_Cell_SPM | Single SPM cell cycling | ~seconds |
| benchmark_Cell_ECM | ECM cell cycling | ~seconds |
| benchmark_Cell_Bucket | Simple bucket cell | ~milliseconds |
| benchmark_LP_cases | Large pack simulations | ~minutes |

## Steps

1. **Build Release** (for accurate timing):
   ```bash
   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --target benchmark_Cell_SPM -j
   ```

2. **Run benchmark**:
   ```bash
   ./build/benchmark/benchmark_Cell_SPM
   ```

3. **Record results** (if `--record`):
   - Append to `benchmark/README.md` with timestamp
   - Format: `| Date | Test | Time | Hardware | Notes |`

## Baseline Performance (v3.0.0)

From documentation:
- 5000 1C CC cycles: < 1 minute
- CC + CV cycles: < 2 minutes
- EPFL battery 1hr CC: ~2 seconds
- EPFL battery + aging: ~3.5 seconds

## Performance Comparison

When comparing builds:
1. Use same hardware/environment
2. Run multiple times, report median
3. Disable other processes for consistent results
4. Compare Debug vs Release impact

## PyBaMM Comparison

For validation against PyBaMM:
```bash
cd benchmark/python
python test_single_cell_SPM.py
```

## Recording Template

Add to `benchmark/README.md`:
```markdown
| 2026-01-29 | Cell_SPM 5000 cycles | 45s | Intel i7-12700 | Release build |
```
