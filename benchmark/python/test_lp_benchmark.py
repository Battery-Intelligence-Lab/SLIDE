# This file is to run liionpack benchmarks given by their team.

from lp_benchmark import *
import time

# ------------------- Basic Benchmark  ----------------------------
benchmark_name = "Basic benchmark"
benchmark_1 = BasicBenchmark()
benchmark_1.setup()
start_time = time.time()
benchmark_1.time_solve_model()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

# ------------------- Small Pack ----------------------------
benchmark_name = "Small Pack"
benchmark = SmallPack()
benchmark.setup()
start_time = time.time()
benchmark.time_discharge_1cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

start_time = time.time()
benchmark.time_discharge_2cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

# ------------------- Medium Pack ----------------------------
benchmark_name = "Medium Pack"
benchmark = MediumPack()
benchmark.setup()
start_time = time.time()
benchmark.time_discharge_1cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

start_time = time.time()
benchmark.time_discharge_2cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

# ------------------- Large Pack ----------------------------
benchmark_name = "Large Pack"
benchmark = LargePack()
benchmark.setup()
start_time = time.time()
benchmark.time_discharge_1cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

start_time = time.time()
benchmark.time_discharge_2cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")

start_time = time.time()
benchmark.time_long_cycle_2cpu()
end_time = time.time()
print(f"{benchmark_name} elapsed time: {(end_time - start_time):.3f} seconds")
