// This executable is a coverage-mapping container, not a test. Coverage.cmake
// whole-archive-links slide_core and deliberately leaves this trivial main
// uninstrumented so every production object is visible to llvm-cov at zero;
// the reporter unions those mappings with separately exported test counts.
int main()
{
  return 0;
}
