"""The gate's fixed vocabulary: failure classes, source suffixes, and schema versions."""

FAILURE_CLASSES = {
    "defensive-only",
    "platform-specific",
    "unreachable-by-construction",
}
SOURCE_SUFFIXES = {".cc", ".cpp", ".cu", ".h", ".hh", ".hpp"}
SESSION_SCHEMA = 1
REPORT_SCHEMA = 1

