add_library(settings INTERFACE)

target_include_directories(settings INTERFACE .)  # #CHECK try to make this PUBLIC so other libraries use headers without relative path.
