# Cell_Bucket unit tests:
add_executable(unit_test_Cell_Bucket Cell_Bucket_test.cpp)
target_link_libraries(unit_test_Cell_Bucket PRIVATE src)
add_test(NAME unit_test_Cell_Bucket COMMAND unit_test_Cell_Bucket WORKING_DIRECTORY bin)

add_executable(unit_test_Cell_ECM Cell_ECM_test.cpp)
target_link_libraries(unit_test_Cell_ECM PRIVATE src)
add_test(NAME unit_test_Cell_ECM COMMAND unit_test_Cell_ECM WORKING_DIRECTORY bin)

add_executable(unit_test_Cell_SPM Cell_SPM_test.cpp)
target_link_libraries(unit_test_Cell_SPM PRIVATE src)
add_test(NAME unit_test_Cell_SPM COMMAND unit_test_Cell_SPM WORKING_DIRECTORY bin)

