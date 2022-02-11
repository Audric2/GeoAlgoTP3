# CMake generated Testfile for 
# Source directory: /home/local.isima.fr/aucatinon/Bureau/TP3/src
# Build directory: /home/local.isima.fr/aucatinon/Bureau/TP3/build/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(compilation_of__genre "/usr/bin/cmake" "--build" "/home/local.isima.fr/aucatinon/Bureau/TP3/build" "--target" "genre")
set_tests_properties(compilation_of__genre PROPERTIES  FIXTURES_SETUP "genre" LABELS "src")
add_test(execution___of__genre "/home/local.isima.fr/aucatinon/Bureau/TP3/build/src/genre")
set_tests_properties(execution___of__genre PROPERTIES  DEPENDS "compilation_of__genre" FIXTURES_REQUIRED "src;genre" LABELS "src" WORKING_DIRECTORY "/home/local.isima.fr/aucatinon/Bureau/TP3/build/src/__exec_test_dir")
add_test(src_SetupFixture "/usr/bin/cmake" "-E" "copy_directory" "/home/local.isima.fr/aucatinon/Bureau/TP3/src" "/home/local.isima.fr/aucatinon/Bureau/TP3/build/src/__exec_test_dir")
set_tests_properties(src_SetupFixture PROPERTIES  FIXTURES_SETUP "src" LABELS "src")
add_test(src_CleanupFixture "/usr/bin/cmake" "-E" "remove_directory" "/home/local.isima.fr/aucatinon/Bureau/TP3/build/src/__exec_test_dir")
set_tests_properties(src_CleanupFixture PROPERTIES  FIXTURES_CLEANUP "src" LABELS "src")
