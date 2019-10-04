# CMake generated Testfile for 
# Source directory: /home/shash/Documents/anon/fullproof/proof/client
# Build directory: /home/shash/Documents/anon/fullproof/proof/build/client
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(client_test "cargo" "test")
set_tests_properties(client_test PROPERTIES  WORKING_DIRECTORY "/home/shash/Documents/anon/fullproof/proof/client")
