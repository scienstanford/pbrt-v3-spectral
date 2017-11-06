# CMake generated Testfile for 
# Source directory: /Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/ptex/src/tests
# Build directory: /Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(wtest "/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests/wtest")
add_test(rtest "/usr/local/Cellar/cmake/3.9.4_1/bin/cmake" "-DOUT=/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests/rtest.out" "-DDATA=/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/ptex/src/tests/rtestok.dat" "-DCMD=/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests/rtest" "-P" "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/ptex/src/tests/compare_test.cmake")
add_test(ftest "/usr/local/Cellar/cmake/3.9.4_1/bin/cmake" "-DOUT=/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests/ftest.out" "-DDATA=/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/ptex/src/tests/ftestok.dat" "-DCMD=/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests/ftest" "-P" "/Users/trishalian/GitRepos/pbrt-v3-spectral/src/ext/ptex/src/tests/compare_test.cmake")
add_test(halftest "/Users/trishalian/GitRepos/pbrt-v3-spectral/bin/src/ext/ptex/src/tests/halftest")
