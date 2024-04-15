# CMake generated Testfile for 
# Source directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/flann
# Build directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/modules/flann
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(opencv_test_flann "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/bin/opencv_test_flann" "--gtest_output=xml:opencv_test_flann.xml")
set_tests_properties(opencv_test_flann PROPERTIES  LABELS "Main;opencv_flann;Accuracy" WORKING_DIRECTORY "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/test-reports/accuracy" _BACKTRACE_TRIPLES "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVUtils.cmake;1795;add_test;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1375;ocv_add_test_from_target;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1133;ocv_add_accuracy_tests;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/flann/CMakeLists.txt;2;ocv_define_module;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/flann/CMakeLists.txt;0;")
