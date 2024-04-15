# CMake generated Testfile for 
# Source directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/ml
# Build directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/modules/ml
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(opencv_test_ml "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/bin/opencv_test_ml" "--gtest_output=xml:opencv_test_ml.xml")
set_tests_properties(opencv_test_ml PROPERTIES  LABELS "Main;opencv_ml;Accuracy" WORKING_DIRECTORY "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/test-reports/accuracy" _BACKTRACE_TRIPLES "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVUtils.cmake;1795;add_test;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1375;ocv_add_test_from_target;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1133;ocv_add_accuracy_tests;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/ml/CMakeLists.txt;2;ocv_define_module;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/ml/CMakeLists.txt;0;")
