# CMake generated Testfile for 
# Source directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/highgui
# Build directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/modules/highgui
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(opencv_test_highgui "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/bin/opencv_test_highgui" "--gtest_output=xml:opencv_test_highgui.xml")
set_tests_properties(opencv_test_highgui PROPERTIES  LABELS "Main;opencv_highgui;Accuracy" WORKING_DIRECTORY "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/test-reports/accuracy" _BACKTRACE_TRIPLES "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVUtils.cmake;1795;add_test;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1375;ocv_add_test_from_target;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/highgui/CMakeLists.txt;309;ocv_add_accuracy_tests;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/highgui/CMakeLists.txt;0;")
