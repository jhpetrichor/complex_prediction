# CMake generated Testfile for 
# Source directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video
# Build directory: /home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/modules/video
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(opencv_test_video "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/bin/opencv_test_video" "--gtest_output=xml:opencv_test_video.xml")
set_tests_properties(opencv_test_video PROPERTIES  LABELS "Main;opencv_video;Accuracy" WORKING_DIRECTORY "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/test-reports/accuracy" _BACKTRACE_TRIPLES "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVUtils.cmake;1795;add_test;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1375;ocv_add_test_from_target;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1133;ocv_add_accuracy_tests;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video/CMakeLists.txt;2;ocv_define_module;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video/CMakeLists.txt;0;")
add_test(opencv_perf_video "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/bin/opencv_perf_video" "--gtest_output=xml:opencv_perf_video.xml")
set_tests_properties(opencv_perf_video PROPERTIES  LABELS "Main;opencv_video;Performance" WORKING_DIRECTORY "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/test-reports/performance" _BACKTRACE_TRIPLES "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVUtils.cmake;1795;add_test;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1274;ocv_add_test_from_target;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1134;ocv_add_perf_tests;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video/CMakeLists.txt;2;ocv_define_module;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video/CMakeLists.txt;0;")
add_test(opencv_sanity_video "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/bin/opencv_perf_video" "--gtest_output=xml:opencv_perf_video.xml" "--perf_min_samples=1" "--perf_force_samples=1" "--perf_verify_sanity")
set_tests_properties(opencv_sanity_video PROPERTIES  LABELS "Main;opencv_video;Sanity" WORKING_DIRECTORY "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/build/test-reports/sanity" _BACKTRACE_TRIPLES "/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVUtils.cmake;1795;add_test;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1275;ocv_add_test_from_target;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/cmake/OpenCVModule.cmake;1134;ocv_add_perf_tests;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video/CMakeLists.txt;2;ocv_define_module;/home/jh/code/complex_predict/dataset/Yeast/complex/opencv-4.x/modules/video/CMakeLists.txt;0;")
