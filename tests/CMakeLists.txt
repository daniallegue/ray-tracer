add_executable(FinalProjectTests 
  "your_tests.cpp"
)

target_compile_features(FinalProjectTests PRIVATE cxx_std_20)
set_project_warnings(FinalProjectTests)
target_link_libraries(FinalProjectTests PRIVATE CGFramework FinalProjectLib Catch2WithMain)
