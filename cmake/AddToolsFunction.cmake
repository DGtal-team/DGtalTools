# DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD must be a global variable between 0 and 99

function(DGtalTools_add_tool tool_file) #optional_avoid_add_test
  string(RANDOM LENGTH 2 ALPHABET "0123456789" _random)
  if (${_random} LESS ${DGTALTOOLS_RANDOMIZED_BUILD_THRESHOLD}  OR ${tool_file} IN_LIST DGTALTOOLS_RANDOMIZED_BUILD_WHITELIST)
    add_executable(${tool_file} ${tool_file}.cpp)
    target_link_libraries (${tool_file} ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})
     install(TARGETS ${tool_file}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
  endif()
endfunction()