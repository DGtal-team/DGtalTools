
add_executable(CumulativeSequenceTest CumulativeSequence.cpp CumulativeSequenceTest.cpp)
#target_link_libraries(CumulativeSequenceTest sequence)

add_executable(LUTBasedNSDistanceTransform LUTBasedNSDistanceTransform.cpp ImageFilter.cpp ImageWriter.cpp BaseDistanceDT.cpp D4DistanceDT.cpp D8DistanceDT.cpp RatioNSDistanceDT.cpp PeriodicNSDistanceDT.cpp CumulativeSequence.cpp)

target_link_libraries (LUTBasedNSDistanceTransform ${DGTAL_LIBRARIES} )

#if (WITH_PNG)
    find_package(png REQUIRED)
    if (PNG_FOUND)
	include_directories(${PNG_INCLUDE_DIRS})
	add_library(pngio PNGImageReader.cpp PNGImageWriter.cpp)
	target_link_libraries(LUTBasedNSDistanceTransform ${PNG_LIBRARIES})
	target_link_libraries(LUTBasedNSDistanceTransform pngio)
    endif (PNG_FOUND)
#endif (WITH_PNG)

install(TARGETS LUTBasedNSDistanceTransform 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)



