get_filename_component(EXERCISENAME ${CMAKE_CURRENT_LIST_DIR} NAME)

add_executable(${EXERCISENAME}
        main.cpp
        mesh_processing.cpp
        mesh_processing.h
        viewer.cpp
        viewer.h)
target_link_libraries(${EXERCISENAME} surface_mesh)
# Lastly, additional libraries may have been built for you.  In addition to linking
# against NanoGUI, we need to link against those as well.
target_link_libraries(${EXERCISENAME} nanogui ${NANOGUI_EXTRA_LIBS})


#================================================================
# Unit Test
#================================================================
add_executable(remesh_test
        remesh_test.cpp
        mesh_processing.cpp
        mesh_processing.h
        viewer.cpp
        viewer.h)
target_link_libraries(remesh_test Catch2::Catch2 nanogui surface_mesh ${NANOGUI_EXTRA_LIBS})

include(CTest)
include(ParseAndAddCatchTests)
ParseAndAddCatchTests(remesh_test)
target_compile_definitions(remesh_test PUBLIC -DUNIT_TEST)