cmake_minimum_required(VERSION 3.1...3.20)
project(main)

include_directories(/usr/include/x86_64-linux-gnu/qt5/QtWidgets/)
#CGAL_Qt5 is needed for the drawing.
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# create_single_source_cgal_program("joe_simpson_visibility.cpp")
create_single_source_cgal_program("main.cpp")


# if(CGAL_Qt5_FOUND)
# #   #required to use basic_viewer
#   add_definitions(-DCGAL_USE_BASIC_VIEWER)
# endif()

# # #create the executable of the application

# if(CGAL_Qt5_FOUND)
#   #link it with the required CGAL libraries
#     # message(STATUS "Found CGAL_Qt5 library.")
#   target_link_libraries(joe_simpson_visibility PUBLIC CGAL::CGAL_Qt5)
# #   target_link_libraries(draw_polygon PUBLIC CGAL::CGAL_Basic_viewer)
# endif()