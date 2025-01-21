target_include_directories( simplexacore PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( simplexacore PRIVATE
   eb.cpp
   eb_box.cpp
   eb_terrain.cpp
   eb_cut_cell.cpp
   eb_aux.cpp
   )
