include(plugin)
openflipper_plugin(TYPES POLYMESH TRIANGLEMESH DEPS Gurobi)
set_source_files_properties(filename.c PROPERTIES LANGUAGE CXX)
add_definitions(-DTRILIBRARY -DANSI_DECLARATORS -DREAL=double -DVOID=int)
