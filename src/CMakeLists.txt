add_executable(bcnutau main.cc models.cc)
target_link_libraries(bcnutau PRIVATE BCNUTAU_project_options BCNUTAU_project_warnings
                              PUBLIC Chili::Chili lhapdf)
