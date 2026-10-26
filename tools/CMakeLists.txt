# Developer tools (see docs/performance.md).
foreach(tool benchmark/operator_timing compare/dump_wave compare/compare_wave)
  get_filename_component(name ${tool} NAME)
  add_executable(spodvr_${name} ${tool}.cpp)
  target_link_libraries(spodvr_${name} PRIVATE spodvr_core spodvr_warnings spodvr_fp)
endforeach()
