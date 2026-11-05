# Developer tools (see docs/performance.md).
foreach(tool benchmark/operator_timing compare/dump_wave compare/compare_wave)
  get_filename_component(name ${tool} NAME)
  add_executable(spodvr_${name} ${tool}.cpp)
  target_link_libraries(spodvr_${name} PRIVATE spodvr_core spodvr_warnings spodvr_fp)
endforeach()

# GPU backend timing (with -DSPODVR_USE_CUDA=ON).
if(SPODVR_USE_CUDA)
  add_executable(spodvr_cuda_operator_timing benchmark/cuda_operator_timing.cpp)
  target_link_libraries(spodvr_cuda_operator_timing PRIVATE spodvr_cuda spodvr_warnings spodvr_fp)
endif()
