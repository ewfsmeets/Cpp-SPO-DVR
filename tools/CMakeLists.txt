# Developer tools (docs/performance.md): spodvr_operator_timing, spodvr_dump_wave, spodvr_compare_wave,
# spodvr_crp_timing.
foreach(tool benchmark/operator_timing compare/dump_wave compare/compare_wave)
  get_filename_component(name ${tool} NAME)
  add_executable(spodvr_${name} ${tool}.cpp)
  target_link_libraries(spodvr_${name} PRIVATE spodvr_core spodvr_warnings spodvr_fp)
endforeach()

# CRP library timing (one thread): links only the library, with its strict FP flags.
add_executable(spodvr_crp_timing benchmark/crp_timing.cpp)
target_link_libraries(spodvr_crp_timing PRIVATE spodvr_crplib spodvr_warnings spodvr_crplib_fp)

# GPU backend timing (with -DSPODVR_USE_CUDA=ON).
if(SPODVR_USE_CUDA)
  add_executable(spodvr_cuda_operator_timing benchmark/cuda_operator_timing.cpp)
  target_link_libraries(spodvr_cuda_operator_timing PRIVATE spodvr_cuda spodvr_warnings spodvr_fp)
endif()
