set(FIPS202_SRCS fips202.c)

# Tests
set(TEST_SRCS test/test_dsig.cpp)


if(APPLE)
  add_compile_options(-Wall)
  add_compile_options(-Wredundant-decls -Wpointer-arith)
  add_compile_options(-I/opt/homebrew/include -I/usr/local/include/)
  add_compile_options(-O3 -fomit-frame-pointer)
  #add_compile_options(-march=native)
  add_link_options(-L/opt/homebrew/lib -L/usr/local/lib)
elseif(UNIX)
  add_compile_options(-Wall)
  add_compile_options(-Wredundant-decls -Wpointer-arith)
  add_compile_options(-I/usr/local/include/)
  add_compile_options(-O3 -fomit-frame-pointer -pthread)
  add_compile_options(-march=native)

  add_link_options(-L/usr/local/lib -pthread)
else()
  message("not tested in windows environment!")
endif()


add_library(fips202_ref ${FIPS202_SRCS})


# Tests
add_executable(test_dsig_exec ${TEST_SRCS})
target_link_libraries(test_dsig_exec ntl gmp fips202_ref)
