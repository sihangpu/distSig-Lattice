set(FIPS202_SRCS fips202.c)

# Tests
set(TEST_SRCS test/test_mulsig.cpp)


if(MSVC)
  add_compile_options(/nologo /O2 /W4 /wd4146 /wd4244)
else()
  add_compile_options(-Wall)
  add_compile_options(-Wmissing-prototypes -Wredundant-decls -Wpointer-arith)
  add_compile_options(-I/opt/homebrew/include -I/usr/local/include)
  add_compile_options(-O3 -fomit-frame-pointer)
  add_link_options(-L/opt/homebrew/lib -L/usr/local/lib)
endif()


add_library(fips202_ref ${FIPS202_SRCS})


# Tests
add_executable(test_mulsig_exec ${TEST_SRCS})
target_link_libraries(test_mulsig_exec ntl gmp mpfr fips202_ref)