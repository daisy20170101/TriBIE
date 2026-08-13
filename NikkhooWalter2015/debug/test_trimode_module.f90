program test_trimode_module
  use nikkhoo_walter
  use, intrinsic :: ieee_arithmetic
  implicit none

  ! Triangle vertices
  real(DP), dimension(3) :: p1, p2, p3

  ! Test points 8 and 9
  real(DP) :: x, y, z
  real(DP) :: ss, ds, ts, mu, lambda
  real(DP), dimension(6) :: stress, strain

  ! Set up triangle vertices
  p1 = [-1.0_DP, -1.0_DP, -5.0_DP]
  p2 = [1.0_DP, -1.0_DP, -5.0_DP]
  p3 = [-1.0_DP, 1.0_DP, -4.0_DP]

  ! Set slip components
  ss = 1.0_DP
  ds = -1.0_DP
  ts = 2.0_DP

  ! Set elastic parameters
  mu = 3.0e10_DP
  lambda = 3.0e10_DP

  print *, '============================================================'
  print *, 'Testing module version (using nikkhoo_walter module)'
  print *, '============================================================'
  print *, ''
  print *, 'This test uses the ACTUAL compiled module, not standalone code.'
  print *, 'If this shows NaN, the module was not properly recompiled.'
  print *, ''

  ! Test point 8
  x = 3.0_DP
  y = -3.0_DP
  z = -6.0_DP

  print *, '========== POINT 8 =========='
  print *, 'Coords: x=', x, ', y=', y, ', z=', z

  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)

  print *, 'Strain(1) Exx = ', strain(1)
  if (ieee_is_nan(strain(1))) then
    print *, '  ✗ FAIL: Got NaN (module still has old code without bounds checking)'
  else
    print *, '  ✓ PASS: Got valid number (module has the fix)'
  end if
  print *, ''

  ! Test point 9
  x = -3.0_DP
  y = 3.0_DP
  z = -3.0_DP

  print *, '========== POINT 9 =========='
  print *, 'Coords: x=', x, ', y=', y, ', z=', z

  call tdstress_hs(x, y, z, p1, p2, p3, ss, ds, ts, mu, lambda, stress, strain)

  print *, 'Strain(1) Exx = ', strain(1)
  if (ieee_is_nan(strain(1))) then
    print *, '  ✗ FAIL: Got NaN (module still has old code without bounds checking)'
  else
    print *, '  ✓ PASS: Got valid number (module has the fix)'
  end if
  print *, ''

  print *, '============================================================'
  print *, ''
  print *, 'If you see NaN above, the problem is:'
  print *, '  1. Old nikkhoo_walter.mod or sub_nikkhoo.o files exist'
  print *, '  2. They are being used instead of newly compiled ones'
  print *, '  3. Or compilation failed silently'
  print *, ''
  print *, 'Solution:'
  print *, '  cd NikkhooWalter2015'
  print *, '  rm -f *.o *.mod test_* debug_*'
  print *, '  gfortran -c sub_nikkhoo.f90'
  print *, '  gfortran -o test_trimode_module test_trimode_module.f90 sub_nikkhoo.o'
  print *, '  ./test_trimode_module'
  print *, ''
  print *, '============================================================'

end program test_trimode_module
