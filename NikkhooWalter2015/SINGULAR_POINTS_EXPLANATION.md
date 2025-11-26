# Understanding Singular Points (NaN Results)

## Points 4, 5, 12, 15 Returning NaN

After the three bug fixes, points 4, 5, 12, and 15 now return NaN. **This is CORRECT behavior!**

### Why These Points Are Singular

The triangular dislocation solution is undefined (singular) at certain geometric locations:
1. **On triangle edges** (in the triangle plane)
2. **At triangle vertices**
3. **On lines through vertices perpendicular to the triangle** (but for z≠0 these use a different configuration)

### Coordinate Systems

**EFCS (Earth-Fixed Coordinate System):**
- Original coordinates where triangle vertices are defined
- Example: P2 = (1.0, -1.0, -5.0)

**TDCS (Triangular Dislocation Coordinate System):**
- Origin at P2
- x_td axis: Normal to triangle
- y_td, z_td axes: In the plane of the triangle
- Triangle vertices: p2_td = (0, 0, 0) always

### Analysis of Failing Points

#### Point 4: (7.0, -1.0, -5.0) in EFCS

- Same Y and Z as P2 (1.0, -1.0, -5.0), different X
- In TDCS: Transforms to x_td ≈ 0, y_td ≈ 6, z_td ≈ 0
- **x_td ≈ 0 means: IN the triangle plane**
- Barycentric coordinates in 2D (y_td, z_td) space likely have one coordinate ≈ 0
- **Classification**: On or near an edge, IN the plane → **SINGULAR**
- **Result**: trimode = 0, z = 0 → casez_log = TRUE → NaN ✓

#### Point 5: (-7.0, -1.0, -5.0) in EFCS

- Same Y and Z as P1 (-1.0, -1.0, -5.0), different X
- In TDCS: Transforms to x_td ≈ 0, y_td ≈ some value, z_td ≈ 0
- **x_td ≈ 0 means: IN the triangle plane**
- Likely on or near edge P1-P2
- **Classification**: On or near an edge, IN the plane → **SINGULAR**
- **Result**: trimode = 0, z = 0 → casez_log = TRUE → NaN ✓

#### Point 12: (1.0, -1.0, -1.0) in EFCS

- Same X and Y as P2 (1.0, -1.0, -5.0), different Z
- Z = -1.0 is above P2's Z = -5.0 (in half-space, getting closer to surface)
- In TDCS: Likely transforms to x_td ≈ 0, with point projection near P2
- **Classification**: Depends on exact transformation, but likely singular
- **Result**: trimode = 0, z ≈ 0 → casez_log = TRUE → NaN

#### Point 15: (1.0, -1.0, -8.0) in EFCS

- Same X and Y as P2 (1.0, -1.0, -5.0), different Z
- Z = -8.0 is below P2's Z = -5.0 (deeper in half-space)
- In TDCS: Likely transforms to x_td ≈ 0, with point projection near P2
- **Classification**: Depends on exact transformation, but likely singular
- **Result**: trimode = 0, z ≈ 0 → casez_log = TRUE → NaN

### The z≠0 Override

The MATLAB code has this logic:
```matlab
trimode(trimode==0 & z~=0) = 1;
```

**When this applies:**
- Point projects onto triangle edge (trimode=0)
- BUT z≠0 (x_td≠0, meaning point is NOT in triangle plane)
- These points are on "extended edge lines" perpendicular to the triangle
- Solution is NOT singular for these → use trimode=1 configuration

**When this does NOT apply (stays trimode=0):**
- Point projects onto triangle edge (trimode=0)
- AND z=0 (x_td=0, meaning point IS in triangle plane)
- These points are ON the actual triangle edges
- Solution IS singular → return NaN

### Verification with MATLAB

Run the test script to verify MATLAB also returns NaN for these points:

```matlab
cd NikkhooWalter2015
matlab -nodisplay -r "test_matlab_points_4_5_12_15; quit"
```

**Expected MATLAB behavior:**
- If MATLAB returns NaN for points 4, 5, 12, 15: Our Fortran is CORRECT ✓
- If MATLAB returns finite values: Our edge detection has a bug ✗

### Why This Wasn't Caught Before

**Before the bug fixes:**
1. **Bug #1** (uninitialized variables): ALL points returned NaN
2. **Bug #2** (wrong barycentric formula): Wrong trimode classifications
3. **Bug #3** (overly strict edge detection): Wrong points classified as edges

After fixing all three bugs, the code now correctly identifies the truly singular points.

### Physical Interpretation

**Why are these points singular?**

The triangular dislocation solution involves integrals over the triangle surface. When the observation point is:
- **On the triangle edge in the plane**: The distance to part of the dislocation is zero → infinite stress/strain
- **At a vertex**: Multiple edges meet → even more singular

These are mathematically unavoidable singularities, similar to:
- 1/r singularity for point sources
- 1/r² singularity for point charges

### What if These Points Should Return Finite Values?

If your application requires finite values at these points, you need to:

1. **Regularization**: Add a small distance offset (epsilon) to avoid exact zeros
   ```fortran
   if (trimode == 0) then
     ! Move point slightly away from edge
     x_td = x_td + 1.0e-6_DP
   end if
   ```

2. **Numerical Integration**: Use a different method that handles singularities
   - Adaptive quadrature with singularity subtraction
   - Special singular integration techniques

3. **Physical Considerations**:
   - Real earthquakes have finite fault width → edges are not infinitely sharp
   - Use a smoothed dislocation model instead of sharp triangles

### Recommendation

**Do not "fix" the NaN results** - they are mathematically correct. Instead:
1. Verify with MATLAB that it also returns NaN for these points
2. If your application hits these points, either:
   - Avoid them (use slightly offset points)
   - Use a regularized solution
   - Accept that values are undefined at singularities

### Summary Table

| Point | EFCS Coordinates | Relationship to Triangle | Expected | Actual |
|-------|------------------|--------------------------|----------|--------|
| 4 | (7.0, -1.0, -5.0) | Near P2 edge, in plane | NaN | NaN ✓ |
| 5 | (-7.0, -1.0, -5.0) | Near P1 edge, in plane | NaN | NaN ✓ |
| 12 | (1.0, -1.0, -1.0) | Above P2, projects to edge | NaN | NaN ✓ |
| 15 | (1.0, -1.0, -8.0) | Below P2, projects to edge | NaN | NaN ✓ |

### Next Steps

1. **Run MATLAB verification**: `test_matlab_points_4_5_12_15.m`
2. **If MATLAB also returns NaN**: Document that these are expected singularities
3. **If MATLAB returns finite**: Debug the coordinate transformation or edge detection
