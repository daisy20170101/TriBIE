# Comparing MATLAB and Fortran Contributions for Points 8 & 9

## What Was Done

I've modified the code to output detailed breakdowns of the three contributions to the half-space solution:
1. **Main Dislocation** - Full-space triangular dislocation
2. **Harmonic Function** - Correction term
3. **Image Dislocation** - Mirror dislocation for half-space boundary condition

## Files Modified

### MATLAB Side
**File**: `TDstressHS.m`
- **Lines 105-108**: Output Main Dislocation contribution
- **Lines 112-115**: Output Harmonic Function contribution
- **Lines 123-126**: Output Image Dislocation contribution
- **Lines 138-142**: Output Total (sum of all three)

### Fortran Side
**File**: `sub_nikkhoo.f90`
- Already has output at lines 59, 71, 87 for the three contributions
- Used by `quick_test.f90` and `test_trimode_module.f90`

## How to Run the Comparison

### Step 1: Run Fortran Test
```bash
cd NikkhooWalter2015
./test_regularization.sh
```

This will output:
```
=== Main Dislocation Contribution ===
Strain: Exx= <value>

=== Harmonic Function Contribution ===
Strain: Exx= <value>

=== Image Dislocation Contribution ===
Strain: Exx= <value>
```

### Step 2: Run MATLAB Test
```matlab
cd NikkhooWalter2015
test_points_8_9_HS
```

This will output:
```
=== MATLAB Main Dislocation Contribution ===
Strain: Exx= <value>
        Eyy= <value>
        Ezz= <value>

=== MATLAB Harmonic Function Contribution ===
Strain: Exx= <value>
        Eyy= <value>
        Ezz= <value>

=== MATLAB Image Dislocation Contribution ===
Strain: Exx= <value>
        Eyy= <value>
        Ezz= <value>

=== MATLAB Total (Main + Harmonic + Image) ===
Strain: Exx= <value>
        Eyy= <value>
        Ezz= <value>
```

## What to Compare

### For Point 8: (3.0, -3.0, -6.0)

Compare each contribution between Fortran and MATLAB:
- **Main Dislocation Exx**: Should match (or identify difference)
- **Harmonic Function Exx**: Should match
- **Image Dislocation Exx**: Should match
- **Total Exx**: Expected = 7.064e-4

### For Point 9: (-3.0, 3.0, -3.0)

Compare each contribution:
- **Main Dislocation Exx**: Should match
- **Harmonic Function Exx**: Should match
- **Image Dislocation Exx**: Should match
- **Total Exx**: Expected = 2.113e-4

## Expected Outcomes

### If MATLAB Returns Finite Values but Fortran Returns NaN

**Likely causes:**
1. **Different coordinate transformations** - Check if MATLAB and Fortran compute the same TDCS coordinates
2. **Different angular dislocation handling** - One may avoid exact singularities through numerical precision
3. **IEEE arithmetic differences** - MATLAB may handle edge cases differently

**Where to look:**
- Compare the TDCS coordinates (x_td, y_td, z_td) between MATLAB and Fortran
- Compare the angular dislocation parameters (W, r-z, etc.) for each of the three angular dislocations
- Check which contribution(s) differ: Main, Harmonic, or Image

### If One Contribution is NaN

If only **Main Dislocation** or **Image Dislocation** returns NaN:
- The issue is in `TDstressFS` (full-space calculation)
- Check the angular dislocation calculations within that contribution
- The three angular dislocations are summed; one may be singular

If **Harmonic Function** returns NaN:
- The issue is in `TDstress_HarFunc`
- Different calculation method, may not have the same singularities

## Next Steps Based on Results

### Scenario 1: All MATLAB contributions are finite
→ The issue is specific to Fortran's numerical handling
→ Need to understand how MATLAB avoids the singularities

### Scenario 2: MATLAB also has NaN in same contribution(s)
→ The singularity is inherent to the geometry
→ Need a different mathematical approach

### Scenario 3: Different contributions are NaN
→ The implementations may differ in subtle ways
→ Need detailed step-by-step comparison of the calculations

## Quick Reference

**Current Fortran behavior** (from previous tests):
```
Point 8:
  Main Dislocation: NaN
  Harmonic Function: -8.262495e-04
  Image Dislocation: NaN
  Total: NaN

Point 9:
  Main Dislocation: NaN
  Harmonic Function: <value>
  Image Dislocation: NaN
  Total: NaN
```

**Expected MATLAB behavior**:
```
Point 8:
  Total: 7.064e-4 (confirmed from previous test)

Point 9:
  Total: 2.113e-4 (confirmed from previous test)
```

The breakdown will reveal:
- Are MATLAB's Main and Image contributions also problematic?
- Is MATLAB somehow getting cancellation between contributions?
- Or does MATLAB avoid the singularities entirely through different numerics?
