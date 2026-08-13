#!/usr/bin/env python3
"""
Script to fix XDMF formatting in 3dtri_BP5.f90 by removing extra leading spaces
"""

import re

def fix_xdmf_formatting():
    # Read the file
    with open('3dtri_BP5.f90', 'r') as f:
        content = f.read()
    
    # Fix the write statements by removing extra leading spaces
    # Pattern: find write statements with extra spaces and fix them
    patterns = [
        # Fix the main XDMF structure
        (r'       write\(99,\*\) \'   <Grid Name="step_\', i, \'" GridType="Uniform">\'', 
         '          write(99,*) \'  <Grid Name="step_\', i, \'" GridType="Uniform">\''),
        
        # Fix the closing tags
        (r'       write\(99,\*\) \' </Grid>\'', 
         '       write(99,*) \' </Grid>\''),
        (r'       write\(99,\*\) \'</Domain>\'', 
         '       write(99,*) \'</Domain>\''),
        (r'       write\(99,\*\) \'</Xdmf>\'', 
         '       write(99,*) \'</Xdmf>\''),
    ]
    
    # Apply fixes
    for pattern, replacement in patterns:
        content = re.sub(pattern, replacement, content)
    
    # Write back to file
    with open('3dtri_BP5.f90', 'w') as f:
        f.write(content)
    
    print("XDMF formatting fixed!")

if __name__ == "__main__":
    fix_xdmf_formatting()
