#!/usr/bin/env python3
"""
Memory Cost Calculator for 3dtri_BP5
Calculates memory requirements for different problem sizes and MPI configurations
"""

def calculate_memory_cost(Nt_all, nmv, nas, ncos, nnul, nsse, n_obv, np1, np2, nprocs):
    """
    Calculate memory cost for 3dtri_BP5
    
    Args:
        Nt_all: Total number of fault elements
        nmv: Number of monitoring points
        nas: Number of average slip intervals
        ncos: Number of cosine slip intervals
        nnul: Number of null slip intervals
        nsse: Number of slow slip events
        n_obv: Number of observation points
        np1: Number of strike profiles
        np2: Number of dip profiles
        nprocs: Number of MPI processes
    
    Returns:
        Dictionary with memory breakdown
    """
    
    # Element sizes in bytes
    DP_SIZE = 8      # Double precision
    INT_SIZE = 4     # Integer
    LOG_SIZE = 1     # Logical
    
    # Calculate local_cells per process (assuming even distribution)
    local_cells = Nt_all // nprocs
    
    # Global arrays (all processes have copies)
    global_arrays = {
        'x_all': Nt_all * DP_SIZE,
        'xi_all': Nt_all * DP_SIZE,
        'z_all': Nt_all * DP_SIZE,
        'yt_all': 2 * Nt_all * DP_SIZE,
        'dydt_all': 2 * Nt_all * DP_SIZE,
        'yt_scale_all': 2 * Nt_all * DP_SIZE,
        'yt0_all': 2 * Nt_all * DP_SIZE,
        'phy1_all': Nt_all * DP_SIZE,
        'phy2_all': Nt_all * DP_SIZE,
        'tau1_all': Nt_all * DP_SIZE,
        'tau2_all': Nt_all * DP_SIZE,
        'slip_all': Nt_all * DP_SIZE,
        'slipinc_all': Nt_all * DP_SIZE,
        'slipds_all': Nt_all * DP_SIZE,
        'slipdsinc_all': Nt_all * DP_SIZE
    }
    
    # Monitoring arrays
    monitoring_arrays = {
        'outs1': nmv * 7 * 10 * DP_SIZE,
        'maxv': nmv * DP_SIZE,
        'maxnum': nmv * INT_SIZE,
        'moment': nmv * DP_SIZE,
        'tmv': nmv * DP_SIZE,
        'Trup': Nt_all * DP_SIZE,
        'rup': Nt_all * LOG_SIZE,
        'area': Nt_all * DP_SIZE
    }
    
    # Slip arrays (largest memory usage)
    slip_arrays = {
        'slipz1_inter': Nt_all * nas * DP_SIZE,
        'slipz1_cos': Nt_all * ncos * DP_SIZE,
        'slipave_inter': Nt_all * nas * DP_SIZE,
        'slipave_cos': Nt_all * ncos * DP_SIZE,
        'v_cos': Nt_all * ncos * DP_SIZE,
        'slip_cos': Nt_all * ncos * DP_SIZE,
        'v_nul': Nt_all * nnul * DP_SIZE,
        'slip_nul': Nt_all * nnul * DP_SIZE,
        'slipz1_tau': Nt_all * nsse * DP_SIZE,
        'slipz1_sse': Nt_all * nsse * DP_SIZE
    }
    
    # Observation arrays
    observation_arrays = {
        'surf1': n_obv * Nt_all * DP_SIZE,
        'surf2': n_obv * Nt_all * DP_SIZE,
        'surf3': n_obv * Nt_all * DP_SIZE,
        'obvs': nmv * 6 * n_obv * DP_SIZE,
        'obvstrk': nmv * 2 * np1 * DP_SIZE,
        'obvdp': nmv * 2 * np2 * DP_SIZE,
        'pstrk': np1 * INT_SIZE,
        'pdp': np2 * INT_SIZE
    }
    
    # Local arrays (per process)
    local_arrays = {
        'phy1': local_cells * DP_SIZE,
        'phy2': local_cells * DP_SIZE,
        'x': local_cells * DP_SIZE,
        'z': local_cells * DP_SIZE,
        'xi': local_cells * DP_SIZE,
        'cca': local_cells * DP_SIZE,
        'ccb': local_cells * DP_SIZE,
        'seff': local_cells * DP_SIZE,
        'xLf': local_cells * DP_SIZE,
        'tau1': local_cells * DP_SIZE,
        'tau2': local_cells * DP_SIZE,
        'tau0': local_cells * DP_SIZE,
        'slip': local_cells * DP_SIZE,
        'slipinc': local_cells * DP_SIZE,
        'slipds': local_cells * DP_SIZE,
        'slipdsinc': local_cells * DP_SIZE,
        'yt': 2 * local_cells * DP_SIZE,
        'dydt': 2 * local_cells * DP_SIZE,
        'yt_scale': 2 * local_cells * DP_SIZE,
        'yt0': 2 * local_cells * DP_SIZE,
        'sr': local_cells * DP_SIZE,
        'vi': local_cells * DP_SIZE
    }
    
    # CRITICAL: Stiffness matrices (largest per-process memory)
    stiffness_arrays = {
        'stiff': local_cells * Nt_all * DP_SIZE,
        'stiff2': local_cells * Nt_all * DP_SIZE
    }
    
    # MPI communication buffers
    mpi_arrays = {
        'sendcounts': nprocs * INT_SIZE,
        'displs': nprocs * INT_SIZE,
        'send_buffer': 2 * Nt_all * DP_SIZE,
        'recv_buffer': 2 * Nt_all * DP_SIZE
    }
    
    # Calculate totals
    total_global = sum(global_arrays.values())
    total_monitoring = sum(monitoring_arrays.values())
    total_slip = sum(slip_arrays.values())
    total_observation = sum(observation_arrays.values())
    total_local = sum(local_arrays.values())
    total_stiffness = sum(stiffness_arrays.values())
    total_mpi = sum(mpi_arrays.values())
    
    # Total memory per process
    total_per_process = (total_global + total_monitoring + total_slip + 
                        total_observation + total_local + total_stiffness + total_mpi)
    
    # Total memory across all processes
    total_all_processes = total_per_process * nprocs
    
    # Convert to MB and GB
    def bytes_to_mb(bytes_val):
        return bytes_val / (1024 * 1024)
    
    def bytes_to_gb(bytes_val):
        return bytes_val / (1024 * 1024 * 1024)
    
    return {
        'global_arrays_mb': bytes_to_mb(total_global),
        'monitoring_arrays_mb': bytes_to_mb(total_monitoring),
        'slip_arrays_mb': bytes_to_mb(total_slip),
        'observation_arrays_mb': bytes_to_mb(total_observation),
        'local_arrays_mb': bytes_to_mb(total_local),
        'stiffness_arrays_mb': bytes_to_mb(total_stiffness),
        'mpi_arrays_mb': bytes_to_mb(total_mpi),
        'total_per_process_mb': bytes_to_mb(total_per_process),
        'total_per_process_gb': bytes_to_gb(total_per_process),
        'total_all_processes_mb': bytes_to_mb(total_all_processes),
        'total_all_processes_gb': bytes_to_gb(total_all_processes),
        'local_cells_per_process': local_cells,
        'stiffness_memory_gb': bytes_to_gb(total_stiffness)
    }

def print_memory_report(params, results):
    """Print formatted memory report"""
    print("=" * 80)
    print("MEMORY COST ANALYSIS FOR 3DTRI_BP5")
    print("=" * 80)
    print(f"Problem Size: {params['Nt_all']:,} total fault elements")
    print(f"MPI Processes: {params['nprocs']}")
    print(f"Local Elements per Process: {results['local_cells_per_process']:,}")
    print()
    
    print("MEMORY BREAKDOWN PER PROCESS:")
    print("-" * 50)
    print(f"Global Arrays:           {results['global_arrays_mb']:8.2f} MB")
    print(f"Monitoring Arrays:       {results['monitoring_arrays_mb']:8.2f} MB")
    print(f"Slip Arrays:             {results['slip_arrays_mb']:8.2f} MB")
    print(f"Observation Arrays:      {results['observation_arrays_mb']:8.2f} MB")
    print(f"Local Arrays:            {results['local_arrays_mb']:8.2f} MB")
    print(f"Stiffness Matrices:      {results['stiffness_arrays_mb']:8.2f} MB")
    print(f"MPI Communication:       {results['mpi_arrays_mb']:8.2f} MB")
    print("-" * 50)
    print(f"TOTAL PER PROCESS:       {results['total_per_process_gb']:8.2f} GB")
    print()
    
    print("MEMORY REQUIREMENTS:")
    print("-" * 50)
    print(f"Per Process:             {results['total_per_process_gb']:8.2f} GB")
    print(f"All Processes:           {results['total_all_processes_gb']:8.2f} GB")
    print()
    
    print("PERFORMANCE NOTES:")
    print("-" * 50)
    print(f"• Stiffness matrices use {results['stiffness_memory_gb']:.2f} GB per process")
    print(f"• This is {results['stiffness_arrays_mb']/results['total_per_process_mb']*100:.1f}% of total memory")
    print(f"• Slip arrays use {results['slip_arrays_mb']:.2f} MB per process")
    print(f"• Local arrays scale with {results['local_cells_per_process']:,} elements per process")
    print()
    
    # Memory recommendations
    if results['total_per_process_gb'] > 1.0:
        print("⚠️  WARNING: High memory usage per process!")
        print("   Consider reducing problem size or increasing MPI processes")
    elif results['total_per_process_gb'] > 0.5:
        print("⚠️  NOTE: Moderate memory usage per process")
        print("   Monitor memory usage during execution")
    else:
        print("✅ Memory usage looks reasonable")
    
    print("=" * 80)

def main():
    """Main function with example calculations"""
    
    # Example 1: Your current parameters
    print("EXAMPLE 1: YOUR CURRENT PARAMETERS")
    params1 = {
        'Nt_all': 77376,
        'nmv': 500,
        'nas': 500,
        'ncos': 500,
        'nnul': 500,
        'nsse': 500,
        'n_obv': 500,
        'np1': 500,
        'np2': 500,
        'nprocs': 25
    }
    
    results1 = calculate_memory_cost(**params1)
    print_memory_report(params1, results1)
    
    print("\n" + "="*80 + "\n")
    
    # Example 2: Larger problem (2x elements)
    print("EXAMPLE 2: LARGER PROBLEM (2x ELEMENTS)")
    params2 = params1.copy()
    params2['Nt_all'] = 154752  # 2x larger
    params2['nprocs'] = 50      # 2x more processes
    
    results2 = calculate_memory_cost(**params2)
    print_memory_report(params2, results2)
    
    print("\n" + "="*80 + "\n")
    
    # Example 3: Smaller problem for testing
    print("EXAMPLE 3: SMALLER PROBLEM FOR TESTING")
    params3 = {
        'Nt_all': 10000,
        'nmv': 100,
        'nas': 100,
        'ncos': 100,
        'nnul': 100,
        'nsse': 100,
        'n_obv': 100,
        'np1': 100,
        'np2': 100,
        'nprocs': 4
    }
    
    results3 = calculate_memory_cost(**params3)
    print_memory_report(params3, results3)

if __name__ == "__main__":
    main()
