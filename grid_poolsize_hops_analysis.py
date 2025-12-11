import configparser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from main import main
import os
import time
import statistics

def run_simulation_with_threshold(threshold_value, n_hops_value, runs=10):
    """Run simulation `runs` times with a given threshold and n_hops, return averaged stats."""
    print(f"[ANALYSIS] Starting simulations for threshold={threshold_value}, n_hops={n_hops_value} (runs={runs})")
    config = configparser.ConfigParser()
    config.read('ConfigFile.ini')

    # Ensure MIXING section exists
    if 'MIXING' not in config:
        config.add_section('MIXING')
    config.set('MIXING', 'threshold', str(threshold_value))
    
    # Update n_hops
    if 'TOPOLOGY' not in config:
        config.add_section('TOPOLOGY')
    config.set('TOPOLOGY', 'n_hops', str(n_hops_value))
    
    # Ensure no external clients when mixes act as clients
    if 'DEFAULT' not in config:
        config['DEFAULT'] = {}
    config['DEFAULT']['n_clients'] = '0'

    temp_config_name = f'temp_config_thresh_{threshold_value}_hops_{n_hops_value}.ini'
    with open(temp_config_name, 'w') as configfile:
        config.write(configfile)

    mean_list = []
    median_list = []
    q25_list = []
    sim_time_list = []

    try:
        # Backup and swap config once for all runs
        if os.path.exists('ConfigFile_backup.ini'):
            os.remove('ConfigFile_backup.ini')
        os.rename('ConfigFile.ini', 'ConfigFile_backup.ini')
        os.rename(temp_config_name, 'ConfigFile.ini')

        for run_idx in range(1, runs + 1):
            print(f"[RUN] threshold={threshold_value}, n_hops={n_hops_value}, run {run_idx}/{runs}")
            start_run = time.time()
            result = main(3)  # keep calling existing main
            run_time = time.time() - start_run

            # robust parsing of result (support list/tuple/dict)
            mean_entropy = 0
            median_entropy = 0
            q25_entropy = 0
            
            if isinstance(result, dict):
                mean_entropy = result.get('mean_entropy', 0)
                median_entropy = result.get('median_entropy', 0)
                q25_entropy = result.get('q25_entropy', 0)
            elif isinstance(result, (list, tuple)) and len(result) >= 4:
                mean_entropy = result[1] or 0
                median_entropy = result[2] or 0
                q25_entropy = result[3] or 0
            else:
                # fallback: try to extract from simulation.Log if main returns simulation object
                try:
                    mean_entropy = float(getattr(result, 'mean_entropy', 0))
                    median_entropy = float(getattr(result, 'median_entropy', 0))
                    q25_entropy = float(getattr(result, 'q25_entropy', 0))
                except Exception:
                    pass

            mean_list.append(mean_entropy)
            median_list.append(median_entropy)
            q25_list.append(q25_entropy)
            sim_time_list.append(run_time)

            print(f"[RUN-RESULT] run {run_idx}: mean={mean_entropy:.4f}, median={median_entropy:.4f}, q25={q25_entropy:.4f}, time={run_time:.2f}s")

    finally:
        # restore original config
        if os.path.exists('ConfigFile.ini'):
            os.remove('ConfigFile.ini')
        if os.path.exists('ConfigFile_backup.ini'):
            os.rename('ConfigFile_backup.ini', 'ConfigFile.ini')
        # remove temp if still present
        if os.path.exists(temp_config_name):
            os.remove(temp_config_name)

    # compute averages
    avg_mean = statistics.mean(mean_list) if mean_list else 0
    avg_median = statistics.mean(median_list) if median_list else 0
    avg_q25 = statistics.mean(q25_list) if q25_list else 0
    avg_time = statistics.mean(sim_time_list) if sim_time_list else 0
    
    # compute standard deviations
    std_mean = statistics.stdev(mean_list) if len(mean_list) > 1 else 0
    std_median = statistics.stdev(median_list) if len(median_list) > 1 else 0
    std_q25 = statistics.stdev(q25_list) if len(q25_list) > 1 else 0

    return {
        'n_hops': n_hops_value,
        'threshold': threshold_value,
        'mean_entropy': avg_mean,
        'median_entropy': avg_median,
        'q25_entropy': avg_q25,
        'simulation_time': avg_time,
        'std_mean': std_mean,
        'std_median': std_median,
        'std_q25': std_q25,
        'raw_means': mean_list,
        'raw_medians': median_list,
        'raw_q25s': q25_list,
        'raw_times': sim_time_list
    }

def create_entropy_statistics_plot_multi_hops(runs_per_threshold=10):
    """Main analysis: run multiple n_hops values, each with multiple pool sizes."""
    
    # Parameters to sweep
    n_hops_values = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    pool_sizes = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    
    print(f"[ANALYSIS] Testing n_hops: {n_hops_values}")
    print(f"[ANALYSIS] Testing pool sizes: {pool_sizes}")
    print(f"[ANALYSIS] Total simulations: {len(n_hops_values)} × {len(pool_sizes)} × {runs_per_threshold} runs")

    all_results = []
    
    # Loop through each n_hops value
    for n_hops in n_hops_values:
        print(f"\n{'='*80}")
        print(f"[PROGRESS] Starting n_hops = {n_hops}")
        print(f"{'='*80}")
        
        hops_results = []
        
        # For each n_hops, loop through pool sizes
        for i, threshold in enumerate(pool_sizes, 1):
            print(f"\n[PROGRESS] n_hops={n_hops}: Running threshold {threshold} ({i}/{len(pool_sizes)})")
            start_time = time.time()
            
            try:
                res = run_simulation_with_threshold(threshold, n_hops, runs=runs_per_threshold)
                hops_results.append(res)
                all_results.append(res)
                
                elapsed = time.time() - start_time
                print(f"[PROGRESS] n_hops={n_hops}, threshold={threshold}: "
                      f"mean={res['mean_entropy']:.4f}±{res['std_mean']:.4f}, "
                      f"time={elapsed:.2f}s")
                
            except Exception as e:
                print(f"[ERROR] n_hops={n_hops}, threshold={threshold} failed: {e}")
                import traceback
                traceback.print_exc()
                
                hops_results.append({
                    'n_hops': n_hops,
                    'threshold': threshold,
                    'mean_entropy': 0,
                    'median_entropy': 0,
                    'q25_entropy': 0,
                    'simulation_time': 0,
                    'std_mean': 0,
                    'std_median': 0,
                    'std_q25': 0
                })
        
        print(f"\n[PROGRESS] Completed n_hops = {n_hops}")

    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    df.to_csv('entropy_vs_pool_size_multi_hops.csv', index=False)
    print(f"\n[ANALYSIS] Raw data saved to 'entropy_vs_pool_size_multi_hops.csv'")
    
    # Create plots
    create_multi_hops_plots(df, n_hops_values, pool_sizes)
    
    # Print summary table
    print_multi_hops_summary(df, n_hops_values)

def create_multi_hops_plots(df, n_hops_values, pool_sizes):
    """Create comprehensive plots for multi-hops analysis"""
    
    # Filter out failed simulations
    valid_df = df[df['mean_entropy'] > 0]
    
    if len(valid_df) == 0:
        print("[WARNING] No valid simulation results to plot!")
        return
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    
    # Plot 1: Separate lines for each n_hops value
    ax1 = plt.subplot(2, 2, 1)
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(n_hops_values)))
    
    for idx, n_hops in enumerate(n_hops_values):
        hops_data = valid_df[valid_df['n_hops'] == n_hops].sort_values('threshold')
        
        if len(hops_data) > 0:
            ax1.plot(hops_data['threshold'], hops_data['mean_entropy'], 
                    'o-', label=f'n_hops={n_hops}', linewidth=2, markersize=6,
                    color=colors[idx])
    
    ax1.set_xlabel('Pool Size (Threshold)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean Entropy (bits)', fontsize=12, fontweight='bold')
    ax1.set_title('Mean Entropy vs Pool Size\n(Different n_hops)', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9, loc='lower right')
    
    # Plot 2: Heatmap of entropy for n_hops vs pool_size
    ax2 = plt.subplot(2, 2, 2)
    
    # Create pivot table for heatmap
    pivot_data = valid_df.pivot_table(
        values='mean_entropy',
        index='n_hops',
        columns='threshold',
        aggfunc='mean'
    )
    
    im = ax2.imshow(pivot_data.values, cmap='RdYlGn', aspect='auto')
    ax2.set_xticks(np.arange(len(pivot_data.columns)))
    ax2.set_yticks(np.arange(len(pivot_data.index)))
    ax2.set_xticklabels(pivot_data.columns, fontsize=8, rotation=45)
    ax2.set_yticklabels(pivot_data.index, fontsize=9)
    ax2.set_xlabel('Pool Size (Threshold)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Hops (n_hops)', fontsize=12, fontweight='bold')
    ax2.set_title('Entropy Heatmap: n_hops vs Pool Size', fontsize=13, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax2)
    cbar.set_label('Mean Entropy (bits)', fontsize=10)
    
    # Plot 3: Max entropy for each n_hops
    ax3 = plt.subplot(2, 2, 3)
    
    max_entropy_per_hops = []
    optimal_threshold_per_hops = []
    
    for n_hops in n_hops_values:
        hops_data = valid_df[valid_df['n_hops'] == n_hops]
        if len(hops_data) > 0:
            max_idx = hops_data['mean_entropy'].idxmax()
            max_entropy_per_hops.append(hops_data.loc[max_idx, 'mean_entropy'])
            optimal_threshold_per_hops.append(hops_data.loc[max_idx, 'threshold'])
        else:
            max_entropy_per_hops.append(0)
            optimal_threshold_per_hops.append(0)
    
    ax3.plot(n_hops_values, max_entropy_per_hops, 'bs-', linewidth=2.5, markersize=8)
    ax3.set_xlabel('Number of Hops (n_hops)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Maximum Mean Entropy (bits)', fontsize=12, fontweight='bold')
    ax3.set_title('Maximum Entropy vs Number of Hops', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_xticks(n_hops_values)
    
    # Annotate optimal points
    for i, (hops, entropy, threshold) in enumerate(zip(n_hops_values, max_entropy_per_hops, optimal_threshold_per_hops)):
        if entropy > 0:
            ax3.annotate(f'@{int(threshold)}', xy=(hops, entropy),
                        xytext=(0, 5), textcoords='offset points',
                        fontsize=8, ha='center')
    
    # Plot 4: Optimal threshold for each n_hops
    ax4 = plt.subplot(2, 2, 4)
    
    ax4.plot(n_hops_values, optimal_threshold_per_hops, 'rs-', linewidth=2.5, markersize=8)
    ax4.set_xlabel('Number of Hops (n_hops)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Optimal Pool Size (Threshold)', fontsize=12, fontweight='bold')
    ax4.set_title('Optimal Pool Size vs Number of Hops', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_xticks(n_hops_values)
    
    plt.suptitle('Grid Topology: Entropy Analysis Across Multiple Hops\n(10 runs per configuration, Mix-as-Client Mode)', 
                 fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig('entropy_analysis_multi_hops.png', dpi=300, bbox_inches='tight')
    plt.savefig('entropy_analysis_multi_hops.pdf', bbox_inches='tight')
    plt.show()
    
    print(f"[ANALYSIS] Plots saved as 'entropy_analysis_multi_hops.png' and '.pdf'")

def print_multi_hops_summary(df, n_hops_values):
    """Print comprehensive summary table"""
    
    valid_df = df[df['mean_entropy'] > 0]
    
    if len(valid_df) == 0:
        print("[ERROR] No valid results to summarize")
        return
    
    print("\n" + "="*100)
    print("ENTROPY ANALYSIS SUMMARY - ACROSS MULTIPLE HOPS")
    print("="*100)
    
    summary_data = []
    
    for n_hops in n_hops_values:
        hops_data = valid_df[valid_df['n_hops'] == n_hops]
        
        if len(hops_data) == 0:
            continue
        
        max_idx = hops_data['mean_entropy'].idxmax()
        max_entropy = hops_data.loc[max_idx, 'mean_entropy']
        optimal_threshold = hops_data.loc[max_idx, 'threshold']
        std_at_optimal = hops_data.loc[max_idx, 'std_mean']
        
        min_entropy = hops_data['mean_entropy'].min()
        avg_entropy = hops_data['mean_entropy'].mean()
        
        summary_data.append({
            'n_hops': n_hops,
            'max_entropy': max_entropy,
            'optimal_threshold': optimal_threshold,
            'std_at_optimal': std_at_optimal,
            'avg_entropy': avg_entropy,
            'min_entropy': min_entropy
        })
        
        print(f"\nn_hops = {n_hops}:")
        print(f"  Maximum Entropy:      {max_entropy:.4f} bits (at pool_size={int(optimal_threshold)})")
        print(f"  Std Dev at optimal:   {std_at_optimal:.4f}")
        print(f"  Average Entropy:      {avg_entropy:.4f} bits")
        print(f"  Minimum Entropy:      {min_entropy:.4f} bits")
        print(f"  Entropy Range:        {min_entropy:.4f} - {max_entropy:.4f} (Δ={max_entropy-min_entropy:.4f})")
    
    print("\n" + "-"*100)
    print("COMPARISON ACROSS HOPS:")
    print("-"*100)
    print(f"{'n_hops':<8} {'Max Entropy':<15} {'Optimal Pool':<15} {'Avg Entropy':<15} {'Range':<15}")
    print("-"*100)
    
    for row in summary_data:
        print(f"{row['n_hops']:<8} {row['max_entropy']:<15.4f} {int(row['optimal_threshold']):<15} "
              f"{row['avg_entropy']:<15.4f} {row['max_entropy']-row['min_entropy']:<15.4f}")
    
    print("="*100)
    
    # Find best overall configuration
    best_idx = max(range(len(summary_data)), key=lambda i: summary_data[i]['max_entropy'])
    best_config = summary_data[best_idx]
    
    print(f"\n[BEST CONFIGURATION] n_hops={best_config['n_hops']}, pool_size={int(best_config['optimal_threshold'])}")
    print(f"                     Maximum Entropy: {best_config['max_entropy']:.4f} bits")

if __name__ == "__main__":
    print("[ENTROPY ANALYSIS] Multi-Hops Entropy vs Pool Size Analysis")
    print("[INFO] This will create plots showing entropy trends across different n_hops values")
    print("[INFO] Each configuration runs 10 times, results are averaged")
    print("[INFO] Expected runtime: 30-60 minutes depending on your system\n")
    
    start_time = time.time()
    
    try:
        create_entropy_statistics_plot_multi_hops(runs_per_threshold=10)
        
        total_time = time.time() - start_time
        print(f"\n[ANALYSIS] Complete! Total runtime: {total_time/60:.2f} minutes")
        print(f"[ANALYSIS] Results saved to:")
        print(f"  - entropy_vs_pool_size_multi_hops.csv (raw data)")
        print(f"  - entropy_analysis_multi_hops.png (plots)")
        print(f"  - entropy_analysis_multi_hops.pdf (plots)")
        
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Analysis stopped by user")
    except Exception as e:
        print(f"\n[ERROR] Analysis failed: {e}")
        import traceback
        traceback.print_exc()