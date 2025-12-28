import configparser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from main import main
import os
import time
import statistics
import seaborn as sns
from itertools import product

def run_simulation_with_params(corrupt_mixes, threshold, n_hops, runs=20):
    """Run simulation multiple times with specific parameters and return averaged stats."""
    print(f"\n[ANALYSIS] Starting simulations for corrupt={corrupt_mixes}, threshold={threshold}, n_hops={n_hops} ({runs} runs)")
    
    config = configparser.ConfigParser()
    config.read('ConfigFile.ini')

    # Update all three parameters
    if 'THREATMODEL' not in config:
        config.add_section('THREATMODEL')
    config.set('THREATMODEL', 'corrupt_mixes', str(corrupt_mixes))
    
    if 'MIXING' not in config:
        config.add_section('MIXING')
    config.set('MIXING', 'threshold', str(threshold))
    
    if 'DEFAULT' not in config:
        config['DEFAULT'] = {}
    config['DEFAULT']['n_hops'] = str(n_hops)
    config['DEFAULT']['n_clients'] = '0'  # Mix-as-client mode

    temp_config_name = f'temp_config_c{corrupt_mixes}_t{threshold}_h{n_hops}.ini'
    with open(temp_config_name, 'w') as configfile:
        config.write(configfile)

    mean_list = []
    median_list = []
    q25_list = []
    sim_time_list = []

    try:
        # Backup and swap config
        if os.path.exists('ConfigFile_backup.ini'):
            os.remove('ConfigFile_backup.ini')
        os.rename('ConfigFile.ini', 'ConfigFile_backup.ini')
        os.rename(temp_config_name, 'ConfigFile.ini')

        for run_idx in range(1, runs + 1):
            print(f"[RUN] corrupt={corrupt_mixes}, threshold={threshold}, n_hops={n_hops}, run {run_idx}/{runs}")
            start_run = time.time()
            
            try:
                result = main(3)
                run_time = time.time() - start_run

                # Parse result
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

                mean_list.append(mean_entropy)
                median_list.append(median_entropy)
                q25_list.append(q25_entropy)
                sim_time_list.append(run_time)

                print(f"[RUN-RESULT] run {run_idx}: mean={mean_entropy:.4f}, median={median_entropy:.4f}, q25={q25_entropy:.4f}, time={run_time:.2f}s")
            
            except Exception as e:
                print(f"[ERROR] Run {run_idx} failed: {e}")
                # Add zeros for failed run
                mean_list.append(0)
                median_list.append(0)
                q25_list.append(0)
                sim_time_list.append(0)

    finally:
        # Restore original config
        if os.path.exists('ConfigFile.ini'):
            os.remove('ConfigFile.ini')
        if os.path.exists('ConfigFile_backup.ini'):
            os.rename('ConfigFile_backup.ini', 'ConfigFile.ini')
        if os.path.exists(temp_config_name):
            os.remove(temp_config_name)

    # Compute statistics
    avg_mean = statistics.mean(mean_list) if mean_list else 0
    avg_median = statistics.mean(median_list) if median_list else 0
    avg_q25 = statistics.mean(q25_list) if q25_list else 0
    avg_time = statistics.mean(sim_time_list) if sim_time_list else 0
    
    std_mean = statistics.stdev(mean_list) if len(mean_list) > 1 else 0
    std_median = statistics.stdev(median_list) if len(median_list) > 1 else 0
    std_q25 = statistics.stdev(q25_list) if len(q25_list) > 1 else 0

    return {
        'corrupt_mixes': corrupt_mixes,
        'threshold': threshold,
        'n_hops': n_hops,
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

def run_parameter_sweep(runs_per_config=20):
    """Main parameter sweep analysis"""
    
    # Define parameter ranges
    corrupt_mixes_values = list(range(10, 101, 5))  # [10, 15, 20, ..., 100]
    threshold_values = [70, 75, 80]
    n_hops_values = [4, 5, 6]
    
    # Calculate total combinations
    total_combinations = len(corrupt_mixes_values) * len(threshold_values) * len(n_hops_values)
    total_runs = total_combinations * runs_per_config
    
    print(f"[PARAMETER SWEEP] Starting comprehensive parameter analysis")
    print(f"[PARAMETERS]")
    print(f"  Corrupt mixes: {corrupt_mixes_values}")
    print(f"  Thresholds: {threshold_values}")
    print(f"  N-hops: {n_hops_values}")
    print(f"[SCALE]")
    print(f"  Total parameter combinations: {total_combinations}")
    print(f"  Runs per combination: {runs_per_config}")
    print(f"  Total simulation runs: {total_runs}")
    print(f"[ESTIMATED TIME] ~{total_runs * 0.5 / 60:.1f} minutes (assuming 30s per run)")
    
    all_results = []
    config_count = 0
    
    # Generate all combinations
    for corrupt, threshold, hops in product(corrupt_mixes_values, threshold_values, n_hops_values):
        config_count += 1
        print(f"\n{'='*80}")
        print(f"[PROGRESS] Configuration {config_count}/{total_combinations}")
        print(f"[PARAMS] corrupt_mixes={corrupt}, threshold={threshold}, n_hops={hops}")
        print(f"{'='*80}")
        
        start_time = time.time()
        
        try:
            result = run_simulation_with_params(corrupt, threshold, hops, runs=runs_per_config)
            all_results.append(result)
            
            elapsed = time.time() - start_time
            print(f"[COMPLETED] Configuration {config_count} in {elapsed:.2f}s")
            print(f"[RESULT] mean_entropy={result['mean_entropy']:.4f}±{result['std_mean']:.4f}")
            
        except Exception as e:
            print(f"[ERROR] Configuration {config_count} failed: {e}")
            import traceback
            traceback.print_exc()
            
            # Add failed result with zeros
            all_results.append({
                'corrupt_mixes': corrupt,
                'threshold': threshold,
                'n_hops': n_hops,
                'mean_entropy': 0,
                'median_entropy': 0,
                'q25_entropy': 0,
                'simulation_time': 0,
                'std_mean': 0,
                'std_median': 0,
                'std_q25': 0
            })
    
    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    
    # Save raw data
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    csv_filename = f'grid_corrupt_mix_results_{timestamp}.csv'
    df.to_csv(csv_filename, index=False)
    print(f"\n[ANALYSIS] Raw data saved to '{csv_filename}'")
    
    # Generate comprehensive plots
    create_analysis_plots(df, corrupt_mixes_values, threshold_values, n_hops_values, timestamp)
    
    # Print summary statistics
    print_parameter_sweep_summary(df)
    
    return df

def create_analysis_plots(df, corrupt_values, threshold_values, n_hops_values, timestamp):
    """Create comprehensive visualization of parameter sweep results"""
    
    # Filter valid results
    valid_df = df[df['mean_entropy'] > 0]
    
    if len(valid_df) == 0:
        print("[WARNING] No valid results to plot!")
        return
    
    # Set up the plot style
    sns.set_style("whitegrid")
    
    # Create a large figure with multiple subplots
    fig = plt.figure(figsize=(18, 12))
    
    # ============ PLOT 1: Entropy vs Corruption (separate lines for threshold) ============
    ax1 = plt.subplot(2, 3, 1)
    for threshold in threshold_values:
        for hops in n_hops_values:
            subset = valid_df[(valid_df['threshold'] == threshold) & (valid_df['n_hops'] == hops)]
            if len(subset) > 0:
                ax1.plot(subset['corrupt_mixes'], subset['mean_entropy'], 
                        'o-', label=f'T={threshold}, H={hops}', linewidth=2, markersize=5)
    
    ax1.set_xlabel('Corrupt Mixes', fontsize=10, fontweight='bold')
    ax1.set_ylabel('Mean Entropy (bits)', fontsize=10, fontweight='bold')
    ax1.set_title('Entropy vs Corruption Level', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=7, loc='best')
    ax1.grid(True, alpha=0.3)
 
    # ============ PLOT 2: Heatmap - Corruption vs Threshold (averaged over n_hops) ============
    ax4 = plt.subplot(2, 3, 2)
    pivot_data = valid_df.pivot_table(
        values='mean_entropy',
        index='corrupt_mixes',
        columns='threshold',
        aggfunc='mean'
    )
    
    im1 = ax4.imshow(pivot_data.values, cmap='RdYlGn', aspect='auto')
    ax4.set_xticks(range(len(pivot_data.columns)))
    ax4.set_yticks(range(len(pivot_data.index)))
    ax4.set_xticklabels(pivot_data.columns, fontsize=8)
    ax4.set_yticklabels(pivot_data.index, fontsize=6)
    ax4.set_xlabel('Pool Threshold', fontsize=10, fontweight='bold')
    ax4.set_ylabel('Corrupt Mixes', fontsize=10, fontweight='bold')
    ax4.set_title('Heatmap: Corruption×Threshold\n(avg over hops)', fontsize=11, fontweight='bold')
    plt.colorbar(im1, ax=ax4, label='Mean Entropy')
    
    # ============ PLOT 3: Heatmap - Corruption vs N-Hops (averaged over threshold) ============
    ax5 = plt.subplot(2, 3, 3)
    pivot_data2 = valid_df.pivot_table(
        values='mean_entropy',
        index='corrupt_mixes',
        columns='n_hops',
        aggfunc='mean'
    )
    
    im2 = ax5.imshow(pivot_data2.values, cmap='RdYlGn', aspect='auto')
    ax5.set_xticks(range(len(pivot_data2.columns)))
    ax5.set_yticks(range(len(pivot_data2.index)))
    ax5.set_xticklabels(pivot_data2.columns, fontsize=8)
    ax5.set_yticklabels(pivot_data2.index, fontsize=6)
    ax5.set_xlabel('Number of Hops', fontsize=10, fontweight='bold')
    ax5.set_ylabel('Corrupt Mixes', fontsize=10, fontweight='bold')
    ax5.set_title('Heatmap: Corruption×Hops\n(avg over thresholds)', fontsize=11, fontweight='bold')
    plt.colorbar(im2, ax=ax5, label='Mean Entropy')
 
    # ============ PLOT 4: Standard Deviation Analysis ============
    ax7 = plt.subplot(2, 3, 4)
    for threshold in threshold_values:
        subset = valid_df[valid_df['threshold'] == threshold]
        if len(subset) > 0:
            ax7.plot(subset['corrupt_mixes'], subset['std_mean'], 
                    'o-', label=f'Threshold={threshold}', linewidth=2, markersize=5)
    
    ax7.set_xlabel('Corrupt Mixes', fontsize=10, fontweight='bold')
    ax7.set_ylabel('Std Dev of Entropy', fontsize=10, fontweight='bold')
    ax7.set_title('Entropy Variability vs Corruption', fontsize=11, fontweight='bold')
    ax7.legend(fontsize=8)
    ax7.grid(True, alpha=0.3)
    
    # ============ PLOT 5: Best Configuration Analysis ============
    ax8 = plt.subplot(2, 3, 5)
    
    # Find best configuration for each corruption level
    best_configs = []
    for corrupt in corrupt_values:
        subset = valid_df[valid_df['corrupt_mixes'] == corrupt]
        if len(subset) > 0:
            best_idx = subset['mean_entropy'].idxmax()
            best_row = subset.loc[best_idx]
            best_configs.append({
                'corrupt': corrupt,
                'best_entropy': best_row['mean_entropy'],
                'best_threshold': best_row['threshold'],
                'best_hops': best_row['n_hops']
            })
    
    if best_configs:
        best_df = pd.DataFrame(best_configs)
        ax8.plot(best_df['corrupt'], best_df['best_entropy'], 'go-', linewidth=2.5, markersize=8)
        ax8.set_xlabel('Corrupt Mixes', fontsize=10, fontweight='bold')
        ax8.set_ylabel('Maximum Achievable Entropy', fontsize=10, fontweight='bold')
        ax8.set_title('Best Entropy vs Corruption\n(optimized threshold & hops)', fontsize=11, fontweight='bold')
        ax8.grid(True, alpha=0.3)
        
        # Annotate some key points
        for i in range(0, len(best_df), 4):  # Annotate every 4th point
            row = best_df.iloc[i]
            ax8.annotate(f"T={int(row['best_threshold'])},H={int(row['best_hops'])}", 
                        xy=(row['corrupt'], row['best_entropy']),
                        xytext=(5, 5), textcoords='offset points', fontsize=7)
    
    # ============ PLOT 6: 3D Surface Plot Visualization ============
    ax9 = plt.subplot(2, 3, 6, projection='3d')
    
    # Sample configurations for cleaner 3D visualization
    sample_df = valid_df[valid_df['n_hops'] == 5]  # Fix n_hops for 3D visualization
    
    if len(sample_df) > 0:
        X = sample_df['corrupt_mixes'].values
        Y = sample_df['threshold'].values
        Z = sample_df['mean_entropy'].values
        
        ax9.scatter(X, Y, Z, c=Z, cmap='viridis', s=50, alpha=0.6)
        ax9.set_xlabel('Corrupt Mixes', fontsize=9)
        ax9.set_ylabel('Threshold', fontsize=9)
        ax9.set_zlabel('Entropy', fontsize=9)
        ax9.set_title('3D: Corruption×Threshold×Entropy\n(n_hops=5)', fontsize=10, fontweight='bold')
    
    # Overall title
    plt.suptitle('Grid Topology (10x10): Corrupted Nodes Analysis\n' + 
                 f'20 runs per configuration',
                 fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    # Save plots
    plot_filename = f'grid_corrupt_mix_analysis_{timestamp}.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.savefig(plot_filename.replace('.png', '.pdf'), bbox_inches='tight')
    
    print(f"[ANALYSIS] Plots saved as '{plot_filename}'")
    plt.show()

def print_parameter_sweep_summary(df):
    """Print comprehensive summary statistics"""
    
    valid_df = df[df['mean_entropy'] > 0]
    
    if len(valid_df) == 0:
        print("[ERROR] No valid results to summarize")
        return
    
    print("\n" + "="*100)
    print("PARAMETER SWEEP SUMMARY")
    print("="*100)
    
    # Overall statistics
    print("\n[OVERALL STATISTICS]")
    print(f"Total configurations tested: {len(df)}")
    print(f"Successful configurations: {len(valid_df)}")
    print(f"Failed configurations: {len(df) - len(valid_df)}")
    print(f"Overall entropy range: {valid_df['mean_entropy'].min():.4f} - {valid_df['mean_entropy'].max():.4f} bits")
    print(f"Overall average entropy: {valid_df['mean_entropy'].mean():.4f} bits")
    
    # Best configuration overall
    best_idx = valid_df['mean_entropy'].idxmax()
    best_config = valid_df.loc[best_idx]
    
    print("\n[BEST CONFIGURATION OVERALL]")
    print(f"  Corrupt mixes: {int(best_config['corrupt_mixes'])}")
    print(f"  Threshold: {int(best_config['threshold'])}")
    print(f"  N-hops: {int(best_config['n_hops'])}")
    print(f"  Mean Entropy: {best_config['mean_entropy']:.4f} ± {best_config['std_mean']:.4f} bits")
    
    # Worst configuration
    worst_idx = valid_df['mean_entropy'].idxmin()
    worst_config = valid_df.loc[worst_idx]
    
    print("\n[WORST CONFIGURATION]")
    print(f"  Corrupt mixes: {int(worst_config['corrupt_mixes'])}")
    print(f"  Threshold: {int(worst_config['threshold'])}")
    print(f"  N-hops: {int(worst_config['n_hops'])}")
    print(f"  Mean Entropy: {worst_config['mean_entropy']:.4f} ± {worst_config['std_mean']:.4f} bits")
    
    # Impact of each parameter
    print("\n[PARAMETER IMPACT ANALYSIS]")
    
    print("\n1. Impact of Corruption Level:")
    for corrupt in sorted(valid_df['corrupt_mixes'].unique()):
        subset = valid_df[valid_df['corrupt_mixes'] == corrupt]
        print(f"  corrupt={corrupt:3d}: avg entropy = {subset['mean_entropy'].mean():.4f} bits "
              f"(range: {subset['mean_entropy'].min():.4f} - {subset['mean_entropy'].max():.4f})")
    
    print("\n2. Impact of Threshold:")
    for threshold in sorted(valid_df['threshold'].unique()):
        subset = valid_df[valid_df['threshold'] == threshold]
        print(f"  threshold={threshold}: avg entropy = {subset['mean_entropy'].mean():.4f} bits "
              f"(range: {subset['mean_entropy'].min():.4f} - {subset['mean_entropy'].max():.4f})")
    
    print("\n3. Impact of N-Hops:")
    for hops in sorted(valid_df['n_hops'].unique()):
        subset = valid_df[valid_df['n_hops'] == hops]
        print(f"  n_hops={hops}: avg entropy = {subset['mean_entropy'].mean():.4f} bits "
              f"(range: {subset['mean_entropy'].min():.4f} - {subset['mean_entropy'].max():.4f})")
    
    # Top 10 configurations
    print("\n[TOP 10 CONFIGURATIONS]")
    top10 = valid_df.nlargest(10, 'mean_entropy')
    print(f"{'Rank':<5} {'Corrupt':<8} {'Threshold':<10} {'Hops':<6} {'Entropy':<12} {'Std Dev':<10}")
    print("-" * 70)
    for i, (idx, row) in enumerate(top10.iterrows(), 1):
        print(f"{i:<5} {int(row['corrupt_mixes']):<8} {int(row['threshold']):<10} "
              f"{int(row['n_hops']):<6} {row['mean_entropy']:<12.4f} {row['std_mean']:<10.4f}")
    
    print("="*100)

if __name__ == "__main__":
    print("[PARAMETER SWEEP] Starting comprehensive parameter sweep analysis")
    print("[INFO] This will test all combinations of:")
    print("  - corrupt_mixes: 10, 15, 20, ..., 100")
    print("  - threshold: 70, 75, 80")
    print("  - n_hops: 4, 5, 6")
    print("[INFO] Each combination runs 20 times for statistical significance")
    print("[WARNING] This will take several hours to complete!")
    
    response = input("\nProceed with full analysis? (y/n): ")
    
    if response.lower() != 'y':
        print("Analysis cancelled.")
        exit(0)
    
    start_time = time.time()
    
    try:
        df_results = run_parameter_sweep(runs_per_config=20)
        
        total_time = time.time() - start_time
        print(f"\n[COMPLETE] Total runtime: {total_time/3600:.2f} hours ({total_time/60:.1f} minutes)")
        print(f"[SUCCESS] Results saved and plots generated")
        
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Analysis stopped by user")
    except Exception as e:
        print(f"\n[ERROR] Analysis failed: {e}")
        import traceback
        traceback.print_exc()