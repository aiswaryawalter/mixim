import configparser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from main import main
import os
import time
import statistics

def run_simulation_with_threshold(threshold_value, runs=10):
    """Run simulation `runs` times with a given threshold and return averaged stats."""
    print(f"[ANALYSIS] Starting simulations for threshold = {threshold_value} (runs={runs})")
    config = configparser.ConfigParser()
    config.read('ConfigFile.ini')

    # Ensure MIXING section exists
    if 'MIXING' not in config:
        config.add_section('MIXING')
    config.set('MIXING', 'threshold', str(threshold_value))
    # ensure no external clients when mixes act as clients
    if 'DEFAULT' not in config:
        config['DEFAULT'] = {}
    config['DEFAULT']['n_clients'] = '0'

    temp_config_name = f'temp_config_thresh_{threshold_value}.ini'
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
            print(f"[RUN] threshold={threshold_value} run {run_idx}/{runs}")
            start_run = time.time()
            result = main(3)  # keep calling existing main; adjust if your main signature differs
            run_time = time.time() - start_run

            # robust parsing of result (support list/tuple/dict)
            mean_entropy = 0; median_entropy = 0; q25_entropy = 0
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

    # compute averages and stddevs
    avg_mean = statistics.mean(mean_list) if mean_list else 0
    avg_median = statistics.mean(median_list) if median_list else 0
    avg_q25 = statistics.mean(q25_list) if q25_list else 0
    avg_time = statistics.mean(sim_time_list) if sim_time_list else 0

    return {
        'threshold': threshold_value,
        'mean_entropy': avg_mean,
        'median_entropy': avg_median,
        'q25_entropy': avg_q25,
        'simulation_time': avg_time,
        'raw_means': mean_list,
        'raw_medians': median_list,
        'raw_q25s': q25_list,
        'raw_times': sim_time_list
    }

def create_entropy_statistics_plot(runs_per_threshold=10):
    """Main analysis: run each pool size `runs_per_threshold` times and plot averaged stats."""
    # ...existing code...
    pool_sizes = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

    results = []
    for i, threshold in enumerate(pool_sizes, 1):
        print(f"\n[PROGRESS] Running threshold {threshold} ({i}/{len(pool_sizes)})")
        start_time = time.time()
        try:
            res = run_simulation_with_threshold(threshold, runs=runs_per_threshold)
            results.append({
                'threshold': res['threshold'],
                'mean_entropy': res['mean_entropy'],
                'median_entropy': res['median_entropy'],
                'q25_entropy': res['q25_entropy'],
                'simulation_time': res['simulation_time']
            })
            print(f"[PROGRESS] threshold={threshold} avg mean={res['mean_entropy']:.4f} time={res['simulation_time']:.2f}s")
        except Exception as e:
            print(f"[ERROR] threshold {threshold} failed: {e}")
            results.append({'threshold': threshold,'mean_entropy':0,'median_entropy':0,'q25_entropy':0,'simulation_time':0})

    df = pd.DataFrame(results)
    df.to_csv('entropy_statistics_analysis.csv', index=False)

    print(f"\n[ANALYSIS] Raw data saved to 'entropy_statistics_analysis.csv'")
    
    # Filter out failed simulations (mean_entropy = 0)
    valid_df = df[df['mean_entropy'] > 0]
    
    if len(valid_df) == 0:
        print("[WARNING] No valid simulation results to plot!")
        return
    
    # Create the single plot
    plt.figure(figsize=(12, 8))
    
    # Plot the three entropy measures
    plt.plot(valid_df['threshold'], valid_df['mean_entropy'], 'bo-', 
             label='Mean Entropy', linewidth=2.5, markersize=8)
    plt.plot(valid_df['threshold'], valid_df['median_entropy'], 'go-', 
             label='Median Entropy', linewidth=2.5, markersize=8)
    plt.plot(valid_df['threshold'], valid_df['q25_entropy'], 'ro-', 
             label='Q25 Entropy', linewidth=2.5, markersize=8)
    
    # Mark maximum entropy points
    max_mean_idx = valid_df['mean_entropy'].idxmax()
    max_mean_threshold = valid_df.loc[max_mean_idx, 'threshold']
    max_mean_value = valid_df.loc[max_mean_idx, 'mean_entropy']
    
    
    # Annotate maximum points
    plt.annotate(f'Max Mean: {max_mean_value:.3f}\n@{max_mean_threshold}', 
                 xy=(max_mean_threshold, max_mean_value),
                 xytext=(10, 20), textcoords='offset points',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', color='blue'))
    
  
    # Formatting
    plt.xlabel('Pool Size (Threshold)', fontsize=14, fontweight='bold')
    plt.ylabel('Entropy (bits)', fontsize=14, fontweight='bold')
    plt.title('Grid Topology: Entropy Statistics vs Pool Size\n(Mix-as-Client Mode)', 
              fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11, loc='lower right')
    
    # Set x-axis ticks
    plt.xticks(valid_df['threshold'], rotation=45)
    
    # Add some padding to y-axis
    y_min, y_max = plt.ylim()
    plt.ylim(y_min - 0.1, y_max + 0.2)
    
    plt.tight_layout()
    plt.savefig('entropy_statistics_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('entropy_statistics_comparison.pdf', bbox_inches='tight')
    plt.show()
    
    print(f"[ANALYSIS] Plot saved as 'entropy_statistics_comparison.png' and '.pdf'")
    
    # Print summary of maximum points
    print("\n" + "="*70)
    print("ENTROPY STATISTICS SUMMARY")
    print("="*70)
    print(f"Maximum Mean Entropy:   {max_mean_value:.4f} bits at pool size {max_mean_threshold}")
    print("="*70)

if __name__ == "__main__":
    print("[ENTROPY STATISTICS ANALYSIS] Starting Pool Size vs Entropy Statistics Analysis")
    print("[INFO] This will create a single plot showing Mean, Median, and Q25 entropy trends")
    print("[INFO] Expected runtime: 10-20 minutes depending on your system\n")
    
    start_time = time.time()
    
    try:
        create_entropy_statistics_plot()
        
        total_time = time.time() - start_time
        print(f"\n[ANALYSIS] Complete! Total runtime: {total_time:.2f} seconds")
        print(f"[ANALYSIS] Results saved to:")
        print(f"  - entropy_statistics_analysis.csv (raw data)")
        print(f"  - entropy_statistics_comparison.png (plot)")
        print(f"  - entropy_statistics_comparison.pdf (plot)")
        
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Analysis stopped by user")
    except Exception as e:
        print(f"\n[ERROR] Analysis failed: {e}")
        import traceback
        traceback.print_exc()