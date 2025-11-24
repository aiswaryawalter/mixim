import configparser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from main import main
import os
import time

def run_simulation_with_threshold(threshold_value):
    """Run simulation with specific threshold value"""
    print(f"[ANALYSIS] Starting simulation with threshold = {threshold_value}")
    
    # Create a temporary config with the specific threshold
    config = configparser.ConfigParser()
    config.read('ConfigFile.ini')
    
    # Update threshold for this run
    config.set('MIXING', 'threshold', str(threshold_value))
    
    # Create a temporary config file
    temp_config_name = f'temp_config_thresh_{threshold_value}.ini'
    with open(temp_config_name, 'w') as configfile:
        config.write(configfile)
    
    try:
        # Temporarily replace the config file
        os.rename('ConfigFile.ini', 'ConfigFile_backup.ini')
        os.rename(temp_config_name, 'ConfigFile.ini')
        
        # Run the simulation
        result = main(3)  # Pass the rate parameter
        
        # Extract the results
        entropy_data = {
            'threshold': threshold_value,
            'entropy_values': result[0] if len(result) > 0 else [],
            'mean_entropy': result[1] if len(result) > 1 else 0,
            'median_entropy': result[2] if len(result) > 2 else 0,
            'q25_entropy': result[3] if len(result) > 3 else 0
        }
        
        print(f"[ANALYSIS] Completed threshold = {threshold_value}, mean entropy = {entropy_data['mean_entropy']:.4f}")
        return entropy_data
        
    finally:
        # Restore original config
        if os.path.exists('ConfigFile.ini'):
            os.remove('ConfigFile.ini')
        if os.path.exists('ConfigFile_backup.ini'):
            os.rename('ConfigFile_backup.ini', 'ConfigFile.ini')
        if os.path.exists(temp_config_name):
            os.remove(temp_config_name)

def create_entropy_statistics_plot():
    """Main analysis function that creates only the entropy statistics plot"""
    
    # Pool sizes to test
    pool_sizes = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    
    print(f"[ANALYSIS] Testing pool sizes: {pool_sizes}")
    print(f"[ANALYSIS] Total simulations to run: {len(pool_sizes)}")
    
    # Store results
    results = []
    
    # Run simulations sequentially
    for i, threshold in enumerate(pool_sizes, 1):
        print(f"\n[PROGRESS] Running simulation {i}/{len(pool_sizes)} (threshold={threshold})")
        start_time = time.time()
        
        try:
            result = run_simulation_with_threshold(threshold)
            result['simulation_time'] = time.time() - start_time
            results.append(result)
            
            print(f"[PROGRESS] Simulation {i} completed in {result['simulation_time']:.2f} seconds")
            
        except Exception as e:
            print(f"[ERROR] Simulation failed for threshold {threshold}: {e}")
            # Add empty result to maintain data structure
            results.append({
                'threshold': threshold,
                'entropy_values': [],
                'mean_entropy': 0,
                'median_entropy': 0,
                'q25_entropy': 0,
                'simulation_time': 0
            })
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Save raw data
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