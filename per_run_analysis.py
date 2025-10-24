import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns

def analyze_temporal_changes(csv_file_path):
    """
    Analyze temporal changes for a single CSV file
    """
    # Load the data
    df = pd.read_csv(csv_file_path)
    
    # Extract metadata from filename
    filename = Path(csv_file_path).stem
    print(f"Analyzing {filename}...")
    print(f"  - Total records: {len(df)}")
    print(f"  - Window indexes: {df['window_index'].min()} to {df['window_index'].max()}")
    
    # Group by window_index to calculate metrics at each time point
    temporal_metrics = []
    
    for window_idx in sorted(df['window_index'].unique()):
        window_data = df[df['window_index'] == window_idx]
        
        # Get simulation time (should be same for all records in this window)
        sim_time = window_data['sim_timestamp'].iloc[0]
        
        # Metric 1: Number of uniquely identified batches
        # (correct_batch_prob = 1 AND correct_batch_is_highest = True)
        uniquely_identified = len(window_data[
            (window_data['correct_batch_prob'] == 1.0) & 
            (window_data['correct_batch_is_highest'] == True)
        ])
        
        # Metric 2: Average anonymity set size
        avg_anonymity_size = window_data['anonymity_set_size'].mean()
        
        # Metric 3: Accuracy (number of batches with correct_batch_is_highest = True)
        accuracy_count = len(window_data[window_data['correct_batch_is_highest'] == True])
        total_batches = len(window_data)
        accuracy_percentage = (accuracy_count / total_batches) * 100 if total_batches > 0 else 0
        
        temporal_metrics.append({
            'window_index': window_idx,
            'sim_timestamp': sim_time,
            'uniquely_identified': uniquely_identified,
            'avg_anonymity_size': avg_anonymity_size,
            'accuracy_count': accuracy_count,
            'accuracy_percentage': accuracy_percentage,
            'total_batches': total_batches
        })
    
    metrics_df = pd.DataFrame(temporal_metrics)
    return metrics_df, filename

def smooth_data(data, window_size=3):
    """
    Apply sliding window averaging to smooth the data
    """
    if len(data) < window_size:
        return data  # Return original if not enough data points
    
    return data.rolling(window=window_size, center=True, min_periods=1).mean()

def create_combined_20client_4batch_plot(all_results, smoothing_window=3):
    """
    Create a combined temporal plot for all 20-client, 4-batch runs
    """
    # Set up the plot style
    plt.style.use('seaborn-v0_8')
    fig, axes = plt.subplots(3, 1, figsize=(16, 12))
    
    # smoothing_text = f" (Smoothed with {smoothing_window}-point moving average)" if smoothing_window > 1 else ""
    fig.suptitle(f'Temporal Analysis: 20 Clients, Batch Size 4\n(Comparing Different Runs)', 
                 fontsize=16, fontweight='bold')
    
    # Color scheme for different runs
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6']  # Red, Blue, Green, Purple
    run_labels = []
    
    # Sort results by filename for consistent ordering
    sorted_results = sorted(all_results.items())
    
    for idx, (filename, (metrics_df, _)) in enumerate(sorted_results):
        # Extract run identifier from filename
        run_id = str(idx + 1)
        run_labels.append(f'Run {run_id}')
    
    metrics = ['uniquely_identified', 'avg_anonymity_size', 'accuracy_percentage']
    metric_titles = [
        'Number of Uniquely Identified Batches Over Time',
        'Average Anonymity Set Size Over Time',
        'Accuracy Over Time (% that Correct Batch has Highest Probability)'
    ]
    metric_ylabels = [
        'Count of Uniquely Identified Batches',
        'Average Anonymity Set Size',
        'Accuracy (%)'
    ]
    
    # Plot each metric
    for metric_idx, (metric_key, metric_title, ylabel) in enumerate(zip(metrics, metric_titles, metric_ylabels)):
        ax = axes[metric_idx]
        
        # Plot line for each run
        for run_idx, (filename, (metrics_df, _)) in enumerate(sorted_results):
            
            # Apply smoothing if requested
            if smoothing_window > 1:
                smoothed_values = smooth_data(metrics_df[metric_key], smoothing_window)
                # Plot both original (light/transparent) and smoothed (bold) lines
                ax.plot(metrics_df['sim_timestamp'], metrics_df[metric_key], 
                       color=colors[run_idx], alpha=0.3, linewidth=1, linestyle='--')
                ax.plot(metrics_df['sim_timestamp'], smoothed_values, 
                       linewidth=3, markersize=4, 
                       color=colors[run_idx], alpha=0.9,
                       label=run_labels[run_idx])
            else:
                # Original plotting without smoothing
                ax.plot(metrics_df['sim_timestamp'], metrics_df[metric_key], 
                       marker='o', linewidth=2, markersize=3, 
                       color=colors[run_idx], alpha=0.8,
                       label=run_labels[run_idx])
        
        # Customize the plot
        ax.set_title(metric_title, fontweight='bold', fontsize=12)
        ax.set_xlabel('Simulated Time', fontweight='bold')
        ax.set_ylabel(ylabel, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Set y-axis limits
        if metric_key == 'accuracy_percentage':
            ax.set_ylim(0, 100)
        else:
            ax.set_ylim(bottom=0)
        
        # Add legend
        ax.legend(loc='best', fontsize=10, framealpha=0.9)
        
        # Add summary statistics
        stats_text = []
        for run_idx, (filename, (metrics_df, _)) in enumerate(sorted_results):
            if smoothing_window > 1:
                smoothed_values = smooth_data(metrics_df[metric_key], smoothing_window)
                mean_val = smoothed_values.mean()
            else:
                mean_val = metrics_df[metric_key].mean()
                
            if metric_key == 'accuracy_percentage':
                stats_text.append(f'{run_labels[run_idx]}: {mean_val:.1f}%')
            else:
                stats_text.append(f'{run_labels[run_idx]}: {mean_val:.1f}')

        if stats_text:
            legend_text = 'Mean values:\n' if smoothing_window > 1 else 'Mean values:\n'
            ax.text(0.02, 0.95, legend_text + '\n'.join(stats_text), 
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                   fontsize=9)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    
    # Save the plot
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)
    
    output_path = diagrams_folder / f'temporal_20client_4batch_combined_smooth_{smoothing_window}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved combined temporal plot: {output_path}")
    
    plt.show()
    plt.close()

def analyze_20client_4batch_files():
    """
    Analyze all 20-client, 4-batch CSV files
    """
    # Define the file patterns for 20-client, 4-batch runs
    file_patterns = [
        "12hr-20client-4batch-19381883_1.csv",
        "12hr-20client-4batch-19496121_1.csv", 
        "12hr-20client-4batch-19496122_1.csv",
        "12hr-20client-4batch-19496123_1.csv"
    ]
    
    current_dir = Path(".")
    all_results = {}
    
    print("Looking for 20-client, 4-batch files...")
    print("=" * 50)
    
    for pattern in file_patterns:
        # Look for files in current directory and subdirectories
        matching_files = list(current_dir.rglob(pattern))
        
        if matching_files:
            csv_file = matching_files[0]  # Use the first match
            print(f"Found: {csv_file}")
            
            try:
                metrics_df, filename = analyze_temporal_changes(csv_file)
                all_results[filename] = (metrics_df, filename)
                print("  - Analysis completed successfully")
            except Exception as e:
                print(f"  - Error analyzing {csv_file.name}: {e}")
        else:
            print(f"  - File not found: {pattern}")
        
        print("-" * 30)
    
    return all_results

def create_summary_comparison(all_results):
    """
    Create a summary comparison across all runs
    """
    if not all_results:
        return
    
    summary_data = []
    
    for filename, (metrics_df, _) in all_results.items():
        # Calculate final window metrics
        final_metrics = metrics_df.iloc[-1] if len(metrics_df) > 0 else None
        
        if final_metrics is not None:
            summary_data.append({
                'filename': filename,
                'final_window': final_metrics['window_index'],
                'final_unique_identified': final_metrics['uniquely_identified'],
                'final_avg_anonymity_size': final_metrics['avg_anonymity_size'],
                'final_accuracy': final_metrics['accuracy_percentage'],
                'avg_unique_identified': metrics_df['uniquely_identified'].mean(),
                'avg_anonymity_size': metrics_df['avg_anonymity_size'].mean(),
                'avg_accuracy': metrics_df['accuracy_percentage'].mean(),
                'total_windows': len(metrics_df),
                'max_sim_time': metrics_df['sim_timestamp'].max()
            })
    
    summary_df = pd.DataFrame(summary_data)
    
    # Save summary
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)
    
    summary_path = diagrams_folder / 'temporal_20client_4batch_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary saved to: {summary_path}")
    
    # Print summary table
    print("\nSUMMARY OF 20-CLIENT, 4-BATCH RUNS:")
    print("=" * 80)
    print(f"{'Run':<15} {'Windows':<8} {'Max Time':<10} {'Avg Unique':<11} {'Avg Anon Size':<13} {'Avg Accuracy':<12}")
    print("-" * 80)
    
    for _, row in summary_df.iterrows():
        run_name = row['filename'].split('_')[-1] if '_' in row['filename'] else 'Unknown'
        print(f"Run {run_name:<11} {row['total_windows']:<8} {row['max_sim_time']:<10.2f} "
              f"{row['avg_unique_identified']:<11.1f} {row['avg_anonymity_size']:<13.1f} "
              f"{row['avg_accuracy']:<12.1f}%")

# Main execution
if __name__ == "__main__":
    print("Analyzing 20-Client, 4-Batch Temporal Changes...")
    print("=" * 50)
    
    # Analyze all files
    all_results = analyze_20client_4batch_files()
    
    if not all_results:
        print("No data found to analyze!")
        exit(1)
    
    print(f"\nFound {len(all_results)} runs to analyze")
    
    # Create plots with different smoothing levels
    smoothing_options = [5]  # 1 = no smoothing, 3 = 3-point average, 5 = 5-point average
    
    for smoothing in smoothing_options:
        print(f"\nCreating plots with smoothing window = {smoothing}")
        create_combined_20client_4batch_plot(all_results, smoothing_window=smoothing)
    
    # Create summary statistics
    create_summary_comparison(all_results)
    
    print(f"\nAnalysis complete!")
    print("Generated files:")
    print("  - diagrams/temporal_20client_4batch_combined_smooth_1.png (no smoothing)")
    print("  - diagrams/temporal_20client_4batch_combined_smooth_3.png (3-point smoothing)")
    print("  - diagrams/temporal_20client_4batch_combined_smooth_5.png (5-point smoothing)")
    print("  - diagrams/temporal_20client_4batch_summary.csv")
    print("\nEach diagram shows 3 graphs:")
    print("  - Top: Number of uniquely identified batches over time (4 lines for different runs)")
    print("  - Middle: Average anonymity set size over time (4 lines for different runs)")
    print("  - Bottom: Accuracy percentage over time (4 lines for different runs)")