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
    parts = filename.split('-')
    n_clients = int(parts[1].replace('client', ''))
    batch_size = int(parts[2].replace('batch', ''))
    
    print(f"Analyzing {filename}...")
    print(f"  - Clients: {n_clients}, Batch Size: {batch_size}")
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
    return metrics_df, n_clients, batch_size, filename
    

def create_combined_temporal_plots(all_results, smoothing_window=3):
    """
    Create combined temporal plots for each client count
    """
    # Group results by client count
    client_groups = {}
    for filename, data in all_results.items():
        metrics_df, n_clients, batch_size, _ = data
        if n_clients not in client_groups:
            client_groups[n_clients] = {}
        client_groups[n_clients][batch_size] = metrics_df
    
    # Create plots for each client count
    for n_clients in sorted(client_groups.keys()):
        create_client_temporal_plot(client_groups[n_clients], n_clients, smoothing_window)

def smooth_data(data, window_size=3):
    """
    Apply sliding window averaging to smooth the data
    
    Args:
        data: pandas Series or numpy array of values
        window_size: size of the sliding window (default: 3)
    
    Returns:
        smoothed data as pandas Series
    """
    if len(data) < window_size:
        return data  # Return original if not enough data points
    
    return data.rolling(window=window_size, center=True, min_periods=1).mean()

def create_client_temporal_plot(batch_data, n_clients, smoothing_window=3):
    """
    Create a temporal plot for one client count showing all batch sizes
    """
    # Set up the plot style
    plt.style.use('seaborn-v0_8')
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))

    smoothing_text = f" (Smoothed with {smoothing_window}-point moving average)" if smoothing_window > 1 else ""
    fig.suptitle(f'Temporal Analysis: {n_clients} Clients\n(Comparing Different Batch Sizes)', 
                 fontsize=16, fontweight='bold')
    
    # Color scheme for different batch sizes
    colors = ['#e74c3c', '#3498db', '#2ecc71']  # Red, Blue, Green for batch sizes 3, 4, 5
    batch_sizes = sorted(batch_data.keys())
    batch_labels = [f'Batch Size {bs}' for bs in batch_sizes]
    
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
        
        # Plot line for each batch size
        for batch_idx, batch_size in enumerate(batch_sizes):
            if batch_size not in batch_data:
                continue
                
            metrics_df = batch_data[batch_size].copy()

            # Apply smoothing
            if smoothing_window > 1:
                smoothed_values = smooth_data(metrics_df[metric_key], smoothing_window)
                # Plot both original (light/transparent) and smoothed (bold) lines
                ax.plot(metrics_df['sim_timestamp'], metrics_df[metric_key], 
                       color=colors[batch_idx], alpha=0.3, linewidth=1, linestyle='--')
                ax.plot(metrics_df['sim_timestamp'], smoothed_values, 
                       marker='o', linewidth=3, markersize=4, 
                       color=colors[batch_idx], alpha=0.9,
                       label=batch_labels[batch_idx])
            else:
                # Original plotting without smoothing
                ax.plot(metrics_df['sim_timestamp'], metrics_df[metric_key], 
                       marker='o', linewidth=2, markersize=4, 
                       color=colors[batch_idx], alpha=0.8,
                       label=batch_labels[batch_idx])
        
        # Customize the plot
        ax.set_title(metric_title, fontweight='bold')
        ax.set_xlabel('Simulated Time')
        ax.set_ylabel(ylabel)
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
        for batch_idx, batch_size in enumerate(batch_sizes):
            if batch_size not in batch_data:
                continue
            metrics_df = batch_data[batch_size]

            if smoothing_window > 1:
                smoothed_values = smooth_data(metrics_df[metric_key], smoothing_window)
                mean_val = smoothed_values.mean()
            else:
                mean_val = metrics_df[metric_key].mean()
                
            if metric_key == 'accuracy_percentage':
                stats_text.append(f'BS{batch_size}: {mean_val:.1f}%')
            else:
                stats_text.append(f'BS{batch_size}: {mean_val:.1f}')

        if stats_text:
            legend_text = 'Smoothed Mean values:\n' if smoothing_window > 1 else 'Mean values:\n'
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
    
    output_path = diagrams_folder / f'temporal_{smoothing_window}_smoothed_{n_clients}_clients.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved combined temporal plot: {output_path}")
    
    plt.show()
    plt.close()

def analyze_all_files(folder_path):
    """
    Analyze all CSV files in the specified folder
    """
    folder = Path(folder_path)
    csv_files = list(folder.glob("*.csv"))
    
    if not csv_files:
        print(f"No CSV files found in {folder_path}")
        return
    
    print(f"Found {len(csv_files)} CSV files to analyze")
    print("=" * 50)
    
    all_results = {}
    
    for csv_file in sorted(csv_files):
        try:
            metrics_df = analyze_temporal_changes(csv_file)
            all_results[csv_file.name] = metrics_df
            print("  - Analysis completed successfully")
            print("-" * 30)
        except Exception as e:
            print(f"  - Error analyzing {csv_file.name}: {e}")
            print("-" * 30)
    
    # Create summary statistics
    create_summary_comparison(all_results)
    
    return all_results

def create_summary_comparison(all_results):
    """
    Create a summary comparison across all experiments
    """
    if not all_results:
        return
    
    summary_data = []
    
    for filename, data in all_results.items():
        metrics_df, n_clients, batch_size, _ = data
        
        # Calculate final window metrics
        final_metrics = metrics_df.iloc[-1] if len(metrics_df) > 0 else None
        
        if final_metrics is not None:
            summary_data.append({
                'filename': filename,
                'n_clients': n_clients,
                'batch_size': batch_size,
                'final_window': final_metrics['window_index'],
                'final_unique_identified': final_metrics['uniquely_identified'],
                'final_avg_anonymity_size': final_metrics['avg_anonymity_size'],
                'final_accuracy': final_metrics['accuracy_percentage'],
                'avg_unique_identified': metrics_df['uniquely_identified'].mean(),
                'avg_anonymity_size': metrics_df['avg_anonymity_size'].mean(),
                'avg_accuracy': metrics_df['accuracy_percentage'].mean(),
                'total_windows': len(metrics_df)
            })
    
    summary_df = pd.DataFrame(summary_data)
    
    # Save summary
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)
    
    summary_path = diagrams_folder / 'temporal_combined_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary saved to: {summary_path}")
    
    # Print summary table
    print("\nSUMMARY OF ALL EXPERIMENTS:")
    print("=" * 80)
    print(f"{'Clients':<8} {'Batch':<6} {'Windows':<8} {'Avg Unique':<11} {'Avg Anon Size':<13} {'Avg Accuracy':<12}")
    print("-" * 80)
    
    for _, row in summary_df.iterrows():
        print(f"{row['n_clients']:<8} {row['batch_size']:<6} {row['total_windows']:<8} "
              f"{row['avg_unique_identified']:<11.1f} {row['avg_anonymity_size']:<13.1f} "
              f"{row['avg_accuracy']:<12.1f}%")

# Main execution
if __name__ == "__main__":
    print("Analyzing Temporal Changes with Combined Plots...")
    print("=" * 50)
    
    # Set the path to your final-logs folder
    logs_folder = "final-logs"  # Change this to your actual path
    
    # Analyze all files
    all_results = analyze_all_files(logs_folder)
    
    if not all_results:
        print("No data found to analyze!")
        exit(1)
    
    # Create plots with different smoothing levels
    smoothing_options = [5]  # 1 = no smoothing, 3 = 3-point average, 5 = 5-point average
    
    for smoothing in smoothing_options:
        print(f"\nCreating plots with smoothing window = {smoothing}")
        create_combined_temporal_plots(all_results, smoothing_window=smoothing)
    
    # Create summary statistics
    create_summary_comparison(all_results)
    
    print(f"\nAnalysis complete!")
    print("Generated files:")
    print("  - diagrams/temporal_combined_10_clients.png")
    print("  - diagrams/temporal_combined_20_clients.png") 
    print("  - diagrams/temporal_combined_30_clients.png")
    print("  - diagrams/temporal_combined_summary.csv")
    print("\nEach diagram shows 3 graphs:")
    print("  - Top: Number of uniquely identified batches over time (3 lines for batch sizes 3,4,5)")
    print("  - Middle: Average anonymity set size over time (3 lines for batch sizes 3,4,5)")
    print("  - Bottom: Accuracy percentage over time (3 lines for batch sizes 3,4,5)")