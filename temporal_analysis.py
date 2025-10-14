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
    
    # Create the visualization
    create_temporal_plot(metrics_df, n_clients, batch_size, filename)
    
    return metrics_df

def create_temporal_plot(metrics_df, n_clients, batch_size, filename):
    """
    Create a line plot showing temporal changes of the three metrics
    """
    # Set up the plot style
    plt.style.use('seaborn-v0_8')
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    fig.suptitle(f'Temporal Analysis: {n_clients} Clients, Batch Size {batch_size}', 
                 fontsize=16, fontweight='bold')
    
    # Color scheme
    colors = ['#e74c3c', '#3498db', '#2ecc71']
    
    # Plot 1: Number of Uniquely Identified Batches
    axes[0].plot(metrics_df['sim_timestamp'], metrics_df['uniquely_identified'], 
                 marker='o', linewidth=2, markersize=4, color=colors[0], alpha=0.8)
    axes[0].set_title('Number of Uniquely Identified Batches Over Time', fontweight='bold')
    axes[0].set_xlabel('Simulation Time')
    axes[0].set_ylabel('Count of Uniquely Identified Batches')
    axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(bottom=0)
    
    # Add statistics annotation
    mean_unique = metrics_df['uniquely_identified'].mean()
    max_unique = metrics_df['uniquely_identified'].max()
    axes[0].text(0.02, 0.95, f'Mean: {mean_unique:.1f}\nMax: {max_unique}', 
                 transform=axes[0].transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 2: Average Anonymity Set Size
    axes[1].plot(metrics_df['sim_timestamp'], metrics_df['avg_anonymity_size'], 
                 marker='s', linewidth=2, markersize=4, color=colors[1], alpha=0.8)
    axes[1].set_title('Average Anonymity Set Size Over Time', fontweight='bold')
    axes[1].set_xlabel('Simulation Time')
    axes[1].set_ylabel('Average Anonymity Set Size')
    axes[1].grid(True, alpha=0.3)
    axes[1].set_ylim(bottom=0)
    
    # Add statistics annotation
    mean_anon = metrics_df['avg_anonymity_size'].mean()
    min_anon = metrics_df['avg_anonymity_size'].min()
    max_anon = metrics_df['avg_anonymity_size'].max()
    axes[1].text(0.02, 0.95, f'Mean: {mean_anon:.1f}\nMin: {min_anon:.1f}\nMax: {max_anon:.1f}', 
                 transform=axes[1].transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 3: Accuracy Percentage
    axes[2].plot(metrics_df['sim_timestamp'], metrics_df['accuracy_percentage'], 
                 marker='^', linewidth=2, markersize=4, color=colors[2], alpha=0.8)
    axes[2].set_title('Accuracy Over Time (% of Correct Highest Probability)', fontweight='bold')
    axes[2].set_xlabel('Simulation Time')
    axes[2].set_ylabel('Accuracy (%)')
    axes[2].grid(True, alpha=0.3)
    axes[2].set_ylim(0, 100)
    
    # Add statistics annotation
    mean_acc = metrics_df['accuracy_percentage'].mean()
    min_acc = metrics_df['accuracy_percentage'].min()
    max_acc = metrics_df['accuracy_percentage'].max()
    axes[2].text(0.02, 0.95, f'Mean: {mean_acc:.1f}%\nMin: {min_acc:.1f}%\nMax: {max_acc:.1f}%', 
                 transform=axes[2].transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    
    # Save the plot
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)  # Create folder if it doesn't exist
    output_path = diagrams_folder / f'temporal_analysis_{filename}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  - Saved plot: {output_path}")
    
    # Show plot (comment out if running in batch)
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
    
    for filename, metrics_df in all_results.items():
        # Extract metadata
        parts = filename.split('-')
        n_clients = int(parts[1].replace('client', ''))
        batch_size = int(parts[2].replace('batch', ''))
        
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
    summary_df.to_csv('temporal_analysis_summary.csv', index=False)
    print(f"Summary saved to: temporal_analysis_summary.csv")
    
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
    # Set the path to your final-logs folder
    logs_folder = "final-logs"  # Change this to your actual path
    
    # Run the analysis
    results = analyze_all_files(logs_folder)
    
    print(f"\nAnalysis complete! Generated {len(results)} temporal analysis plots.")
    print("Each plot shows:")
    print("  - Top panel: Number of uniquely identified batches over time")
    print("  - Middle panel: Average anonymity set size over time") 
    print("  - Bottom panel: Accuracy percentage over time")