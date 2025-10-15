import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob

def load_and_analyze_client_impact():
    """
    Analyze how number of clients affects anonymity metrics for different batch sizes
    """
    # Define the file patterns for each batch size
    batch_size_patterns = {
        3: {
            10: "*10client-3batch*.csv",
            20: "*20client-3batch*.csv", 
            30: "*30client-3batch*.csv"
        },
        4: {
            10: "*10client-4batch*.csv",
            20: "*20client-4batch*.csv",
            30: "*30client-4batch*.csv"
        },
        5: {
            10: "*10client-5batch*.csv",
            20: "*20client-5batch*.csv",
            30: "*30client-5batch*.csv"
        }
    }
    
    logs_folder = Path("final-logs")
    results = {}
    
    # Process each batch size
    for batch_size, client_patterns in batch_size_patterns.items():
        print(f"\nProcessing Batch Size {batch_size}...")
        
        batch_results = {}
        
        # Process each client count for this batch size
        for n_clients, pattern in client_patterns.items():
            print(f"  Looking for files matching: {pattern}")
            
            # Find matching files
            matching_files = list(logs_folder.glob(pattern))
            
            if not matching_files:
                print(f"    WARNING: No files found for {n_clients} clients, batch size {batch_size}")
                continue
            
            # Use the first matching file
            csv_file = matching_files[0]
            print(f"    Found: {csv_file.name}")
            
            # Load and process the file
            df = pd.read_csv(csv_file)
            
            # Get data from the largest window index only
            max_window = df['window_index'].max()
            df_max_window = df[df['window_index'] == max_window]
            
            print(f"    Max window index: {max_window}")
            print(f"    Records at max window: {len(df_max_window)}")
            
            if len(df_max_window) == 0:
                print(f"    WARNING: No data at max window for {n_clients} clients")
                continue
            
            # Calculate the three metrics
            # 1. Number of uniquely identified batches
            uniquely_identified = len(df_max_window[
                (df_max_window['correct_batch_prob'] == 1.0) & 
                (df_max_window['correct_batch_is_highest'] == True)
            ])
            
            # 2. Average anonymity set size
            avg_anonymity_size = df_max_window['anonymity_set_size'].mean()
            
            # 3. Accuracy percentage
            accuracy_count = len(df_max_window[df_max_window['correct_batch_is_highest'] == True])
            total_batches = len(df_max_window)
            accuracy_percentage = (accuracy_count / total_batches) * 100 if total_batches > 0 else 0
            
            # Store results
            batch_results[n_clients] = {
                'uniquely_identified': uniquely_identified,
                'avg_anonymity_size': avg_anonymity_size,
                'accuracy_percentage': accuracy_percentage,
                'total_batches': total_batches,
                'filename': csv_file.name
            }
            
            print(f"    Metrics: Unique={uniquely_identified}, Avg_Anon={avg_anonymity_size:.2f}, Accuracy={accuracy_percentage:.1f}%")
        
        results[batch_size] = batch_results
    
    return results

def create_client_impact_plots(results):
    """
    Create line plots showing how metrics change with number of clients for each batch size
    """
    # Set up the plotting style
    plt.style.use('seaborn-v0_8')
    
    # Create subplots for each batch size
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle('Impact of Number of Clients on Anonymity Metrics', 
                 fontsize=16, fontweight='bold')
    
    colors = ['#e74c3c', '#3498db', '#2ecc71']  # Red, Blue, Green
    metrics = ['uniquely_identified', 'avg_anonymity_size', 'accuracy_percentage']
    metric_titles = [
        'Number of Uniquely Identified Batches',
        'Average Anonymity Set Size', 
        'Accuracy Percentage (%)'
    ]
    
    batch_sizes = [3, 4, 5]
    
    for batch_idx, batch_size in enumerate(batch_sizes):
        print(f"\nCreating plots for Batch Size {batch_size}")
        
        if batch_size not in results or not results[batch_size]:
            print(f"  No data available for batch size {batch_size}")
            continue
        
        batch_data = results[batch_size]
        client_counts = sorted(batch_data.keys())
        
        if len(client_counts) < 2:
            print(f"  Not enough client counts for batch size {batch_size}")
            continue
        
        # Plot each metric
        for metric_idx, (metric_key, metric_title) in enumerate(zip(metrics, metric_titles)):
            ax = axes[metric_idx, batch_idx]
            
            # Extract values for this metric
            values = [batch_data[n_clients][metric_key] for n_clients in client_counts]
            
            # Create line plot
            ax.plot(client_counts, values, marker='o', linewidth=3, markersize=8, 
                   color=colors[metric_idx], alpha=0.8, markerfacecolor='white', 
                   markeredgewidth=2, markeredgecolor=colors[metric_idx])
            
            # Customize the plot
            ax.set_title(f'{metric_title}\n(Batch Size {batch_size})', fontweight='bold', fontsize=11)
            ax.set_xlabel('Number of Clients', fontweight='bold')
            ax.set_ylabel(metric_title.replace(' (%)', ''), fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.set_xticks(client_counts)
            
            # Set y-axis limits based on metric type
            if metric_key == 'accuracy_percentage':
                ax.set_ylim(0, 100)
            elif metric_key == 'uniquely_identified':
                ax.set_ylim(bottom=0)
            else:  # avg_anonymity_size
                ax.set_ylim(bottom=0)
            
            # Add value annotations on points
            for x, y in zip(client_counts, values):
                if metric_key == 'accuracy_percentage':
                    ax.annotate(f'{y:.1f}%', (x, y), textcoords="offset points", 
                               xytext=(0,10), ha='center', fontweight='bold', fontsize=9)
                elif metric_key == 'avg_anonymity_size':
                    ax.annotate(f'{y:.1f}', (x, y), textcoords="offset points", 
                               xytext=(0,10), ha='center', fontweight='bold', fontsize=9)
                else:  # uniquely_identified
                    ax.annotate(f'{int(y)}', (x, y), textcoords="offset points", 
                               xytext=(0,10), ha='center', fontweight='bold', fontsize=9)
            
            # # Add trend analysis
            # if len(client_counts) >= 2:
            #     # Calculate trend
            #     if values[-1] > values[0]:
            #         trend = "↗"
            #     elif values[-1] < values[0]:
            #         trend = "↘"
            #     else:
            #         trend = "→"
                
            #     ax.text(0.02, 0.95, f'Trend: {trend}', transform=ax.transAxes, 
            #            verticalalignment='top', fontweight='bold', fontsize=10,
            #            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    
    # Save the plots
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)
    
    output_path = diagrams_folder / 'client_impact_analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved comprehensive plot: {output_path}")
    
    plt.show()
    
    return fig

def create_summary_table(results):
    """
    Create a summary table of all results
    """
    print("\n" + "="*80)
    print("SUMMARY TABLE: Impact of Number of Clients on Anonymity Metrics")
    print("="*80)
    
    # Create summary data
    summary_data = []
    
    for batch_size in [3, 4, 5]:
        if batch_size not in results:
            continue
            
        batch_data = results[batch_size]
        
        for n_clients in sorted(batch_data.keys()):
            data = batch_data[n_clients]
            summary_data.append({
                'Batch Size': batch_size,
                'Clients': n_clients,
                'Unique Identified': data['uniquely_identified'],
                'Avg Anonymity Size': data['avg_anonymity_size'],
                'Accuracy (%)': data['accuracy_percentage'],
                'Total Batches': data['total_batches'],
                'File': data['filename']
            })
    
    # Convert to DataFrame and display
    df_summary = pd.DataFrame(summary_data)
    
    if not df_summary.empty:
        print(df_summary.to_string(index=False, float_format='%.2f'))
        
        # Save summary table
        diagrams_folder = Path('diagrams')
        diagrams_folder.mkdir(exist_ok=True)
        
        summary_path = diagrams_folder / 'client_impact_summary.csv'
        df_summary.to_csv(summary_path, index=False)
        print(f"\nSummary table saved: {summary_path}")
    
    return df_summary

def analyze_trends(results):
    """
    Analyze trends for each batch size
    """
    print("\n" + "="*60)
    print("TREND ANALYSIS")
    print("="*60)
    
    for batch_size in [3, 4, 5]:
        if batch_size not in results or len(results[batch_size]) < 2:
            continue
        
        print(f"\nBatch Size {batch_size}:")
        batch_data = results[batch_size]
        client_counts = sorted(batch_data.keys())
        
        # Analyze each metric
        metrics = [
            ('uniquely_identified', 'Uniquely Identified Batches'),
            ('avg_anonymity_size', 'Average Anonymity Set Size'),
            ('accuracy_percentage', 'Accuracy Percentage')
        ]
        
        for metric_key, metric_name in metrics:
            values = [batch_data[n_clients][metric_key] for n_clients in client_counts]
            
            if len(values) >= 2:
                change = values[-1] - values[0]
                percent_change = (change / values[0] * 100) if values[0] != 0 else 0
                
                if abs(percent_change) < 5:
                    trend_desc = "relatively stable"
                elif percent_change > 0:
                    trend_desc = f"increases by {percent_change:.1f}%"
                else:
                    trend_desc = f"decreases by {abs(percent_change):.1f}%"
                
                print(f"  {metric_name}: {trend_desc}")
                print(f"    Values: {' → '.join([f'{v:.1f}' for v in values])}")

# Main execution
if __name__ == "__main__":
    print("Analyzing Client Impact on Anonymity Metrics...")
    print("="*50)
    
    # Load and analyze data
    results = load_and_analyze_client_impact()
    
    if not results:
        print("No data found to analyze!")
        exit(1)
    
    # Create comprehensive plots
    create_client_impact_plots(results)
    
    # Create summary table
    create_summary_table(results)
    
    # Analyze trends
    analyze_trends(results)
    
    print(f"\nAnalysis complete!")
    print("Generated files:")
    print("  - diagrams/client_impact_analysis.png (main plots)")
    print("  - diagrams/client_impact_summary.csv (summary table)")