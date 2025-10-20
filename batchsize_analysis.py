import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob

def load_and_analyze_batch_size_impact():
    """
    Analyze how batch size affects anonymity metrics for different client counts
    """
    # Define the file patterns for each client count
    client_count_patterns = {
        10: {
            3: "*10client-3batch*.csv",
            4: "*10client-4batch*.csv",
            5: "*10client-5batch*.csv"
        },
        20: {
            3: "*20client-3batch*.csv",
            4: "*20client-4batch*.csv", 
            5: "*20client-5batch*.csv"
        },
        30: {
            3: "*30client-3batch*.csv",
            4: "*30client-4batch*.csv",
            5: "*30client-5batch*.csv"
        }
    }
    
    logs_folder = Path("final-logs")
    results = {}
    
    # Process each client count
    for n_clients, batch_patterns in client_count_patterns.items():
        print(f"\nProcessing Client Count {n_clients}...")
        
        client_results = {}
        
        # Process each batch size for this client count
        for batch_size, pattern in batch_patterns.items():
            print(f"  Looking for files matching: {pattern}")
            
            # Find matching files
            matching_files = list(logs_folder.glob(pattern))
            
            if not matching_files:
                print(f"    WARNING: No files found for batch size {batch_size}, {n_clients} clients")
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
                print(f"    WARNING: No data at max window for batch size {batch_size}")
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
            client_results[batch_size] = {
                'uniquely_identified': uniquely_identified,
                'avg_anonymity_size': avg_anonymity_size,
                'accuracy_percentage': accuracy_percentage,
                'total_batches': total_batches,
                'filename': csv_file.name
            }
            
            print(f"    Metrics: Unique={uniquely_identified}, Avg_Anon={avg_anonymity_size:.2f}, Accuracy={accuracy_percentage:.1f}%")
        
        results[n_clients] = client_results
    
    return results

def create_batch_size_impact_plots(results):
    """
    Create line plots showing how metrics change with batch size for each client count
    """
    # Set up the plotting style
    plt.style.use('seaborn-v0_8')
    
    # Create subplots for each client count
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Impact of Batch Size on Anonymity Metrics', 
                 fontsize=16, fontweight='bold')
    
    colors = ['#e74c3c', '#3498db', '#2ecc71']  # Red, Blue, Green
    client_counts = [10, 20, 30]
    client_labels = ['10 Clients', '20 Clients', '30 Clients']
    metrics = ['uniquely_identified', 'avg_anonymity_size', 'accuracy_percentage']
    metric_titles = [
        'Number of Uniquely Identified Batches',
        'Average Anonymity Set Size', 
        'Accuracy Percentage (%)'
    ]

    # Plot each metric
    for metric_idx, (metric_key, metric_title) in enumerate(zip(metrics, metric_titles)):
        ax = axes[metric_idx]
        
        # Plot line for each client count
        for client_idx, n_clients in enumerate(client_counts):
            if n_clients not in results or not results[n_clients]:
                print(f"  No data available for client count {n_clients}")
                continue
            
            client_data = results[n_clients]
            batch_sizes = sorted(client_data.keys())
            
            if len(batch_sizes) < 2:
                print(f"  Not enough batch sizes for client count {n_clients}")
                continue
            
            # Extract values for this metric and client count
            values = [client_data[batch_size][metric_key] for batch_size in batch_sizes]
            
            # Create line plot for this client count
            ax.plot(batch_sizes, values, marker='o', linewidth=3, markersize=8, 
                   color=colors[client_idx], alpha=0.8, markerfacecolor='white', 
                   markeredgewidth=2, markeredgecolor=colors[client_idx],
                   label=client_labels[client_idx])
            
            # Add value annotations on points
            for x, y in zip(batch_sizes, values):
                if metric_key == 'accuracy_percentage':
                    ax.annotate(f'{y:.1f}%', (x, y), textcoords="offset points", 
                               xytext=(0,10), ha='center', fontweight='bold', fontsize=8,
                               color=colors[client_idx])
                elif metric_key == 'avg_anonymity_size':
                    ax.annotate(f'{y:.1f}', (x, y), textcoords="offset points", 
                               xytext=(0,10), ha='center', fontweight='bold', fontsize=8,
                               color=colors[client_idx])
                else:  # uniquely_identified
                    ax.annotate(f'{int(y)}', (x, y), textcoords="offset points", 
                               xytext=(0,10), ha='center', fontweight='bold', fontsize=8,
                               color=colors[client_idx])
        
        # Customize the plot
        ax.set_title(f'{metric_title}', fontweight='bold', fontsize=12)
        ax.set_xlabel('Batch Size', fontweight='bold', fontsize=11)
        ax.set_ylabel(metric_title.replace(' (%)', ''), fontweight='bold', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xticks([3, 4, 5])  # Ensure all batch sizes are shown
        
        # Set y-axis limits based on metric type
        if metric_key == 'accuracy_percentage':
            ax.set_ylim(0, 100)
        elif metric_key == 'uniquely_identified':
            ax.set_ylim(bottom=0)
        else:  # avg_anonymity_size
            ax.set_ylim(bottom=0)
        
        # Add legend
        ax.legend(loc='best', fontsize=10, framealpha=0.9)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    
    # Save the plots
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)
    
    output_path = diagrams_folder / 'batch_size_impact_combined.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved combined plot: {output_path}")
        
    plt.show()
    
    return fig

def create_summary_table(results):
    """
    Create a summary table of all results
    """
    print("\n" + "="*80)
    print("SUMMARY TABLE: Impact of Batch Size on Anonymity Metrics")
    print("="*80)
    
    # Create summary data
    summary_data = []
    
    for n_clients in [10, 20, 30]:
        if n_clients not in results:
            continue
            
        client_data = results[n_clients]
        
        for batch_size in sorted(client_data.keys()):
            data = client_data[batch_size]
            summary_data.append({
                'Clients': n_clients,
                'Batch Size': batch_size,
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
        
        summary_path = diagrams_folder / 'batch_size_impact_summary.csv'
        df_summary.to_csv(summary_path, index=False)
        print(f"\nSummary table saved: {summary_path}")
    
    return df_summary

def analyze_trends(results):
    """
    Analyze trends for each client count
    """
    print("\n" + "="*60)
    print("TREND ANALYSIS")
    print("="*60)
    
    for n_clients in [10, 20, 30]:
        if n_clients not in results or len(results[n_clients]) < 2:
            continue
        
        print(f"\n{n_clients} Clients:")
        client_data = results[n_clients]
        batch_sizes = sorted(client_data.keys())
        
        # Analyze each metric
        metrics = [
            ('uniquely_identified', 'Uniquely Identified Batches'),
            ('avg_anonymity_size', 'Average Anonymity Set Size'),
            ('accuracy_percentage', 'Accuracy Percentage')
        ]
        
        for metric_key, metric_name in metrics:
            values = [client_data[batch_size][metric_key] for batch_size in batch_sizes]
            
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

def create_individual_client_plots(results):
    """
    Create individual plots for each client count (3 separate diagrams)
    """
    diagrams_folder = Path('diagrams')
    diagrams_folder.mkdir(exist_ok=True)
    
    for n_clients in [10, 20, 30]:
        if n_clients not in results or not results[n_clients]:
            continue
        
        client_data = results[n_clients]
        batch_sizes = sorted(client_data.keys())
        
        if len(batch_sizes) < 2:
            continue
        
        # Create individual plot for this client count
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f'Impact of Batch Size on Anonymity Metrics\n{n_clients} Clients (Data from Largest Window Index)', 
                     fontsize=14, fontweight='bold')
        
        colors = ['#e74c3c', '#3498db', '#2ecc71']
        metrics = ['uniquely_identified', 'avg_anonymity_size', 'accuracy_percentage']
        metric_titles = [
            'Number of Uniquely Identified Batches',
            'Average Anonymity Set Size', 
            'Accuracy Percentage (%)'
        ]
        
        for metric_idx, (metric_key, metric_title) in enumerate(zip(metrics, metric_titles)):
            ax = axes[metric_idx]
            
            # Extract values for this metric
            values = [client_data[batch_size][metric_key] for batch_size in batch_sizes]
            
            # Create line plot
            ax.plot(batch_sizes, values, marker='o', linewidth=4, markersize=10, 
                   color=colors[metric_idx], alpha=0.8, markerfacecolor='white', 
                   markeredgewidth=3, markeredgecolor=colors[metric_idx])
            
            # Customize the plot
            ax.set_title(f'{metric_title}', fontweight='bold', fontsize=12)
            ax.set_xlabel('Batch Size', fontweight='bold', fontsize=11)
            ax.set_ylabel(metric_title.replace(' (%)', ''), fontweight='bold', fontsize=11)
            ax.grid(True, alpha=0.3)
            ax.set_xticks(batch_sizes)
            
            # Set y-axis limits
            if metric_key == 'accuracy_percentage':
                ax.set_ylim(0, 100)
            elif metric_key == 'uniquely_identified':
                ax.set_ylim(bottom=0)
            else:
                ax.set_ylim(bottom=0)
            
            # Add value annotations
            for x, y in zip(batch_sizes, values):
                if metric_key == 'accuracy_percentage':
                    ax.annotate(f'{y:.1f}%', (x, y), textcoords="offset points", 
                               xytext=(0,15), ha='center', fontweight='bold', fontsize=11)
                elif metric_key == 'avg_anonymity_size':
                    ax.annotate(f'{y:.1f}', (x, y), textcoords="offset points", 
                               xytext=(0,15), ha='center', fontweight='bold', fontsize=11)
                else:
                    ax.annotate(f'{int(y)}', (x, y), textcoords="offset points", 
                               xytext=(0,15), ha='center', fontweight='bold', fontsize=11)
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        
        # Save individual plot
        output_path = diagrams_folder / f'batch_size_impact_{n_clients}_clients.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved individual plot: {output_path}")
        
        plt.show()

# Main execution
if __name__ == "__main__":
    print("Analyzing Batch Size Impact on Anonymity Metrics...")
    print("="*50)
    
    # Load and analyze data
    results = load_and_analyze_batch_size_impact()
    
    if not results:
        print("No data found to analyze!")
        exit(1)
    
    # Create comprehensive plots
    create_batch_size_impact_plots(results)
    
    # Create individual plots for each client count
    # create_individual_client_plots(results)
    
    # Create summary table
    create_summary_table(results)
    
    # Analyze trends
    analyze_trends(results)
    
    print(f"\nAnalysis complete!")
    print("Generated files:")
    print("  - diagrams/batch_size_impact_analysis.png (comprehensive 3x3 grid)")
    print("  - diagrams/batch_size_impact_10_clients.png (individual plot)")
    print("  - diagrams/batch_size_impact_20_clients.png (individual plot)")
    print("  - diagrams/batch_size_impact_30_clients.png (individual plot)")
    print("  - diagrams/batch_size_impact_summary.csv (summary table)")