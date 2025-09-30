import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import ast
from pathlib import Path

# Set style for better-looking plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def load_and_process_data(csv_file):
    """Load and process the batch logs CSV file"""
    df = pd.read_csv(csv_file)
    
    # Convert string representations of sets/dicts to actual objects
    df['anonymity_set'] = df['anonymity_set'].apply(ast.literal_eval)
    df['batch_prob'] = df['batch_prob'].apply(ast.literal_eval)
    
    # Add derived columns
    df['uniquely_identified'] = df['anonymity_set_size'] == 1
    df['correct_batch_rank'] = df.apply(get_correct_batch_rank, axis=1)
    
    return df

def get_correct_batch_rank(row):
    """Get the rank of the correct batch in the probability distribution"""
    batch_probs = row['batch_prob']
    correct_prob = row['correct_batch_prob']
    
    # Sort probabilities in descending order
    sorted_probs = sorted(batch_probs.values(), reverse=True)
    
    # Find rank (1-indexed)
    for rank, prob in enumerate(sorted_probs, 1):
        if abs(prob - correct_prob) < 1e-10:  # Handle floating point precision
            return rank
    return len(sorted_probs)  # Fallback

def plot_1_fraction_uniquely_identified(df):
    """1. Fraction of Batches Uniquely Identified - Pie Chart"""
    unique_counts = df['uniquely_identified'].value_counts()
    
    plt.figure(figsize=(8, 6))
    labels = ['Not Uniquely Identified', 'Uniquely Identified']
    colors = ['#2ecc71', '#e74c3c']  # Green for good, red for bad
    
    plt.pie(unique_counts.values, labels=labels, autopct='%1.1f%%', 
            colors=colors, startangle=90)
    plt.title('Fraction of Batches Uniquely Identified\n(Lower is Better for Anonymity)', 
              fontsize=14, fontweight='bold')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('logs_analysis/1_fraction_uniquely_identified.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Unique identification rate: {unique_counts.get(True, 0) / len(df) * 100:.1f}%")

def plot_2_anonymity_set_distribution(df):
    """2. Anonymity Set Size Distribution - Histogram"""
    plt.figure(figsize=(10, 6))
    
    plt.hist(df['anonymity_set_size'], bins=range(1, df['anonymity_set_size'].max() + 2), 
             alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('Anonymity Set Size', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.title('Distribution of Anonymity Set Sizes\n(Higher values indicate better anonymity)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Add statistics
    mean_size = df['anonymity_set_size'].mean()
    plt.axvline(mean_size, color='red', linestyle='--', 
                label=f'Mean: {mean_size:.1f}')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('logs_analysis/2_anonymity_set_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_3_correct_batch_probability(df):
    """3. Probability Assigned to the Correct Batch - Scatter/Line Plot"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Scatter plot
    ax1.scatter(df.index, df['correct_batch_prob'], alpha=0.6, s=30)
    ax1.set_xlabel('Batch Index', fontsize=12)
    ax1.set_ylabel('Probability Assigned to Correct Batch', fontsize=12)
    ax1.set_title('Probability Assigned to Correct Batch Over Time\n(Lower is Better for Anonymity)', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label='50% threshold')
    ax1.legend()
    
    # Line plot with moving average
    window_size = max(1, len(df) // 20)  # 5% of data points
    moving_avg = df['correct_batch_prob'].rolling(window=window_size, center=True).mean()
    
    ax2.plot(df.index, df['correct_batch_prob'], alpha=0.3, color='lightblue', label='Individual batches')
    ax2.plot(df.index, moving_avg, color='darkblue', linewidth=2, label=f'Moving average (window={window_size})')
    ax2.set_xlabel('Batch Index', fontsize=12)
    ax2.set_ylabel('Probability Assigned to Correct Batch', fontsize=12)
    ax2.set_title('Correct Batch Probability Trend', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('logs_analysis/3_correct_batch_probability.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_4_accuracy_highest_probability(df):
    """4. Accuracy: Correct Batch Has Highest Probability - Bar Chart"""
    # Calculate accuracy over time windows
    window_size = max(1, len(df) // 10)  # 10 windows
    windows = []
    accuracies = []
    
    for i in range(0, len(df), window_size):
        window_data = df.iloc[i:i+window_size]
        accuracy = window_data['correct_batch_is_highest'].mean()
        windows.append(i + window_size//2)  # Middle of window
        accuracies.append(accuracy)
    
    plt.figure(figsize=(12, 6))
    plt.bar(windows, accuracies, width=window_size*0.8, alpha=0.7, color='coral')
    plt.xlabel('Batch Index (Window Centers)', fontsize=12)
    plt.ylabel('Accuracy (Fraction)', fontsize=12)
    plt.title('Accuracy: Fraction of Batches Where Correct Batch Has Highest Probability\n(Lower is Better for Anonymity)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add overall accuracy line
    overall_accuracy = df['correct_batch_is_highest'].mean()
    plt.axhline(y=overall_accuracy, color='red', linestyle='--', 
                label=f'Overall Accuracy: {overall_accuracy:.1%}')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('logs_analysis/4_accuracy_highest_probability.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_5_temporal_changes(df):
    """5. Temporal Changes - Line Plots"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # 1. Fraction uniquely identified over time
    window_size = max(1, len(df) // 20)
    unique_frac = df['uniquely_identified'].rolling(window=window_size, center=True).mean()
    
    ax1.plot(df['sim_timestamp'], unique_frac, color='red', linewidth=2)
    ax1.set_xlabel('Simulation Time', fontsize=10)
    ax1.set_ylabel('Fraction Uniquely Identified', fontsize=10)
    ax1.set_title('Temporal: Unique Identification Rate', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Average anonymity set size over time
    avg_anon_size = df['anonymity_set_size'].rolling(window=window_size, center=True).mean()
    
    ax2.plot(df['sim_timestamp'], avg_anon_size, color='green', linewidth=2)
    ax2.set_xlabel('Simulation Time', fontsize=10)
    ax2.set_ylabel('Average Anonymity Set Size', fontsize=10)
    ax2.set_title('Temporal: Anonymity Set Size', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # 3. Correct batch probability over time
    avg_correct_prob = df['correct_batch_prob'].rolling(window=window_size, center=True).mean()
    
    ax3.plot(df['sim_timestamp'], avg_correct_prob, color='blue', linewidth=2)
    ax3.set_xlabel('Simulation Time', fontsize=10)
    ax3.set_ylabel('Average Correct Batch Probability', fontsize=10)
    ax3.set_title('Temporal: Correct Batch Probability', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # 4. Accuracy over time
    accuracy = df['correct_batch_is_highest'].rolling(window=window_size, center=True).mean()
    
    ax4.plot(df['sim_timestamp'], accuracy, color='orange', linewidth=2)
    ax4.set_xlabel('Simulation Time', fontsize=10)
    ax4.set_ylabel('Accuracy (Highest Prob)', fontsize=10)
    ax4.set_title('Temporal: Adversary Accuracy', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('logs_analysis/5_temporal_changes.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_6_client_impact(df):
    """6. Impact of Number of Clients"""
    client_groups = df.groupby('n_clients').agg({
        'uniquely_identified': 'mean',
        'anonymity_set_size': 'mean',
        'correct_batch_prob': 'mean',
        'correct_batch_is_highest': 'mean'
    }).reset_index()
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # 1. Unique identification vs clients
    ax1.plot(client_groups['n_clients'], client_groups['uniquely_identified'], 
             marker='o', linewidth=2, markersize=8, color='red')
    ax1.set_xlabel('Number of Clients', fontsize=10)
    ax1.set_ylabel('Fraction Uniquely Identified', fontsize=10)
    ax1.set_title('Unique Identification vs Clients', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Anonymity set size vs clients
    ax2.plot(client_groups['n_clients'], client_groups['anonymity_set_size'], 
             marker='s', linewidth=2, markersize=8, color='green')
    ax2.set_xlabel('Number of Clients', fontsize=10)
    ax2.set_ylabel('Average Anonymity Set Size', fontsize=10)
    ax2.set_title('Anonymity Set Size vs Clients', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # 3. Correct batch probability vs clients
    ax3.plot(client_groups['n_clients'], client_groups['correct_batch_prob'], 
             marker='^', linewidth=2, markersize=8, color='blue')
    ax3.set_xlabel('Number of Clients', fontsize=10)
    ax3.set_ylabel('Average Correct Batch Probability', fontsize=10)
    ax3.set_title('Correct Batch Probability vs Clients', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # 4. Accuracy vs clients
    ax4.plot(client_groups['n_clients'], client_groups['correct_batch_is_highest'], 
             marker='d', linewidth=2, markersize=8, color='orange')
    ax4.set_xlabel('Number of Clients', fontsize=10)
    ax4.set_ylabel('Adversary Accuracy', fontsize=10)
    ax4.set_title('Adversary Accuracy vs Clients', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('logs_analysis/6_client_impact.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_7_batch_size_impact(df):
    """7. Impact of Batch Size"""
    batch_groups = df.groupby('batch_size').agg({
        'uniquely_identified': 'mean',
        'anonymity_set_size': 'mean',
        'correct_batch_prob': 'mean',
        'correct_batch_is_highest': 'mean'
    }).reset_index()
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # 1. Unique identification vs batch size
    ax1.plot(batch_groups['batch_size'], batch_groups['uniquely_identified'], 
             marker='o', linewidth=2, markersize=8, color='red')
    ax1.set_xlabel('Batch Size', fontsize=10)
    ax1.set_ylabel('Fraction Uniquely Identified', fontsize=10)
    ax1.set_title('Unique Identification vs Batch Size', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # 2. Anonymity set size vs batch size
    ax2.plot(batch_groups['batch_size'], batch_groups['anonymity_set_size'], 
             marker='s', linewidth=2, markersize=8, color='green')
    ax2.set_xlabel('Batch Size', fontsize=10)
    ax2.set_ylabel('Average Anonymity Set Size', fontsize=10)
    ax2.set_title('Anonymity Set Size vs Batch Size', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # 3. Correct batch probability vs batch size
    ax3.plot(batch_groups['batch_size'], batch_groups['correct_batch_prob'], 
             marker='^', linewidth=2, markersize=8, color='blue')
    ax3.set_xlabel('Batch Size', fontsize=10)
    ax3.set_ylabel('Average Correct Batch Probability', fontsize=10)
    ax3.set_title('Correct Batch Probability vs Batch Size', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # 4. Accuracy vs batch size
    ax4.plot(batch_groups['batch_size'], batch_groups['correct_batch_is_highest'], 
             marker='d', linewidth=2, markersize=8, color='orange')
    ax4.set_xlabel('Batch Size', fontsize=10)
    ax4.set_ylabel('Adversary Accuracy', fontsize=10)
    ax4.set_title('Adversary Accuracy vs Batch Size', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('logs_analysis/7_batch_size_impact.png', dpi=300, bbox_inches='tight')
    plt.show()

def generate_summary_statistics(df):
    """Generate summary statistics"""
    print("="*60)
    print("BATCH TRACKING ANALYSIS SUMMARY")
    print("="*60)
    
    print(f"Total batches analyzed: {len(df)}")
    print(f"Simulation time range: {df['sim_timestamp'].min():.2f} - {df['sim_timestamp'].max():.2f}")
    print(f"Number of clients: {df['n_clients'].unique()}")
    print(f"Batch sizes: {sorted(df['batch_size'].unique())}")
    
    print("\n" + "="*40)
    print("ANONYMITY METRICS")
    print("="*40)
    
    unique_rate = df['uniquely_identified'].mean()
    print(f"Unique identification rate: {unique_rate:.1%}")
    print(f"Average anonymity set size: {df['anonymity_set_size'].mean():.2f}")
    print(f"Median anonymity set size: {df['anonymity_set_size'].median():.1f}")
    
    print(f"\nCorrect batch probability:")
    print(f"  Mean: {df['correct_batch_prob'].mean():.3f}")
    print(f"  Median: {df['correct_batch_prob'].median():.3f}")
    print(f"  Min: {df['correct_batch_prob'].min():.3f}")
    print(f"  Max: {df['correct_batch_prob'].max():.3f}")
    
    accuracy = df['correct_batch_is_highest'].mean()
    print(f"\nAdversary accuracy (highest prob): {accuracy:.1%}")
    
    print("\n" + "="*40)
    print("PRIVACY ASSESSMENT")
    print("="*40)
    
    if unique_rate < 0.1:
        privacy_level = "EXCELLENT"
    elif unique_rate < 0.3:
        privacy_level = "GOOD"
    elif unique_rate < 0.5:
        privacy_level = "MODERATE"
    else:
        privacy_level = "POOR"
    
    print(f"Overall privacy level: {privacy_level}")
    print(f"Reason: {unique_rate:.1%} of batches are uniquely identified")

def main():
    """Main function to run all analyses"""
    # Create output directory
    Path('logs_analysis').mkdir(exist_ok=True)
    
    # Load data
    csv_file = 'batch_logs_19202858_6.csv'  # Update with your file path
    df = load_and_process_data(csv_file)
    
    print("Generating visualizations...")
    
    # Generate all plots
    plot_1_fraction_uniquely_identified(df)
    plot_2_anonymity_set_distribution(df)
    plot_3_correct_batch_probability(df)
    plot_4_accuracy_highest_probability(df)
    plot_5_temporal_changes(df)
    plot_6_client_impact(df)
    plot_7_batch_size_impact(df)
    
    # Generate summary
    generate_summary_statistics(df)
    
    print(f"\nAll visualizations saved to 'logs_analysis/' directory")
    print("Analysis complete!")

if __name__ == "__main__":
    main()