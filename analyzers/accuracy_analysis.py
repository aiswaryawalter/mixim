import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from BatchLogAnalyzer import BatchLogAnalyzer
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class AccuracyAnalysis(BatchLogAnalyzer):
    
    def get_accuracy_data(self, df):
        """Extract accuracy data from the DataFrame"""
        if df.empty:
            return [], [], []
        
        batch_indices = df['out_batch_id'].tolist()
        correct_batch_is_highest = df['correct_batch_is_highest'].tolist()
        window_indices = df['window_index'].tolist()
        
        # Convert boolean values to integers (True=1, False=0)
        accuracy_values = [1 if val else 0 for val in correct_batch_is_highest]
        
        print(f"Total batches analyzed: {len(batch_indices)}")
        print(f"Correct batch has highest probability: {sum(accuracy_values)} times")
        print(f"Overall accuracy rate: {np.mean(accuracy_values):.4f} ({np.mean(accuracy_values)*100:.1f}%)")
        
        return batch_indices, accuracy_values, window_indices
    
    def calculate_accuracy_over_time(self, batch_indices, accuracy_values):
        """Calculate cumulative accuracy over time"""
        cumulative_correct = np.cumsum(accuracy_values)
        cumulative_total = np.arange(1, len(accuracy_values) + 1)
        cumulative_accuracy = cumulative_correct / cumulative_total
        
        return batch_indices, cumulative_accuracy
    
    def calculate_window_accuracy(self, df):
        """Calculate accuracy per window"""
        if df.empty:
            return {}, {}
        
        window_accuracy = {}
        window_counts = {}
        
        for window in df['window_index'].unique():
            window_data = df[df['window_index'] == window]
            correct_count = window_data['correct_batch_is_highest'].sum()
            total_count = len(window_data)
            accuracy = correct_count / total_count if total_count > 0 else 0
            
            window_accuracy[window] = accuracy
            window_counts[window] = total_count
            
        return window_accuracy, window_counts
    
    def plot_accuracy_over_time(self, folder_data_dict):
        """Plot cumulative accuracy over time (batch progression)"""
        self.setup_plot_style(figsize=(15, 8))
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        for i, (folder_name, (batch_indices, accuracy_values, _)) in enumerate(folder_data_dict.items()):
            if not batch_indices:
                continue
            
            # Calculate cumulative accuracy
            batch_nums, cum_accuracy = self.calculate_accuracy_over_time(batch_indices, accuracy_values)
            
            plt.plot(batch_nums, cum_accuracy, 
                    label=f'{folder_name} (final: {cum_accuracy[-1]:.3f})', 
                    color=colors[i % len(colors)], linewidth=2, marker='o', markersize=4)
        
        plt.xlabel('Batch Index', fontsize=14)
        plt.ylabel('Cumulative Accuracy Rate', fontsize=14)
        plt.title('Adversary Success Rate Over Time\n(Cumulative Accuracy of Highest Probability Prediction)', 
                 fontsize=16, fontweight='bold')
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1.1)
        
        # Add horizontal lines for reference
        plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% accuracy')
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='Perfect accuracy')
        
        plt.tight_layout()
        filename = self.generate_filename("accuracy_over_time")
        self.save_plot(filename)
        plt.show()
    
    def plot_window_accuracy(self, folder_data_dict):
        """Plot accuracy per window as bar chart"""
        self.setup_plot_style(figsize=(15, 8))
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'orange', 'plum']
        
        all_windows = set()
        folder_window_data = {}
        
        # Collect all window data
        for folder_name, (_, _, _) in folder_data_dict.items():
            # Need to recalculate per window from original data
            df = self.load_log_files(folder_name)
            window_data = self.extract_window_data(df, None)  # Use largest window
            window_accuracy, window_counts = self.calculate_window_accuracy(df)
            
            folder_window_data[folder_name] = window_accuracy
            all_windows.update(window_accuracy.keys())
        
        all_windows = sorted(all_windows)
        width = 0.35
        x_positions = np.arange(len(all_windows))
        
        # Plot bars for each folder
        for i, (folder_name, window_accuracy) in enumerate(folder_window_data.items()):
            accuracies = [window_accuracy.get(window, 0) for window in all_windows]
            
            plt.bar(x_positions + i * width, accuracies, width,
                   label=folder_name, color=colors[i % len(colors)], 
                   alpha=0.8, edgecolor='black', linewidth=0.5)
        
        plt.xlabel('Window Index', fontsize=14)
        plt.ylabel('Accuracy Rate', fontsize=14)
        plt.title('Adversary Success Rate by Window\n(Fraction of Batches Where Correct Batch Has Highest Probability)', 
                 fontsize=16, fontweight='bold')
        plt.xticks(x_positions + width/2, [str(w) for w in all_windows])
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3, axis='y')
        plt.ylim(0, 1.1)
        
        # Add reference line
        plt.axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label='50% accuracy')
        
        plt.tight_layout()
        filename = self.generate_filename("accuracy_by_window")
        self.save_plot(filename)
        plt.show()
    
    def plot_accuracy_distribution(self, folder_data_dict):
        """Plot distribution of accuracy values (0 or 1)"""
        self.setup_plot_style(figsize=(10, 6))
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'orange', 'plum']
        
        x_labels = ['Incorrect\n(Highest ≠ Correct)', 'Correct\n(Highest = Correct)']
        x_positions = [0, 1]
        width = 0.35
        
        for i, (folder_name, (_, accuracy_values, _)) in enumerate(folder_data_dict.items()):
            if not accuracy_values:
                continue
            
            # Count 0s and 1s
            incorrect_count = accuracy_values.count(0)
            correct_count = accuracy_values.count(1)
            counts = [incorrect_count, correct_count]
            
            plt.bar([pos + i * width for pos in x_positions], counts, width,
                   label=f'{folder_name} (acc: {np.mean(accuracy_values):.3f})', 
                   color=colors[i % len(colors)], alpha=0.8, 
                   edgecolor='black', linewidth=0.5)
        
        plt.xlabel('Prediction Outcome', fontsize=14)
        plt.ylabel('Number of Batches', fontsize=14)
        plt.title('Distribution of Adversary Prediction Outcomes', fontsize=16, fontweight='bold')
        plt.xticks([pos + width/2 for pos in x_positions], x_labels)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        filename = self.generate_filename("accuracy_distribution")
        self.save_plot(filename)
        plt.show()
    
    def print_detailed_accuracy_stats(self, folder_data_dict):
        """Print detailed accuracy statistics"""
        print("\n" + "="*60)
        print("DETAILED ACCURACY STATISTICS")
        print("="*60)
        
        for folder_name, (batch_indices, accuracy_values, window_indices) in folder_data_dict.items():
            if not accuracy_values:
                continue
            
            total_batches = len(accuracy_values)
            correct_predictions = sum(accuracy_values)
            accuracy_rate = np.mean(accuracy_values)
            
            print(f"\n{folder_name}:")
            print(f"  Total batches: {total_batches}")
            print(f"  Correct predictions: {correct_predictions}")
            print(f"  Incorrect predictions: {total_batches - correct_predictions}")
            print(f"  Success rate: {accuracy_rate:.4f} ({accuracy_rate*100:.1f}%)")
            
            # Security interpretation
            if accuracy_rate > 0.5:
                risk_level = "HIGH RISK"
                interpretation = "Adversary performs better than random guessing"
            elif accuracy_rate == 0.5:
                risk_level = "MODERATE RISK"
                interpretation = "Adversary performs at random level"
            else:
                risk_level = "LOW RISK"
                interpretation = "Adversary performs worse than random"
            
            print(f"  Security assessment: {risk_level}")
            print(f"  Interpretation: {interpretation}")
            
            # Calculate confidence interval (for large samples)
            if total_batches >= 30:
                std_error = np.sqrt(accuracy_rate * (1 - accuracy_rate) / total_batches)
                ci_lower = accuracy_rate - 1.96 * std_error
                ci_upper = accuracy_rate + 1.96 * std_error
                print(f"  95% Confidence Interval: [{ci_lower:.3f}, {ci_upper:.3f}]")
    
    def analyze_adversary_accuracy(self, folders=['12_logs', '24_logs'], window_index=None):
        """Main analysis method for adversary accuracy"""
        print("="*60)
        print("ADVERSARY ACCURACY ANALYSIS")
        print("Measuring success rate when always picking highest probability batch")
        print("="*60)
        
        folder_data = {}
        for folder in folders:
            print(f"\nProcessing folder: {folder}")
            df = self.load_log_files(folder)
            if df.empty:
                continue
            
            window_data = self.extract_window_data(df, window_index)
            if window_data.empty:
                continue
            
            batch_indices, accuracy_values, window_indices = self.get_accuracy_data(window_data)
            folder_data[folder] = (batch_indices, accuracy_values, window_indices)
        
        if folder_data:
            self.print_detailed_accuracy_stats(folder_data)
            self.plot_accuracy_over_time(folder_data)
            self.plot_window_accuracy(folder_data)
            self.plot_accuracy_distribution(folder_data)
        
        return folder_data
    
    def compare_accuracy_across_windows(self, folders=['12_logs', '24_logs']):
        """Compare accuracy across all windows"""
        print("\n" + "="*60)
        print("ACCURACY COMPARISON ACROSS ALL WINDOWS")
        print("="*60)
        
        for folder in folders:
            print(f"\nAnalyzing all windows in {folder}:")
            df = self.load_log_files(folder)
            if df.empty:
                continue
            
            window_accuracy, window_counts = self.calculate_window_accuracy(df)
            
            for window in sorted(window_accuracy.keys()):
                accuracy = window_accuracy[window]
                count = window_counts[window]
                print(f"  Window {window}: {accuracy:.3f} accuracy ({count} batches)")

def main():
    analyzer = AccuracyAnalysis()
    
    # Analyze largest window
    analyzer.analyze_adversary_accuracy()
    
    # Compare across all windows
    analyzer.compare_accuracy_across_windows()

if __name__ == "__main__":
    main()