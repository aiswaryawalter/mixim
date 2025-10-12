import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from BatchLogAnalyzer import BatchLogAnalyzer
import matplotlib.pyplot as plt
import numpy as np

class BatchProbabilityAnalysis(BatchLogAnalyzer):
    
    def get_correct_batch_probabilities(self, df):
        """Extract correct batch probabilities"""
        if df.empty:
            return [], []
        
        batch_indices = df['out_batch_id'].tolist()
        correct_probs = df['correct_batch_prob'].tolist()
        
        print(f"Probabilities range: {min(correct_probs):.4f} to {max(correct_probs):.4f}")
        print(f"Mean probability: {np.mean(correct_probs):.4f}")
        
        return batch_indices, correct_probs
    
    def plot_scatter(self, folder_data_dict):
        """Plot scatter plot of probabilities vs batch index"""
        self.setup_plot_style()
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        markers = ['o', 's', '^', 'D', 'v']
        
        for i, (folder_name, (batch_indices, correct_probs)) in enumerate(folder_data_dict.items()):
            if not batch_indices:
                continue
            
            plt.scatter(batch_indices, correct_probs, 
                       alpha=0.7, label=f'{folder_name} (n={len(batch_indices)})', 
                       color=colors[i % len(colors)], s=50, 
                       marker=markers[i % len(markers)], 
                       edgecolors='black', linewidth=0.5)
        
        plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='Random Guess')
        plt.xlabel('Outgoing Batch Index', fontsize=14)
        plt.ylabel('Probability Assigned to Correct Batch', fontsize=14)
        plt.title('Correct Batch Probability vs Batch Index', fontsize=16, fontweight='bold')
        plt.legend(loc='best', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = self.generate_filename("batch_probability_scatter")
        self.save_plot(filename)
        plt.show()
    
    def plot_distribution(self, folder_data_dict):
        """Plot distribution of probabilities"""
        self.setup_plot_style(figsize=(12, 8))
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'orange', 'plum']
        
        for i, (folder_name, (batch_indices, correct_probs)) in enumerate(folder_data_dict.items()):
            if not correct_probs:
                continue
            
            plt.hist(correct_probs, bins=20, alpha=0.7, 
                    label=f'{folder_name} (n={len(correct_probs)})', 
                    color=colors[i % len(colors)], edgecolor='black', density=True)
        
        plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.8, label='Random Guess')
        plt.xlabel('Probability Assigned to Correct Batch', fontsize=14)
        plt.ylabel('Density', fontsize=14)
        plt.title('Distribution of Correct Batch Probabilities', fontsize=16, fontweight='bold')
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filename = self.generate_filename("batch_probability_distribution")
        self.save_plot(filename)
        plt.show()
    
    def analyze_batch_probabilities(self, folders=['12_logs', '24_logs'], window_index=None):
        """Main analysis method"""
        print("="*60)
        print("BATCH PROBABILITY ANALYSIS")
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
            
            batch_indices, correct_probs = self.get_correct_batch_probabilities(window_data)
            folder_data[folder] = (batch_indices, correct_probs)
        
        if folder_data:
            self.plot_scatter(folder_data)
            self.plot_distribution(folder_data)
        
        return folder_data

def main():
    analyzer = BatchProbabilityAnalysis()
    analyzer.analyze_batch_probabilities()

if __name__ == "__main__":
    main()