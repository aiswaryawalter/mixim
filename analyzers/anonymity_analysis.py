import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from BatchLogAnalyzer import BatchLogAnalyzer
import matplotlib.pyplot as plt
import numpy as np

class AnonymityAnalysis(BatchLogAnalyzer):
    
    def get_anonymity_set_sizes(self, df):
        """Extract anonymity set sizes from the DataFrame"""
        if df.empty:
            return []
        
        anonymity_sizes = df['anonymity_set_size'].tolist()
        print(f"Anonymity set sizes range: {min(anonymity_sizes)} to {max(anonymity_sizes)}")
        print(f"Mean anonymity set size: {np.mean(anonymity_sizes):.2f}")
        return anonymity_sizes
    
    def plot_anonymity_distribution(self, folder_data_dict):
        """Plot anonymity set size distribution"""
        self.setup_plot_style()
        
        colors = ['skyblue', 'lightcoral', 'lightgreen', 'orange', 'plum']
        all_sizes = []
        for sizes in folder_data_dict.values():
            all_sizes.extend(sizes)
        
        if not all_sizes:
            print("No data to plot")
            return
        
        bin_range = range(min(all_sizes), max(all_sizes) + 2)
        width = 0.35
        x_positions = np.arange(len(bin_range) - 1)
        
        for i, (folder_name, anonymity_sizes) in enumerate(folder_data_dict.items()):
            if not anonymity_sizes:
                continue
            
            counts, _ = np.histogram(anonymity_sizes, bins=bin_range)
            plt.bar(x_positions + i * width, counts, width, 
                   alpha=0.8, label=f'{folder_name} (n={len(anonymity_sizes)})', 
                   color=colors[i % len(colors)], edgecolor='black', linewidth=0.5)
        
        plt.xlabel('Anonymity Set Size', fontsize=14)
        plt.ylabel('Count', fontsize=14)
        plt.title('Anonymity Set Size Distribution Across Batches', fontsize=16, fontweight='bold')
        plt.xticks(x_positions + width/2, [str(i) for i in bin_range[:-1]])
        plt.legend(loc='upper right', fontsize=12)
        plt.grid(True, alpha=0.3, axis='y')
        
        # Add statistics
        stats_text = ""
        for folder_name, anonymity_sizes in folder_data_dict.items():
            if anonymity_sizes:
                mean_size = np.mean(anonymity_sizes)
                std_size = np.std(anonymity_sizes)
                median_size = np.median(anonymity_sizes)
                stats_text += f"{folder_name}: μ={mean_size:.1f}, σ={std_size:.1f}, median={median_size}\n"
        
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                fontsize=11, verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        
        filename = self.generate_filename("anonymity_distribution")
        self.save_plot(filename)
        plt.show()
    
    def analyze_anonymity_sets(self, folders=['12_logs', '24_logs'], window_index=None):
        """Main method to analyze anonymity set sizes"""
        print("="*60)
        print("ANONYMITY SET SIZE ANALYSIS")
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
            
            anonymity_sizes = self.get_anonymity_set_sizes(window_data)
            folder_data[folder] = anonymity_sizes
        
        if folder_data:
            self.plot_anonymity_distribution(folder_data)
        
        return folder_data

def main():
    analyzer = AnonymityAnalysis()
    analyzer.analyze_anonymity_sets()

if __name__ == "__main__":
    main()