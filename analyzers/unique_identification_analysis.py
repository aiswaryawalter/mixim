import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from BatchLogAnalyzer import BatchLogAnalyzer
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class UniqueIdentificationAnalysis(BatchLogAnalyzer):
    
    def get_unique_identification_data(self, df):
        """Extract unique identification data from the DataFrame"""
        if df.empty:
            return [], []
        
        anonymity_sizes = df['anonymity_set_size'].tolist()
        batch_indices = df['out_batch_id'].tolist()
        
        # Identify uniquely identified batches (anonymity_set_size == 1)
        uniquely_identified = [1 if size == 1 else 0 for size in anonymity_sizes]
        
        total_batches = len(anonymity_sizes)
        unique_count = sum(uniquely_identified)
        unique_fraction = unique_count / total_batches if total_batches > 0 else 0
        
        print(f"Total batches: {total_batches}")
        print(f"Uniquely identified batches: {unique_count}")
        print(f"Fraction uniquely identified: {unique_fraction:.4f} ({unique_fraction*100:.1f}%)")
        
        return batch_indices, uniquely_identified
    
    def plot_unique_identification_pie(self, folder_data_dict):
        """Plot pie charts for unique identification fractions"""
        num_folders = len(folder_data_dict)
        if num_folders == 1:
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            axes = [ax]
        else:
            fig, axes = plt.subplots(1, num_folders, figsize=(8*num_folders, 8))
            if num_folders == 1:
                axes = [axes]
        
        colors = ['#ff9999', '#66b3ff']  # Red for unique, Blue for anonymous
        
        for i, (folder_name, (batch_indices, uniquely_identified)) in enumerate(folder_data_dict.items()):
            if not batch_indices:
                continue
            
            ax = axes[i] if num_folders > 1 else axes[0]
            
            unique_count = sum(uniquely_identified)
            anonymous_count = len(uniquely_identified) - unique_count
            total_count = len(uniquely_identified)
            
            # Data for pie chart
            sizes = [unique_count, anonymous_count]
            labels = [f'Uniquely Identified\n({unique_count}/{total_count})', 
                     f'Anonymous\n({anonymous_count}/{total_count})']
            percentages = [unique_count/total_count*100, anonymous_count/total_count*100]
            
            # Create pie chart
            wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                                             startangle=90, textprops={'fontsize': 12})
            
            # Enhance the appearance
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
                autotext.set_fontsize(14)
            
            ax.set_title(f'{folder_name}\nUnique Identification Analysis', 
                        fontsize=16, fontweight='bold', pad=20)
        
        plt.tight_layout()
        filename = self.generate_filename("unique_identification_pie")
        self.save_plot(filename)
        plt.show()
    
    def plot_unique_identification_comparison(self, folder_data_dict):
        """Plot bar chart comparing unique identification across folders"""
        self.setup_plot_style(figsize=(12, 8))
        
        folder_names = []
        unique_fractions = []
        anonymous_fractions = []
        total_counts = []
        
        for folder_name, (batch_indices, uniquely_identified) in folder_data_dict.items():
            if not batch_indices:
                continue
            
            unique_count = sum(uniquely_identified)
            total_count = len(uniquely_identified)
            unique_fraction = unique_count / total_count
            anonymous_fraction = 1 - unique_fraction
            
            folder_names.append(folder_name)
            unique_fractions.append(unique_fraction)
            anonymous_fractions.append(anonymous_fraction)
            total_counts.append(total_count)
        
        x_positions = np.arange(len(folder_names))
        width = 0.6
        
        # Create stacked bar chart
        bars1 = plt.bar(x_positions, unique_fractions, width, 
                       label='Uniquely Identified', color='#ff9999', alpha=0.8)
        bars2 = plt.bar(x_positions, anonymous_fractions, width, 
                       bottom=unique_fractions, label='Anonymous', color='#66b3ff', alpha=0.8)
        
        # Add value labels on bars
        for i, (unique_frac, anon_frac, total) in enumerate(zip(unique_fractions, anonymous_fractions, total_counts)):
            # Label for unique identification portion
            if unique_frac > 0.05:  # Only show if segment is large enough
                plt.text(i, unique_frac/2, f'{unique_frac:.2%}\n({int(unique_frac*total)})', 
                        ha='center', va='center', fontweight='bold', color='white')
            
            # Label for anonymous portion
            if anon_frac > 0.05:  # Only show if segment is large enough
                plt.text(i, unique_frac + anon_frac/2, f'{anon_frac:.2%}\n({int(anon_frac*total)})', 
                        ha='center', va='center', fontweight='bold', color='white')
        
        plt.xlabel('Dataset', fontsize=14)
        plt.ylabel('Fraction of Batches', fontsize=14)
        plt.title('Comparison of Unique Identification Rates\nAcross Datasets', fontsize=16, fontweight='bold')
        plt.xticks(x_positions, folder_names)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3, axis='y')
        plt.ylim(0, 1)
        
        # Add total batch counts as text
        for i, total in enumerate(total_counts):
            plt.text(i, 1.02, f'n={total}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        filename = self.generate_filename("unique_identification_comparison")
        self.save_plot(filename)
        plt.show()
    
    def plot_anonymity_distribution_breakdown(self, folder_data_dict):
        """Plot detailed breakdown of anonymity set sizes"""
        self.setup_plot_style(figsize=(14, 8))
        
        colors = ['#ff4444', '#ff9999', '#66b3ff', '#4444ff', '#9999ff']
        
        for folder_name, (batch_indices, _) in folder_data_dict.items():
            if not batch_indices:
                continue
            
            # Get original anonymity sizes
            df = self.load_log_files(folder_name)
            window_data = self.extract_window_data(df, None)
            anonymity_sizes = window_data['anonymity_set_size'].tolist()
            
            # Count occurrences of each anonymity set size
            size_counts = pd.Series(anonymity_sizes).value_counts().sort_index()
            total_batches = len(anonymity_sizes)
            
            # Create detailed breakdown
            sizes = size_counts.index.tolist()
            counts = size_counts.values.tolist()
            percentages = [count/total_batches*100 for count in counts]
            
            # Plot as horizontal bar chart
            plt.figure(figsize=(10, 6))
            bars = plt.barh(range(len(sizes)), percentages, color=colors[:len(sizes)], alpha=0.8)
            
            # Add value labels
            for i, (size, count, pct) in enumerate(zip(sizes, counts, percentages)):
                plt.text(pct + 1, i, f'{count} batches ({pct:.1f}%)', 
                        va='center', fontsize=11)
                
                # Highlight unique identification (size == 1)
                if size == 1:
                    bars[i].set_color('#ff0000')
                    bars[i].set_alpha(1.0)
                    plt.text(pct/2, i, 'UNIQUELY\nIDENTIFIED', 
                            ha='center', va='center', fontweight='bold', 
                            color='white', fontsize=10)
            
            plt.xlabel('Percentage of Batches', fontsize=14)
            plt.ylabel('Anonymity Set Size', fontsize=14)
            plt.title(f'Detailed Anonymity Set Size Distribution\n{folder_name} (Total: {total_batches} batches)', 
                     fontsize=16, fontweight='bold')
            plt.yticks(range(len(sizes)), [f'Size {size}' for size in sizes])
            plt.grid(True, alpha=0.3, axis='x')
            
            # Add privacy assessment
            unique_fraction = percentages[0] if sizes[0] == 1 else 0
            if unique_fraction > 50:
                assessment = "CRITICAL: Majority uniquely identified"
                color = 'red'
            elif unique_fraction > 25:
                assessment = "HIGH RISK: Many batches uniquely identified"
                color = 'orange'
            elif unique_fraction > 10:
                assessment = "MODERATE RISK: Some unique identification"
                color = 'yellow'
            else:
                assessment = "LOW RISK: Good anonymity protection"
                color = 'green'
            
            plt.text(0.98, 0.98, f'Privacy Assessment:\n{assessment}', 
                    transform=plt.gca().transAxes, fontsize=12, 
                    ha='right', va='top', fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
            
            plt.tight_layout()
            filename = self.generate_filename(f"anonymity_breakdown_{folder_name}")
            self.save_plot(filename)
            plt.show()
    
    def calculate_privacy_metrics(self, folder_data_dict):
        """Calculate detailed privacy metrics"""
        print("\n" + "="*60)
        print("UNIQUE IDENTIFICATION PRIVACY METRICS")
        print("="*60)
        
        summary_data = []
        
        for folder_name, (batch_indices, uniquely_identified) in folder_data_dict.items():
            if not batch_indices:
                continue
            
            total_batches = len(uniquely_identified)
            unique_count = sum(uniquely_identified)
            anonymous_count = total_batches - unique_count
            unique_fraction = unique_count / total_batches
            
            # Get original anonymity sizes for more detailed analysis
            df = self.load_log_files(folder_name)
            window_data = self.extract_window_data(df, None)
            anonymity_sizes = window_data['anonymity_set_size'].tolist()
            
            # Calculate additional metrics
            avg_anonymity = np.mean(anonymity_sizes)
            median_anonymity = np.median(anonymity_sizes)
            max_anonymity = max(anonymity_sizes)
            std_anonymity = np.std(anonymity_sizes)
            
            # Privacy score (0 = worst, 1 = best)
            privacy_score = 1 - unique_fraction
            
            # Risk assessment
            if unique_fraction > 0.5:
                risk_level = "CRITICAL"
            elif unique_fraction > 0.25:
                risk_level = "HIGH"
            elif unique_fraction > 0.1:
                risk_level = "MODERATE"
            else:
                risk_level = "LOW"
            
            print(f"\n{folder_name}:")
            print(f"  Total batches: {total_batches}")
            print(f"  Uniquely identified: {unique_count} ({unique_fraction:.2%})")
            print(f"  Anonymous: {anonymous_count} ({1-unique_fraction:.2%})")
            print(f"  Average anonymity set size: {avg_anonymity:.2f}")
            print(f"  Median anonymity set size: {median_anonymity}")
            print(f"  Maximum anonymity set size: {max_anonymity}")
            print(f"  Anonymity std deviation: {std_anonymity:.2f}")
            print(f"  Privacy score: {privacy_score:.3f}")
            print(f"  Risk level: {risk_level}")
            
            # Store for summary table
            summary_data.append({
                'Dataset': folder_name,
                'Total_Batches': total_batches,
                'Uniquely_Identified': unique_count,
                'Unique_Fraction': unique_fraction,
                'Avg_Anonymity': avg_anonymity,
                'Privacy_Score': privacy_score,
                'Risk_Level': risk_level
            })
        
        # Create summary DataFrame
        summary_df = pd.DataFrame(summary_data)
        print("\n" + "="*60)
        print("SUMMARY TABLE")
        print("="*60)
        print(summary_df.round(3).to_string(index=False))
        
        return summary_df
    
    def analyze_unique_identification(self, folders=['12_logs', '24_logs'], window_index=None):
        """Main analysis method for unique identification"""
        print("="*60)
        print("UNIQUE IDENTIFICATION ANALYSIS")
        print("Measuring fraction of batches with anonymity_set_size == 1")
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
            
            batch_indices, uniquely_identified = self.get_unique_identification_data(window_data)
            folder_data[folder] = (batch_indices, uniquely_identified)
        
        if folder_data:
            summary_df = self.calculate_privacy_metrics(folder_data)
            self.plot_unique_identification_pie(folder_data)
            self.plot_unique_identification_comparison(folder_data)
            self.plot_anonymity_distribution_breakdown(folder_data)
        
        return folder_data, summary_df

def main():
    analyzer = UniqueIdentificationAnalysis()
    analyzer.analyze_unique_identification()

if __name__ == "__main__":
    main()