import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from BatchLogAnalyzer import BatchLogAnalyzer
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

class TemporalAnalysis(BatchLogAnalyzer):
    
    def extract_temporal_data(self, df):
        """Extract temporal metrics from the DataFrame"""
        if df.empty:
            return {}
        
        temporal_data = {
            'batch_indices': df['out_batch_id'].tolist(),
            'anonymity_sizes': df['anonymity_set_size'].tolist(),
            'correct_probs': df['correct_batch_prob'].tolist(),
            'correct_is_highest': df['correct_batch_is_highest'].tolist(),
            'window_indices': df['window_index'].tolist(),
            'sim_timestamps': df['sim_timestamp'].tolist() if 'sim_timestamp' in df.columns else None
        }
        
        print(f"Extracted temporal data for {len(temporal_data['batch_indices'])} batches")
        return temporal_data
    
    def calculate_rolling_metrics(self, temporal_data, window_size=5):
        """Calculate rolling averages for metrics"""
        batch_indices = temporal_data['batch_indices']
        anonymity_sizes = temporal_data['anonymity_sizes']
        correct_probs = temporal_data['correct_probs']
        correct_is_highest = [1 if x else 0 for x in temporal_data['correct_is_highest']]
        
        # Convert to pandas Series for rolling calculations
        anonymity_series = pd.Series(anonymity_sizes, index=batch_indices)
        prob_series = pd.Series(correct_probs, index=batch_indices)
        accuracy_series = pd.Series(correct_is_highest, index=batch_indices)
        
        # Calculate rolling means
        rolling_anonymity = anonymity_series.rolling(window=window_size, min_periods=1).mean()
        rolling_prob = prob_series.rolling(window=window_size, min_periods=1).mean()
        rolling_accuracy = accuracy_series.rolling(window=window_size, min_periods=1).mean()
        
        # Calculate fraction uniquely identified (anonymity_size = 1)
        uniquely_identified = pd.Series([1 if size == 1 else 0 for size in anonymity_sizes], index=batch_indices)
        rolling_unique = uniquely_identified.rolling(window=window_size, min_periods=1).mean()
        
        return {
            'batch_indices': batch_indices,
            'rolling_anonymity': rolling_anonymity.values,
            'rolling_prob': rolling_prob.values,
            'rolling_accuracy': rolling_accuracy.values,
            'rolling_unique_fraction': rolling_unique.values,
            'raw_anonymity': anonymity_sizes,
            'raw_prob': correct_probs,
            'raw_accuracy': correct_is_highest
        }
    
    def plot_temporal_evolution(self, folder_data_dict, window_size=5):
        """Plot how metrics evolve over time"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        # Plot 1: Anonymity Set Size Over Time
        ax1 = axes[0, 0]
        for i, (folder_name, temporal_data) in enumerate(folder_data_dict.items()):
            if not temporal_data:
                continue
            
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            
            # Plot raw data as scatter
            ax1.scatter(rolling_data['batch_indices'], rolling_data['raw_anonymity'], 
                       alpha=0.3, s=20, color=colors[i % len(colors)])
            
            # Plot rolling average as line
            ax1.plot(rolling_data['batch_indices'], rolling_data['rolling_anonymity'], 
                    label=f'{folder_name} (rolling avg)', color=colors[i % len(colors)], linewidth=2)
        
        ax1.set_xlabel('Batch Index')
        ax1.set_ylabel('Anonymity Set Size')
        ax1.set_title('Anonymity Set Size Evolution Over Time')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Fraction Uniquely Identified Over Time
        ax2 = axes[0, 1]
        for i, (folder_name, temporal_data) in enumerate(folder_data_dict.items()):
            if not temporal_data:
                continue
            
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            
            ax2.plot(rolling_data['batch_indices'], rolling_data['rolling_unique_fraction'], 
                    label=f'{folder_name}', color=colors[i % len(colors)], linewidth=2, marker='o', markersize=4)
        
        ax2.set_xlabel('Batch Index')
        ax2.set_ylabel('Fraction Uniquely Identified')
        ax2.set_title('Fraction of Batches with Anonymity Set Size = 1')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)
        
        # Plot 3: Correct Batch Probability Over Time
        ax3 = axes[1, 0]
        for i, (folder_name, temporal_data) in enumerate(folder_data_dict.items()):
            if not temporal_data:
                continue
            
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            
            # Plot raw probabilities as scatter
            ax3.scatter(rolling_data['batch_indices'], rolling_data['raw_prob'], 
                       alpha=0.3, s=20, color=colors[i % len(colors)])
            
            # Plot rolling average
            ax3.plot(rolling_data['batch_indices'], rolling_data['rolling_prob'], 
                    label=f'{folder_name} (rolling avg)', color=colors[i % len(colors)], linewidth=2)
        
        ax3.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='Random Guess')
        ax3.set_xlabel('Batch Index')
        ax3.set_ylabel('Correct Batch Probability')
        ax3.set_title('Correct Batch Probability Evolution')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Adversary Success Rate Over Time
        ax4 = axes[1, 1]
        for i, (folder_name, temporal_data) in enumerate(folder_data_dict.items()):
            if not temporal_data:
                continue
            
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            
            ax4.plot(rolling_data['batch_indices'], rolling_data['rolling_accuracy'], 
                    label=f'{folder_name}', color=colors[i % len(colors)], linewidth=2, marker='s', markersize=4)
        
        ax4.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='Random Guess')
        ax4.set_xlabel('Batch Index')
        ax4.set_ylabel('Adversary Success Rate')
        ax4.set_title('Adversary Success Rate Evolution')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        ax4.set_ylim(0, 1)
        
        plt.tight_layout()
        filename = self.generate_filename("temporal_evolution")
        self.save_plot(filename)
        plt.show()
    
    def plot_trend_analysis(self, folder_data_dict, window_size=5):
        """Plot trend analysis with statistical significance"""
        self.setup_plot_style(figsize=(16, 8))
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        # Plot 1: Anonymity trends
        ax1 = axes[0]
        for i, (folder_name, temporal_data) in enumerate(folder_data_dict.items()):
            if not temporal_data:
                continue
            
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            x = np.array(rolling_data['batch_indices'])
            y = np.array(rolling_data['rolling_anonymity'])
            
            # Plot the data
            ax1.plot(x, y, label=f'{folder_name}', color=colors[i % len(colors)], linewidth=2, marker='o', markersize=4)
            
            # Calculate and plot trend line
            if len(x) > 2:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                trend_line = slope * x + intercept
                
                # Plot trend line
                linestyle = '-' if p_value < 0.05 else '--'
                alpha = 0.8 if p_value < 0.05 else 0.5
                ax1.plot(x, trend_line, color=colors[i % len(colors)], linestyle=linestyle, alpha=alpha, linewidth=1)
                
                # Add trend info to label
                trend_direction = "↗" if slope > 0 else "↘" if slope < 0 else "→"
                significance = "*" if p_value < 0.05 else ""
                ax1.text(0.02, 0.98 - i*0.1, f'{folder_name}: {trend_direction} slope={slope:.4f}{significance}', 
                        transform=ax1.transAxes, fontsize=10, color=colors[i % len(colors)])
        
        ax1.set_xlabel('Batch Index')
        ax1.set_ylabel('Anonymity Set Size (Rolling Average)')
        ax1.set_title('Anonymity Size Trends Over Time\n(* = statistically significant)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Privacy degradation trends
        ax2 = axes[1]
        for i, (folder_name, temporal_data) in enumerate(folder_data_dict.items()):
            if not temporal_data:
                continue
            
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            x = np.array(rolling_data['batch_indices'])
            y = np.array(rolling_data['rolling_accuracy'])  # Higher = worse privacy
            
            # Plot the data
            ax2.plot(x, y, label=f'{folder_name}', color=colors[i % len(colors)], linewidth=2, marker='s', markersize=4)
            
            # Calculate and plot trend line
            if len(x) > 2:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                trend_line = slope * x + intercept
                
                # Plot trend line
                linestyle = '-' if p_value < 0.05 else '--'
                alpha = 0.8 if p_value < 0.05 else 0.5
                ax2.plot(x, trend_line, color=colors[i % len(colors)], linestyle=linestyle, alpha=alpha, linewidth=1)
                
                # Add trend info to label
                if slope > 0:
                    trend_direction = "↗ (degrading)"
                elif slope < 0:
                    trend_direction = "↘ (improving)"
                else:
                    trend_direction = "→ (stable)"
                significance = "*" if p_value < 0.05 else ""
                ax2.text(0.02, 0.98 - i*0.1, f'{folder_name}: {trend_direction} slope={slope:.4f}{significance}', 
                        transform=ax2.transAxes, fontsize=10, color=colors[i % len(colors)])
        
        ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='Random Guess')
        ax2.set_xlabel('Batch Index')
        ax2.set_ylabel('Adversary Success Rate (Rolling Average)')
        ax2.set_title('Privacy Degradation Trends Over Time\n(* = statistically significant)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)
        
        plt.tight_layout()
        filename = self.generate_filename("temporal_trends")
        self.save_plot(filename)
        plt.show()
    
    def calculate_temporal_statistics(self, folder_data_dict, window_size=5):
        """Calculate detailed temporal statistics"""
        print("\n" + "="*60)
        print("TEMPORAL EVOLUTION STATISTICS")
        print("="*60)
        
        for folder_name, temporal_data in folder_data_dict.items():
            if not temporal_data:
                continue
            
            print(f"\n{folder_name}:")
            rolling_data = self.calculate_rolling_metrics(temporal_data, window_size)
            
            # Anonymity trends
            x = np.array(rolling_data['batch_indices'])
            anonymity_y = np.array(rolling_data['rolling_anonymity'])
            accuracy_y = np.array(rolling_data['rolling_accuracy'])
            unique_y = np.array(rolling_data['rolling_unique_fraction'])
            
            if len(x) > 2:
                # Anonymity size trend
                anon_slope, _, anon_r, anon_p, _ = stats.linregress(x, anonymity_y)
                print(f"  Anonymity Size Trend:")
                print(f"    Slope: {anon_slope:.6f} (per batch)")
                print(f"    R²: {anon_r**2:.4f}")
                print(f"    P-value: {anon_p:.4f}")
                print(f"    Interpretation: {'Improving' if anon_slope > 0 else 'Degrading' if anon_slope < 0 else 'Stable'}")
                
                # Adversary success trend
                acc_slope, _, acc_r, acc_p, _ = stats.linregress(x, accuracy_y)
                print(f"  Adversary Success Trend:")
                print(f"    Slope: {acc_slope:.6f} (per batch)")
                print(f"    R²: {acc_r**2:.4f}")
                print(f"    P-value: {acc_p:.4f}")
                print(f"    Interpretation: {'Privacy degrading' if acc_slope > 0 else 'Privacy improving' if acc_slope < 0 else 'Stable'}")
                
                # Unique identification trend
                unique_slope, _, unique_r, unique_p, _ = stats.linregress(x, unique_y)
                print(f"  Unique Identification Trend:")
                print(f"    Slope: {unique_slope:.6f} (per batch)")
                print(f"    R²: {unique_r**2:.4f}")
                print(f"    P-value: {unique_p:.4f}")
                print(f"    Interpretation: {'More uniquely identified' if unique_slope > 0 else 'Less uniquely identified' if unique_slope < 0 else 'Stable'}")
                
                # Overall assessment
                print(f"  Overall Privacy Assessment:")
                privacy_indicators = []
                if anon_p < 0.05:
                    privacy_indicators.append("anonymity" + (" improving" if anon_slope > 0 else " degrading"))
                if acc_p < 0.05:
                    privacy_indicators.append("adversary success" + (" increasing" if acc_slope > 0 else " decreasing"))
                if unique_p < 0.05:
                    privacy_indicators.append("unique identification" + (" increasing" if unique_slope > 0 else " decreasing"))
                
                if privacy_indicators:
                    print(f"    Significant trends: {', '.join(privacy_indicators)}")
                else:
                    print(f"    No statistically significant trends detected")
    
    def analyze_temporal_changes(self, folders=['12_logs', '24_logs'], window_index=None, rolling_window=5):
        """Main analysis method for temporal changes"""
        print("="*60)
        print("TEMPORAL ANALYSIS - PRIVACY METRICS EVOLUTION")
        print("="*60)
        
        folder_data = {}
        for folder in folders:
            print(f"\nProcessing folder: {folder}")
            df = self.load_log_files(folder)
            if df.empty:
                continue
            
            if window_index is None:
                # Use all data for temporal analysis
                temporal_data = self.extract_temporal_data(df)
            else:
                window_data = self.extract_window_data(df, window_index)
                temporal_data = self.extract_temporal_data(window_data)
            
            folder_data[folder] = temporal_data
        
        if folder_data:
            self.calculate_temporal_statistics(folder_data, rolling_window)
            self.plot_temporal_evolution(folder_data, rolling_window)
            self.plot_trend_analysis(folder_data, rolling_window)
        
        return folder_data
    
    def analyze_cross_window_temporal(self, folders=['12_logs', '24_logs']):
        """Analyze temporal changes across different windows"""
        print("\n" + "="*60)
        print("CROSS-WINDOW TEMPORAL ANALYSIS")
        print("="*60)
        
        for folder in folders:
            print(f"\nAnalyzing temporal patterns across windows in {folder}:")
            df = self.load_log_files(folder)
            if df.empty:
                continue
            
            window_metrics = {}
            for window in sorted(df['window_index'].unique()):
                window_data = df[df['window_index'] == window]
                
                # Calculate window-level metrics
                avg_anonymity = window_data['anonymity_set_size'].mean()
                avg_prob = window_data['correct_batch_prob'].mean()
                success_rate = window_data['correct_batch_is_highest'].mean()
                unique_fraction = (window_data['anonymity_set_size'] == 1).mean()
                
                window_metrics[window] = {
                    'avg_anonymity': avg_anonymity,
                    'avg_prob': avg_prob,
                    'success_rate': success_rate,
                    'unique_fraction': unique_fraction
                }
                
                print(f"  Window {window}: Anonymity={avg_anonymity:.2f}, Success={success_rate:.3f}, Unique={unique_fraction:.3f}")
            
            # Plot window-level trends
            windows = sorted(window_metrics.keys())
            if len(windows) > 2:
                self.plot_window_trends(folder, windows, window_metrics)
    
    def plot_window_trends(self, folder_name, windows, window_metrics):
        """Plot trends across windows"""
        self.setup_plot_style(figsize=(12, 8))
        
        anonymity_values = [window_metrics[w]['avg_anonymity'] for w in windows]
        success_values = [window_metrics[w]['success_rate'] for w in windows]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot anonymity across windows
        ax1.plot(windows, anonymity_values, 'bo-', linewidth=2, markersize=8)
        ax1.set_xlabel('Window Index')
        ax1.set_ylabel('Average Anonymity Set Size')
        ax1.set_title(f'Anonymity Evolution Across Windows\n({folder_name})')
        ax1.grid(True, alpha=0.3)
        
        # Plot success rate across windows
        ax2.plot(windows, success_values, 'ro-', linewidth=2, markersize=8)
        ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, label='Random Guess')
        ax2.set_xlabel('Window Index')
        ax2.set_ylabel('Adversary Success Rate')
        ax2.set_title(f'Privacy Degradation Across Windows\n({folder_name})')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)
        
        plt.tight_layout()
        filename = self.generate_filename(f"window_trends_{folder_name}")
        self.save_plot(filename)
        plt.show()

def main():
    analyzer = TemporalAnalysis()
    
    # Analyze temporal changes
    analyzer.analyze_temporal_changes()
    
    # Analyze cross-window temporal patterns
    analyzer.analyze_cross_window_temporal()

if __name__ == "__main__":
    main()