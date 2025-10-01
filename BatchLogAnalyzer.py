import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import numpy as np
from pathlib import Path
from datetime import datetime

class BatchLogAnalyzer:
    def __init__(self, base_path="/Users/aiswaryawalter/Desktop/mixim"):
        self.base_path = base_path
        self.diagrams_path = os.path.join(base_path, "diagrams")
        
        # Create diagrams folder if it doesn't exist
        os.makedirs(self.diagrams_path, exist_ok=True)
        print(f"Diagrams will be saved to: {self.diagrams_path}")
    
    def get_timestamp(self):
        """Generate timestamp for unique filenames"""
        return datetime.now().strftime("%Y%m%d_%H%M%S")
    
    def generate_filename(self, base_name, extension=".png"):
        """Generate unique filename with timestamp"""
        timestamp = self.get_timestamp()
        return f"{base_name}_{timestamp}{extension}"
        
    def load_log_files(self, folder_name):
        """Load all CSV files from a specific log folder"""
        folder_path = os.path.join(self.base_path, folder_name)
        csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
        
        all_data = []
        for file in csv_files:
            try:
                df = pd.read_csv(file)
                df['source_file'] = os.path.basename(file)
                df['folder'] = folder_name
                all_data.append(df)
                print(f"Loaded {len(df)} rows from {os.path.basename(file)}")
            except Exception as e:
                print(f"Error loading {file}: {e}")
        
        if all_data:
            combined_df = pd.concat(all_data, ignore_index=True)
            print(f"Total loaded from {folder_name}: {len(combined_df)} rows")
            return combined_df
        else:
            print(f"No data loaded from {folder_name}")
            return pd.DataFrame()
    
    def extract_window_data(self, df, window_index=None):
        """Extract data from a specific window or the largest window."""
        if df.empty:
            return df
            
        if window_index is None:
            max_window = df['window_index'].max()
            print(f"Using largest window index: {max_window}")
            window_data = df[df['window_index'] == max_window].copy()
        else:
            window_data = df[df['window_index'] == window_index].copy()
            print(f"Using specified window index: {window_index}")
        
        print(f"Extracted {len(window_data)} rows from window {window_data['window_index'].iloc[0] if not window_data.empty else 'N/A'}")
        return window_data
    
    def save_plot(self, filename, dpi=300):
        """Save plot with consistent settings"""
        save_path = os.path.join(self.diagrams_path, filename)
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
        return save_path
    
    def setup_plot_style(self, figsize=(14, 8)):
        """Setup consistent plot styling"""
        plt.style.use('default')
        plt.figure(figsize=figsize)
        return plt.gca()