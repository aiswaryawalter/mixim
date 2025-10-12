import os
import sys

# Add analyzers directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'analyzers'))

from analyzers.anonymity_analysis import AnonymityAnalysis
from analyzers.probability_analysis import BatchProbabilityAnalysis
from analyzers.accuracy_analysis import AccuracyAnalysis
from analyzers.temporal_analysis import TemporalAnalysis
from analyzers.unique_identification_analysis import UniqueIdentificationAnalysis



# Import other analyses as you create them

def main():
    print("="*80)
    print("RUNNING ALL BATCH LOG ANALYSES")
    print("="*80)
    
    folders = ['12_logs', '24_logs']
    
    # 1. Anonymity Set Size Analysis
    print("\n1. ANONYMITY SET SIZE ANALYSIS")
    print("-" * 50)
    anonymity_analyzer = AnonymityAnalysis()
    anonymity_analyzer.analyze_anonymity_sets(folders)
    
    # 2. Batch Probability Analysis
    print("\n2. BATCH PROBABILITY ANALYSIS")
    print("-" * 50)
    probability_analyzer = BatchProbabilityAnalysis()
    probability_analyzer.analyze_batch_probabilities(folders)

    # 3. Adversary Accuracy Analysis
    print("\n3. ADVERSARY ACCURACY ANALYSIS")
    print("-" * 50)
    accuracy_analyzer = AccuracyAnalysis()
    accuracy_analyzer.analyze_adversary_accuracy(folders)

    # 4. Temporal Analysis
    print("\n4. TEMPORAL ANALYSIS")
    print("-" * 50)
    temporal_analyzer = TemporalAnalysis()
    temporal_analyzer.analyze_temporal_changes(folders)

    # 5. Unique Identification Analysis
    print("\n5. UNIQUE IDENTIFICATION ANALYSIS")
    print("-" * 50)
    unique_analyzer = UniqueIdentificationAnalysis()
    unique_analyzer.analyze_unique_identification(folders)
    
    print("\n" + "="*80)
    print("ALL ANALYSES COMPLETED")
    print("="*80)
    
    # List generated diagrams
    diagrams_path = os.path.join(os.path.dirname(__file__), "diagrams")
    if os.path.exists(diagrams_path):
        diagram_files = [f for f in os.listdir(diagrams_path) if f.endswith('.png')]
        print(f"\nGenerated {len(diagram_files)} diagrams:")
        for file in sorted(diagram_files)[-10:]:  # Show last 10
            print(f"  - {file}")

if __name__ == "__main__":
    main()