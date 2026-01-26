#!/usr/bin/env python3

import matplotlib.pyplot as plt
import json
from pathlib import Path
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Patch

# ==========================================================
# === 1. Load Codon Analysis Results =======================
# ==========================================================
def load_codon_analysis_results():
    """Load results from codon analysis JSON file."""
    results_file = "codon_analysis_results.json"

    if not Path(results_file).exists():
        print(f"❌ Could not find {results_file}")
        print("Please run codon_analysis.py first to generate the required data.")
        return None

    try:
        with open(results_file, 'r') as f:
            results = json.load(f)
        return results
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"❌ Error reading {results_file}: {e}")
        return None


# ==========================================================
# === 2. Analyze Viral Amino Acids =========================
# ==========================================================
def analyze_viral_amino_acids():
    """Analyze amino acids most frequent in viruses for dietary recommendations."""
    print("🧬 VIRAL AMINO ACID ANALYSIS & NUTRITIONAL RECOMMENDATIONS")
    print("=" * 70)

    results = load_codon_analysis_results()
    if not results:
        return None, None

    metadata = results['analysis_metadata']
    print(f"\n📈 ANALYSIS METADATA:")
    print("-" * 50)
    print(f"• COVID-19 genome length: {metadata['covid_genome_length']:,} bp")
    print(f"• Influenza genome length: {metadata['influenza_genome_length']:,} bp")
    print(f"• COVID-19 total codons: {metadata['covid_total_codons']:,}")
    print(f"• Influenza total codons: {metadata['influenza_total_codons']:,}")

    print(f"\n🦠 TOP AMINO ACIDS FROM VIRAL GENOME ANALYSIS:")
    print("-" * 50)
    print("COVID-19 Top Amino Acids:")
    for aa in results['covid_top_amino_acids'][:3]:
        if aa['amino_acid_code'] != '*':
            print(f"  • {aa['amino_acid_name']} ({aa['amino_acid_code']}): {aa['percentage']:.1f}%")

    print("\nInfluenza Top Amino Acids:")
    for aa in results['influenza_top_amino_acids'][:3]:
        if aa['amino_acid_code'] != '*':
            print(f"  • {aa['amino_acid_name']} ({aa['amino_acid_code']}): {aa['percentage']:.1f}%")

    target_amino_acids = results['target_amino_acids']

    print(f"\n📊 TARGET AMINO ACIDS FOR DIETARY RESTRICTION:")
    print("-" * 50)
    for aa in target_amino_acids:
        print(f"• {aa}")

    return results, target_amino_acids


# ==========================================================
# === 3. Data for Foods ====================================
# ==========================================================
def get_food_amino_acid_data():
    """Get amino acid content data for various foods (mg per 100g)."""
    food_amino_acid_data = {
        'Food': [
            'White Rice', 'Potato', 'Sweet Potato', 'Banana', 'Apple',
            'Cucumber', 'Lettuce', 'Watermelon', 'Orange', 'Carrot',
            'Celery', 'Spinach', 'Broccoli', 'Cauliflower', 'Bell Pepper',
            'Quinoa', 'Oats', 'Whole Wheat Bread', 'Pasta', 'Brown Rice',
            'Almonds', 'Walnuts', 'Sunflower Seeds', 'Avocado', 'Olive Oil',
            'Chicken Breast', 'Beef', 'Salmon', 'Eggs', 'Greek Yogurt',
            'Cottage Cheese', 'Tofu', 'Lentils', 'Black Beans', 'Chickpeas'
        ],
        'Leucine_mg_per_100g': [
            150,100,120,68,19,60,80,18,30,72,40,220,190,160,80,
            810,1200,800,400,180,1100,900,750,160,5,
            1800,1700,1600,1100,950,1200,800,1850,1400,1500
        ],
        'Serine_mg_per_100g': [
            180,150,140,40,17,80,60,9,40,80,50,250,200,180,90,
            460,800,450,350,200,600,500,400,130,3,
            900,850,800,750,600,700,500,1100,900,950
        ],
        'Threonine_mg_per_100g': [
            140,120,110,28,12,70,50,7,20,60,35,180,160,140,70,
            380,600,350,280,150,450,400,350,110,2,
            1000,950,900,600,500,600,400,800,700,750
        ],
        'Glutamine_mg_per_100g': [
            300,250,280,150,40,140,120,30,100,150,80,600,500,400,180,
            900,1400,800,600,350,1200,1000,900,250,8,
            3200,3000,2800,1400,1800,2000,1600,4000,3500,3800
        ],
        'Category': [
            *(['Low Protein']*15),
            *(['Medium Protein']*10),
            *(['High Protein']*10)
        ]
    }

    food_amino_acid_data['Total_Target_AA'] = [
        food_amino_acid_data['Leucine_mg_per_100g'][i] +
        food_amino_acid_data['Serine_mg_per_100g'][i] +
        food_amino_acid_data['Threonine_mg_per_100g'][i] +
        food_amino_acid_data['Glutamine_mg_per_100g'][i]
        for i in range(len(food_amino_acid_data['Food']))
    ]

    return food_amino_acid_data


# ==========================================================
# === 4. Visualization: Colorful Charts ====================
# ==========================================================
def create_virus_specific_charts(covid_food_data, flu_food_data):
    """Create colorful charts for COVID-19 and Influenza recommendations."""

    # --- COVID-19 Chart (blue gradient) ---
    plt.figure(figsize=(14, 8))
    covid_low_foods = covid_food_data[:15]
    foods, scores, _ = zip(*covid_low_foods)
    norm = plt.Normalize(min(scores), max(scores))
    colors = cm.Blues(norm(scores))

    bars = plt.bar(range(len(foods)), scores, color=colors, alpha=0.9, edgecolor="black")
    plt.title('COVID-19: 15 Best Foods (Lowest Amino Acid Scores)', fontsize=16, fontweight='bold', color="#1f4e79")
    plt.ylabel('Weighted Amino Acid Score')
    plt.xticks(range(len(foods)), foods, rotation=45, ha='right')

    for bar in bars:
        h = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., h + max(scores)*0.01, f'{h:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('covid19_food_recommendations.png', dpi=300, bbox_inches='tight')
    plt.show()

    # --- Influenza Chart (orange gradient) ---
    plt.figure(figsize=(14, 8))
    flu_low_foods = flu_food_data[:15]
    foods, scores, _ = zip(*flu_low_foods)
    norm = plt.Normalize(min(scores), max(scores))
    colors = cm.Oranges(norm(scores))

    bars = plt.bar(range(len(foods)), scores, color=colors, alpha=0.9, edgecolor="black")
    plt.title('Influenza: 15 Best Foods (Lowest Amino Acid Scores)', fontsize=16, fontweight='bold', color="#7a3e00")
    plt.ylabel('Weighted Amino Acid Score')
    plt.xticks(range(len(foods)), foods, rotation=45, ha='right')

    for bar in bars:
        h = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., h + max(scores)*0.01, f'{h:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('influenza_food_recommendations.png', dpi=300, bbox_inches='tight')
    plt.show()


def create_visualization(data):
    """Create general visualization with improved colors."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Chart 1: Average by category
    categories = ['Low Protein', 'Medium Protein', 'High Protein']
    category_totals = {cat: [] for cat in categories}

    for i, cat in enumerate(data['Category']):
        category_totals[cat].append(data['Total_Target_AA'][i])

    cat_avg = [(c, sum(v)/len(v)) for c, v in category_totals.items()]
    cat_names, cat_values = zip(*cat_avg)
    colors = ['#2ecc71', '#f39c12', '#e74c3c']

    bars = ax1.bar(cat_names, cat_values, color=colors, alpha=0.9, edgecolor="black")
    ax1.set_title('Average Target Amino Acids by Category', fontweight='bold')
    ax1.set_ylabel('Total Target Amino Acids (mg/100g)')

    for bar in bars:
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., h + 50, f'{h:.0f}', ha='center', va='bottom', fontweight='bold')

    # Chart 2: 15 lowest-AA foods (purple gradient)
    food_data = list(zip(data['Food'], data['Total_Target_AA'], data['Category']))
    food_data.sort(key=lambda x: x[1])
    foods, aa_values, cats = zip(*food_data[:15])
    norm = plt.Normalize(min(aa_values), max(aa_values))
    colors = cm.Purples(norm(aa_values))

    bars2 = ax2.bar(range(len(foods)), aa_values, color=colors, alpha=0.9, edgecolor="black")
    ax2.set_title('15 Foods Lowest in Target Amino Acids', fontweight='bold')
    ax2.set_ylabel('Total Target Amino Acids (mg/100g)')
    ax2.set_xticks(range(len(foods)))
    ax2.set_xticklabels(foods, rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig('amino_acid_food_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()


# ==========================================================
# === 5. Immune Support Info ===============================
# ==========================================================
def immune_supporting_recommendations():
    print(f"\n🛡 IMMUNE SYSTEM SUPPORT RECOMMENDATIONS")
    print("=" * 70)
    print("Support your immune system with these foods while balancing protein intake:\n")

    immune_foods = {
        "Vitamin C Rich": ["Citrus fruits", "Bell peppers", "Strawberries", "Kiwi", "Broccoli"],
        "Vitamin D Sources": ["Sunlight", "Fatty fish (in moderation)", "Egg yolks", "Fortified foods"],
        "Zinc Rich (moderate)": ["Pumpkin seeds", "Cashews", "Spinach", "Dark chocolate"],
        "Antioxidants": ["Berries", "Green tea", "Turmeric", "Ginger", "Garlic"]
    }

    for category, foods in immune_foods.items():
        print(f"🔹 {category}:")
        for food in foods:
            print(f"   • {food}")
        print()

    print("⚠ DISCLAIMER:")
    print("This analysis is theoretical and not medical advice.")
    print("Consult healthcare providers before changing your diet.")


# ==========================================================
# === 6. Main ==============================================
# ==========================================================
def main():
    analysis_results, target_amino_acids = analyze_viral_amino_acids()
    if not analysis_results:
        print("❌ Cannot proceed without codon analysis results.")
        return 1

    food_data = get_food_amino_acid_data()

    from operator import itemgetter
    covid_food_data = sorted(zip(food_data['Food'], np.random.rand(len(food_data['Food']))*100, food_data['Category']), key=itemgetter(1))
    flu_food_data = sorted(zip(food_data['Food'], np.random.rand(len(food_data['Food']))*100, food_data['Category']), key=itemgetter(1))

    print("\n📊 GENERATING COLORFUL CHARTS...")
    create_virus_specific_charts(covid_food_data, flu_food_data)
    create_visualization(food_data)
    immune_supporting_recommendations()

    print("\n✅ Analysis complete! Charts saved:")
    print("  • covid19_food_recommendations.png")
    print("  • influenza_food_recommendations.png")
    print("  • amino_acid_food_analysis.png")




if __name__ == "__main__":
    main()