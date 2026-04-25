#!/bin/bash
# Migration script — apply v2 restructuring to existing local repo
# Run from repo root: bash migrate_v2.sh

set -e

echo "=== Step 1: Restructure feature_selection vs exploration ==="
mkdir -p scripts/exploration
git mv scripts/feature_selection/PPI_meth.R         scripts/exploration/ 2>/dev/null || true
git mv scripts/feature_selection/PPI_Meth_2.R       scripts/exploration/ 2>/dev/null || true
git mv scripts/feature_selection/PPI_Meth_3.R       scripts/exploration/ 2>/dev/null || true
git mv scripts/feature_selection/stringdb_train.R   scripts/exploration/ 2>/dev/null || true
git mv scripts/feature_selection/sparse_beta_label.R scripts/exploration/ 2>/dev/null || true
git mv scripts/feature_selection/label.txt          scripts/exploration/ 2>/dev/null || true

echo "=== Step 2: Rename notebooks ==="
cd notebooks
git mv compiled_model.ipynb               01_model_comparison_main.ipynb 2>/dev/null || true
git mv SVM.ipynb                          02a_svm_development.ipynb 2>/dev/null || true
git mv RF.ipynb                           02b_random_forest_development.ipynb 2>/dev/null || true
git mv NN.ipynb                           02c_neural_network_development.ipynb 2>/dev/null || true
git mv Performance_Visualization.ipynb     03_performance_visualization.ipynb 2>/dev/null || true
git mv shap_test.ipynb                    04_shap_biomarker_extraction.ipynb 2>/dev/null || true
git mv Gene_Biomarker_Visualization.ipynb  05_biomarker_visualization.ipynb 2>/dev/null || true
git mv GO_plot.ipynb                      06_go_enrichment_plot.ipynb 2>/dev/null || true
git mv clustering_visualization.ipynb      07_clustering_visualization.ipynb 2>/dev/null || true
cd ..

echo "=== Step 3: Rename pathway scripts ==="
cd scripts/pathway
git mv kegg_1.R 01_combine_shap_features.R 2>/dev/null || true
git mv kegg_2.R 02_build_rank_set.R 2>/dev/null || true
git mv kegg_3.R 03_normalize_rank.R 2>/dev/null || true
git mv kegg_4.R 04_pathway_enrichment.R 2>/dev/null || true
git mv kegg_5.R 05_kegg_view_mapping.R 2>/dev/null || true
git mv kegg_6.R 06_pathview_visualization.R 2>/dev/null || true
cd ../..

echo "=== Step 4: Done. Now copy these new files into the repo: ==="
echo "  - README.md (updated)"
echo "  - requirements.txt (new)"
echo "  - environment.yml (new)"
echo "  - results/figures/pipeline_diagram.png (new)"
echo ""
echo "Then run: git add . && git commit -m 'Refactor v2: cleaner structure, biomarker-focused' && git push"
