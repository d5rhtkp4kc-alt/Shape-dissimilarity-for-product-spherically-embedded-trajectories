---------------- Data -----------------

Data in the 'energy' folder is the global Energy Consumption Composition (ECC) dataset from the Energy Institute, and is downloaded via https://www.energyinst.org/statistical-review/resources-and-data-downloads.

Data for male Life-Table Death Counts (LTDC) should be downloaded from the Human Mortality Database via https://www.mortality.org. 

Data in the 'gait' folder is the the Multivariate Gait (MultiGait) dataset, and is downloaded via https://archive.ics.uci.edu/dataset/760/multivariate+gait+data.

---------------- Figures -----------------

Run 'visual_plot.R' for Figures 1, 2 and 3.

---------------- Experiments -----------------

STEP I. Results for dissimilarity matrices.

Synthetic Data: First run 'sim_run.R' which reports the results for k-NN using dissimilarity matrices and generates synthetic datasets for neural-network-based models.
Real-World Data: First run 'real_data.R' which reports the results for k-NN using dissimilarity matrices and generates real-world datasets for neural-network-based models.

STEP II. Results for neural-network-based models.

Run 'cnn.ipynb' for CNN after generating required datasets in STEP I.
Run 'manifold.ipynb' for ManifoldNet after generating required datasets in STEP I.
