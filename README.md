# ME-ICA_project
Task-Based Functional Connectivity for ME-ICA:
Can ME-ICA generate a consistent subject-level connectivity matrix for individuals with high task-correlated motion that is robust enough for clinical diagnostic use, such as in chronic stroke patients? 

Using openneuro for the dataset here: https://github.com/OpenNeuroDatasets/ds007661 / https://openneuro.org/datasets/ds007661/versions/1.0.0 

To see the full and final code --> Go into the "notebooks" file and check out and download the R markdown (html/Rmd) files

Our aim is to address the following: --> 
1. Generate and examine group-level functional connectivity matrix from 8 subjects during a hand-grasp task by examining connectivity between the primary motor cortex (precentral gyrus/M1) and the rest of the brain. 
2. Identify the top ten ROI’s showing the strongest functional connectivity with M1.
3. To establish a healthy reference pattern of motor-network connectivity that may inform future comparisons with clinical populations affected by motor dysfunction, such as Parkinson’s disease and stroke.

Project idea by: John Lio, Arshiya Khurmi, Kelei Xiao, Pirinthiya Thayaparan

Project Definition

Background

Our group
We are a group of undergraduate students conducting summer projects at the Centre for Addiction and Mental Health (CAMH) as a part of the Institute of Medical Science at the University of Toronto’s Summer Undergraduate Research Program (IMS SURP). Collectively, our summer projects focus on the various neurobiological mechanisms and outcomes associated with psychiatric conditions.

Project background
Motor task functional magnetic resonance imaging (fMRI) is widely used to investigate motor function in neurological disorders such as Stroke and Parkinson's disease. However, these populations often exhibit increased head motion during motor tasks, producing signals which can obscure true neural activity, subsequently reducing the reliability of neuroimaging findings. A particular challenge is task-correlated motion, which occurs when head movement coincides with task performance, such as during hand-grasp movements, making it difficult to distinguish motion-related artifacts from task-associated neural activity. Multi-Echo Independent Component Analysis (ME-ICA) addresses this issue by combining multi-echo fMRI acquisition with independent component analysis (ICA) to separate blood oxygenation level-dependent (BOLD) signals from non-BOLD sources of noise, including head motion and physiological fluctuations. By removing noise-related components while preserving neural signals, ME-ICA may improve the accuracy and interpretability of motor-task fMRI analyses in populations prone to excessive movement.

Tools
This project utilizes a variety of tools to support development and open access research practices -->
RStudio: Core environment for building connectivity matrices and seed-based analysis
Terminal: Unlocking the dataset from openneuro, and pushing the project onto Github
Github: Communication and storage of project files

Data
This project uses fMRI data from a study completed by Reddy et al., which aims to compare ME-ICA with other denoising tools in addressing noise associated with task-correlated motion. (https://doi.org/10.1101/2023.07.19.549746)

Deliverables
This project will produce a GitHub repository containing all figures, connectivity matrices, ROI analyses, and R markdown files.

Results

Progress overview
The final presentation of this project was delivered on June 5, 2026. All planned deliverables were completed and are available in the project’s GitHub repository.

Tools we learned during this project
RStudio & R: Throughout the course, we were grateful to have had hands-on experience using RStudio in processing neuroimaging data to produce connectivity matrices and ROI analyses.
Terminal & Bash: We primarily learned how to use datalad to store larger datasets, which was particularly helpful as our computer storage was limited.
Github & Git: We learned how to use git push, git pull, and git commit to keep our project organized.

Deliverables
Our main deliverables are available in a GitHub repository containing all figures, connectivity matrices, ROI analyses, and R markdown files. (https://github.com/brainhack-school2026/ME-ICA_project)

Conclusion

Acknowledgements
We would like to thank BrainHack School 2026 for providing us with immense amounts of support throughout our project. Particularly we are grateful for the TA team, who were particularly compassionate in their feedback and mentorship as we navigated the learning experience.

References
Esteban, O., Markiewicz, C.J., Blair, R.W. et al. fMRIPrep: a robust preprocessing pipeline for functional MRI. Nat Methods 16, 111–116 (2019). https://doi.org/10.1038/s41592-018-0235-4
Neha A Reddy, Kristina M Zvolanek, and Molly G Bright (2026). Motor Task Multi-Echo fMRI - fMRIPrep Derivatives. OpenNeuro. [Dataset] doi: https://doi.org/10.18112/openneuro.ds007661.v1.0.0
Neha A. Reddy, Kristina M. Zvolanek, Stefano Moia, César Caballero-Gaudes, Molly G. Bright bioRxiv 2023.07.19.549746; doi: https://doi.org/10.1101/2023.07.19.549746
Reddy, N. A., Zvolanek, K. M., Moia, S., Caballero-Gaudes, C., & Bright, M. G. (2023). Denoising task-correlated head motion from motor-task fMRI data with multi-echo ICA. bioRxiv. https://doi.org/10.1101/2023.07.19.549746

