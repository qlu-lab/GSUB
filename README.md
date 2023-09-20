# GSUB
GSUB (Gwas-by-SUBtraction) is a command line tool to perform GWAS-by-subtraction. It uses closed-form solutions to estimate the associations between each SNP and the two latent components for two traits (i.e., $\gamma_1$ and $\gamma_2$ in the diagram below). Its inputs are two GWAS summary statistics for two traits, and the outputs are two summary statistics for the two latent factors. The $\gamma_2$ is shown in red is because it is the parameter of interest: the SNP effect on the latent factor for trait 2 after regressing out the factor for trait 1.

![GSUB workflow](https://github.com/qlu-lab/GSUB/blob/main/figures/GSUB_general_workflow.png)

## Manual
Please see the [wiki](https://github.com/qlu-lab/GSUB/wiki) for the detailed manual of `GSUB`.

## Version History
* September 19, 2023: Initial release.

## Citation
If you use GSUB, please cite:
Coming soon...

## Contact
For questions and comments, please open a Github issue (preferred) or contact Yuchang Wu at [ywu423@wisc.edu](mailto:ywu423@wisc.edu) or Stephen Dorn at [svdorn@wisc.edu](mailto:svdorn@wisc.edu).
