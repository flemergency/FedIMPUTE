# FedIMPUTE

### R package for preprint [FedIMPUTE: Privacy-Preserving Missing Value Imputation for Multi-site Heterogeneous Electronic Health Records](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4930174).

<p align="left">
  <img height="420" src="Workflow.jpg">
</p>

We propose FedIMPUTE, a communication-efficient federated learning (FL) based approach for missing value imputation (MVI). Our method enables multiple sites to collaboratively perform MVI in a privacy-preserving manner, addressing challenges of data-sharing constraints and population heterogeneity. We begin by conducting MVI locally at each participating site, followed by the application of various FL strategies, ranging from basic to advanced, to federate local MVI models without sharing site-specific data. The federated model is then broadcast and used by each site for MVI. We evaluate FedIMPUTE using both simulation studies and a real-world application on electronic health records (EHRs) to predict emergency department (ED) outcomes as a proof of concept. Simulation studies show that FedIMPUTE outperforms all baseline MVI methods under comparison, improving downstream prediction performance and effectively handling data heterogeneity across sites. By using ED datasets from three hospitals within the Duke University Health System (DUHS), FedIMPUTE achieves the lowest mean squared error (MSE) among benchmark MVI methods, indicating superior imputation accuracy. Additionally, FedIMPUTE provides good downstream prediction performance, outperforming or matching other benchmark methods. FedIMPUTE enhances the performance of downstream risk prediction tasks, particularly for sites with high missing data rates and small sample sizes. It is easy to implement and communication-efficient, requiring sites to share only non-patient-level summary statistics.
