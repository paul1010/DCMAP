# DCMAP

Optimal Clustering with Dependent Costs in Bayesian Networks
https://doi.org/10.1613/jair.1.23620 

Paul Pao-Yen Wu
https://orcid.org/0000-0001-5960-8203
Fabrizio Ruggeri
https://orcid.org/0000-0003-0615-4417
Kerrie Mengersen
https://orcid.org/0000-0001-8625-9168
Abstract

Background: Clustering of nodes in Bayesian Networks (BNs) and related graphical models such as Dynamic BNs (DBNs) has been demonstrated to enhance computational efficiency and improve model learning. It typically involves partitioning the underlying Directed Acyclic Graph (DAG) into cliques or optimising for some cost or criteria.

Objectives: Given a DAG, we focus on a critical but understudied aspect of optimal clustering involving cost dependency. This is where inference outcomes and hence clustering costs depend on both nodes within a cluster and the mapping of clusters that are connected by at least one arc.

Methods: We propose a novel algorithm called Dependent Cluster MAPping (DCMAP) which can, given an arbitrary, positive cost function, iteratively and rapidly find near-optimal, then optimal cluster mappings.

Results: DCMAP is shown analytically to be optimal in terms of finding all of the least cost cluster mapping solutions and with no more iterations than an equally informed algorithm. Demonstrated on a complex systems seagrass DBN with 9.91 × 109 and 1.51 × 1021 possible cluster mappings for 25 and 50 node configurations, it took 856 and 1569 iterations on average to find the first optimal solution, respectively.

Conclusions: The effectiveness of DCMAP enables future research in BN learning using optimisation, such as through enhancing computational efficiency or minimising entropy for learning. This is critically important as computation of marginal distributions or updating model parameters is NP-hard.


