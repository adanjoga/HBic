# HBic Algorithm
**``HBic``** is a biclustering algorithm for heterogeneous and missing data.``HBic`` handles mixed-type data, including numeric, binary, and categorical attributes. This is the source code of ``HBic`` developed in `MATLAB R2020b`. In addition, the ``Python`` implementation is available at [**``py-hbic``**](https://github.com/ClementChauvet/py-hbic/).

**``HBic``** natively handles mixed datasets with multiple mixed-data types. Some of the main characteristics of  data ``HBic`` are:
+ A fitness function is proposed for evaluating biclusters with mixed-type attributes and missing values.
+ A model selection approach determines the most relievable biclusters based on their similarity.
+  ``HBic`` automatically identifies the number of biclusters or takes this parameter as input if this knowledge is available.

----

**``HBic``** is described in detail in our paper:
```
Adán José-García, Julie Jacques, Clement Chauvet, Vincent Sobanski, and Clarisse Dhaenes  
HBIC: A Biclustering Approach for Heterogeneous Datasets
To be published in the 27TH EUROPEAN CONFERENCE ON ARTIFICIAL INTELLIGENCE.
https://www.ecai2024.eu/
```

----

**``Getting Started``**

>**HBic** was developed with MATLAB. To try the algorithm, look at the scripts `demo_heterogeneous_data.m` and `demo_numerical_data.m`.

----
**``Contact us``**
```
Adán José-García (adan.jose-garcia@univ-lille.fr)
```
