# HBic Algorithm
**``HBic``** is a biclustering algorithm for heterogeneous and missing data.``HBic`` handles mixed-type data including numerical, binary, and categorical attributes. This is the source code of ``HBic`` developed in `MATLAB R2020b`. In addition, the ``Python`` implementation is available at [**``py-hbic``**](https://github.com/ClementChauvet/py-hbic/).

**``HBic``** natively handles mixed datasets having multiple mixed-data types. Some of the main characteristics of  data ``HBic`` are:
+ A fitness function is proposed for evaluating biclusters with mixed-type attributes and missing values.
+ A model selection approach determines the most relievable biclusters based on their similarity.
+  ``HBic`` automatically identifies the number of biclusters, or takes as input this parameter if this knowledge is available.

----

**``HBic``** is described in detail in our paper:
```
Adán José-García, Julie Jacques, Clement Chauvet, Vincent Sobanski, and Clarisse Dhaenes  
A Biclustering Approach for Heterogeneous and Missing Data
Submitted to: Pattern Recognition.
```

----

**``Getting Started``**

>**HBic** was developed with MATLAB. To try the algorithm look at the scripts `demo_heterogeneous_data.m` and `demo_numerical_data.m`.

----
**``Contact me``**
```
Adán José-García (adan.jose-garcia@univ-lille.fr)
```
