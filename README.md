# CausalPath
This is a tool for pathway analysis of proteomics and phosphoproteomics datasets. CausalPath aims to identify all pathway relations that can explain observed correlations in experiments. 

We provide this tool as a web service in the below link:
[todo]

Alternatively, users can locally run this Java code on their datasets, and visualize the results in [ChiBE](https://github.com/PathwayCommons/chibe). To get a copy of the project and build it, first make sure you have [git](https://git-scm.com/downloads), [maven](https://maven.apache.org/download.cgi) and [JDK8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) in your system, then while in your projects folder, do:
```
git pull https://github.com/PathwayAndDataAnalysis/causalpath.git
cd causalpath
mvn clean install
mvn assembly:single
```

