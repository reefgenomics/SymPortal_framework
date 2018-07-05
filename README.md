# README
Welcome to the SymPortal GitHub repository. SymPortal is an analytical tool used to resolve _Symbiodinium_ spp. taxa using NGS data of the ITS2 marker gene.

Our methodology augments tried-and-tested intragenomic resolution theory with the power of NGS to resolve between symbiont taxa at a level far surpassing alternative methodologies using this marker. 
We employ novel logic to identify within-sample informative intragenomic sequences, which we have termed defining intragenomic variants (DIVs), and use combinations of these DIVs to identify ITS2 type profiles representative of putative _Symbiodinium_ taxa.

The SymPortal analytical framework consists of two parts: 
the SymPortal analysis and the SQL database with which the analysis is integrated. 

The SymPortal analytical framework may be run either locally, by running the python scripts housed on this GitHub repository 
or remotely, through submission of data to [SymPortal.org](http:symportal.org). 
Analyses run remotely through SymPortal.org will have access to the latest version of the SymPortal database. 
This remotely hosted PostgreSQL database will contain sequencing information from all samples previously analysed that were submitted via SymPortal.org. 
Analyses run via [SymPortal.org](http:symportal.org) will therefore have access to the greatest possible resolving power and enable comparability to other datasets previously submitted to SymPortal.org. 
If running the SymPortal framework locally, it will be necessary for the user to populate the local database with which the local SymPortal analysis will integrate. 
As such, the ability to resolve ITS2 type profiles will be contingent on the extent of sequencing information housed in the local database when running analyses locally.
Further documentation may be found at the SymPortal GitHub repository.

Please visit the [SymPortal wiki](https://github.com/SymPortal/SymPortal_framework/wiki) for further details on setting up a local instance of the SymPortal analytical framework.
