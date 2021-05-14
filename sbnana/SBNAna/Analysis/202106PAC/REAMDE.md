# 2021 May PAC Documentation 

Please edit this readme file with the following information as you push your changes: 
- Instructions (including datasets used, versioning, etc.) to reproduce this analysis. 
- For fitter scripts which are not yet in this repository, please point to relevant code here. 
- [DocDB](https://sbn-docdb.fnal.gov/cgi-bin/sso/DocumentDatabase) references or references to other documentation where relevant.

*Please help other developers by merging to develop every day - couple of days in the weeks leading up to the PAC.*


### Instructions for Event Selection Contributions


*Standard Cuts*

If your cuts and vars are standard i.e. "Cut all events with 0 slices" please leave them
in [`SBNAna/Cuts`](../../Cuts) and [`SBNAna/Vars`](../../Vars).

*Final cuts for the PAC* 

The version of the analysis specific cuts you use here to `SBNAna/Cuts/NumuCuts202106.h` etc. 

*Macros*

Please commit the scripts you use to produce plots for the 2021 June PAC to the [`SBNAna/Analysis/202106PAC`](.) directory. Fitter scripts can also go here. 


