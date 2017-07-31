# med4way
### A Stata command for the 4-way decomposition using parametric regression models

- Current version: `2.2.0` 
- Release date: `31jul2017`

---

### Description

`med4way` uses parametric regression models to estimate the components of the 4-way decomposition of the total effect of a treatment on a outcome in the presence of a mediator with which the exposure may interact. This decomposition breaks down the total effect of the treatment on the outcome into components due to mediation alone, to interaction alone, to both mediation and interaction, and to neither mediation nor interaction.

`med4way` provides standard errors and confidence intervals for the estimated components using the delta method (default) or the bootstrap.

`med4way` allows continuous, binary, count or survival outcomes, and continuous or binary mediators. 

Note: the 4-way decomposition holds without any assumptions about confounding. However, to interpret each of the components causally does require assumptions about confounding. See VanderWeele (2014) for a detailed exposition of those assumptions.


### Installation

- To install the current version of `med4way` directly from GitHub, run:
```Stata
net install med4way, from("https://raw.githubusercontent.com/anddis/med4way/master/")
```
from within a web-aware Stata (version 13+).

- For older versions of Stata, download and extract the [zip file](https://github.com/anddis/med4way/archive/master.zip) and then run:
```Stata
net install med4way, from(mydir)
```
from within Stata, where *mydir* is the directory that containes the extracted files.

- After installation, see the help file:
```Stata
help med4way
```
- To download in the current working directory the datasets needed to run the example code in the help file, type:
```Stata
net get med4way, from("https://raw.githubusercontent.com/anddis/med4way/master/")
```

### Authors

Andrea Discacciati (1), Andrea Bellavia (2,3), Linda Valeri (4,5)

*(1) Unit of Biostatistics, Karolinska Institutet, Stockholm, Sweden (2) Department of Environmental Health, Harvard T.H. Chan School of Public Health, Boston, MA (3) Department of Biostatistics, Harvard T.H. Chan School of Public Health, Boston, MA (4) Department of Psychiatry, Harvard Medical School, Boston, MA (5) Psychiatric Biostatistics Laboratory, McLean Hospital, Belmont, MA*

### References

VanderWeele, T.J., 2014. A unification of mediation and interaction: a 4-way decomposition. Epidemiology (Cambridge, Mass.), 25(5), p.749.