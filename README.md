

# Monkeypox Knowledge Graph

The Monkeypox KG is built using viral and human proteins reported in different resources. Additionally, the KG represents chemicals tested against Monkeypox and their targets, associated biological processes, molecular functions, diseases and side effects. 

**KG status**

Version 1 stats:

* Number of Nodes: 8235
* Number of Edges: 40422

Version 2 stats (2nd September) :

* Number of Nodes: 9129
* Number of Edges: 44568

# Running the KG 

### Pre-requisite: Installation of python 3 and jupyter notebook (Please skip it if you already have these in your pc)

We recommend installing anaconda (https://www.anaconda.com/) which will automatically install python and jupyter notebook. 

a) Next step would be to create a dedicated environment within anaconda to run the scripts. To do this, open anaconda. You'll be able to see HOME at the left of the screen.

b) Underneath HOME, you'll see Environments. Click it. 

c) At the bottom of the page, you'll see 'Create' with a '+' sign. Click it. 

d) You can now give a name for you environment, say: mpox. Switch the version of python to 3.10 and click create. It will take few minutes.

e) Once done, you can now click 'Home' and you should see in the middle of screen that 'Applications on' will be 'mpox'. However, please remember that next time you open anaconda, you will have to select 'mpox' again because the original environment will be set to 'base (root)'.

Please check the graph.ipynb for understanding step wise process of KG generation and KG statistics. The KG has been exported to formats such as graphml, sif and so on for visualizations in other platforms. For example, the graphml file can be imported to Cytoscape directly. These files are located in 'data\export' folder.

# Dissemination and deployment

1) The KG is registered in EMBL-EBI's BioModels (https://www.ebi.ac.uk/biomodels/MODEL2208040001). 

2) A manuscript for this work is currently under review in a peer-reviewed journal. However, you can read the manuscript in bioRxiv (https://www.biorxiv.org/content/10.1101/2022.08.02.502453v1)

3) The KG is hosted publicly at the NDEx server: https://doi.org/10.18119/N9SG7D 
