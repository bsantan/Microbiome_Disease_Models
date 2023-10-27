# Microbiome Disease Models

This repository hosts models of known influence of the gut microbiome in disease. All concepts are rooted in biological ontologies such that the content can be easily integrated into other knowledge bases, reviewed by experts, and accurate in depicting biological content.

## Getting Started

These instructions will provide the necessary environments, programs with installation instructions, and input files in order to create a visualization of any given mechanism in this repository using Cytoscape. 

### Dependencies
The following dependencies are listed in the environment.yml file, and installed in the installation step. This software has only been tested on Unix based OS systems, not Windows.
```
Python>=3.8.3
py4cytoscape>=1.3.0
```

## Installation

```
git clone https://github.com/bsantan/Microbiome_Disease_Models
```

First install mamba, which will be used to create the environment. To create an environment with all dependencies and activate the environment, run the following commands:

```
conda install mamba

mamba env create -f environment.yml
conda activate Microbiome_Disease_Models
```

Ensure that Cytoscape (any version later than 3.8.0) is up and running before continuing.

## Running the Script

From the /scripts folder execute: 

```
python visualize_mechanism.py --input-mech INPUTMECH --labels-file LABELSFILE
```

InputMech: "<directory>/<input_mech_filename>.csv"
LabelsFile: "<nodes>.csv"

## Expected Outputs
 
The visualize_mechanism.py script will always generate the following files with the same filename string as found in the original InputMech:
  
### Subgraph Attributes

A .noa file which specifies the category of each node as specified in the LabelsFile. These are expected to be in biolink format.
  
```
Node|Attribute
insulin receptor (human)|Protein
Parkinson disease|Disease
```

### Subgraph Visualization
  
A .png file generated in cytoscape with all nodes found in the path search, colored by the node category.
