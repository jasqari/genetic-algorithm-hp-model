# HP Protein Folding Model using Genetic Algorithm
<p align="center">
 <img src="https://github.com/jasqari/genetic-algorithm-hp-model/assets/44480584/17b02515-f662-4639-8b13-5a470f1fa0c5"/>
</p>

## Overview
This repository contains Python codes that collectively aim to find optimal foldings for protein sequences using the concepts of 2D hydrophobic-polar protein folding model and genetic algorithm.
* `hp_model.py` contains the tools necessary for the 2D HP protein folding model,
* `genetic_algorithm.py` is an implementation of the genetic algorithm suitable for the problem at hand,
* `main.py` is a driver program that utilizes all the tools to output solutions.

## Requirements
```
pip3 install -r requirements.txt
```

## Usage
Have a look at the options and their default values:
```
python3 main.py -h
```

Provide a protein sequence either by typing it in or by choosing from `sample_seqs.txt`:
```
python3 main.py LHAKGMHG
```
```
python3 main.py 5
```
Take advantage of `sample_seqs.txt` by augmenting it with any long sequences you may want to input.

Control the search for the optimal folding using the following parameters:
```
python3 main.py 5 --pop_size 1000 --num_gens 500 --epsilon 1e-5 --energy_func custodio
```
The quality of the search results relies heavily on the values of these parameters, especially `--pop_size` and `--energy_func`.
For shorter sequences lower the `--pop_size` as the algorithm won't be able to generate enough unique solutions and as a result, the search process will be stopped. Alternatively, change the algorithm to include duplicate solutions in the population, which may in turn hurt the quality of search results.
No need to mention, that for a longer sequence, a large enough population size and a lot of patience are necessary.
`berger` energy function [1](https://dl.acm.org/doi/10.1145/279069.279080) is relatively simpler than `custodio` [2](https://www.scielo.br/j/gmb/a/nbhndMBBG5z4y39hyZ7463f/?format=html&lang=en) and might work better in some scenarios. However, in most cases, `custodio` outperformed `berger` as a fitness function. Moreover, you may customize the default weights involved in `custodio` according to your needs and further enhance the search results.

Finally, log the optimization process and graphically display the (locally) optimal protein folding:
```
python3 main.py 3 --log --display
```

Solutions will be saved to `optimal_folding.txt`.

