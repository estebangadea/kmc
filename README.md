# Electrochemical kinetic Monte Carlo

kMC is a fix for LAMMPS molecular dynamics that mimics the reaction kinetics of the Butler-Volmer equation

## Installation

Just copy fix_kmc.cpp and .h to your MC folder inside src directory of your LAMMPS build (last version supported lammps-23Jun2022).

## Usage

```python
#set your reactant
variable r_type equal 1

#set your product
variable p_type equal 2

#set your electrode
variable e_type equal 3

#set your pre-exponential factor (Smalley *et. al.* 1995)
variable pexp equal 6e-7

#set your electric potential in correct units
variable pot equal 0.1

#how often will the kMC steps be applied
variable nevery equal 50

#set the fix
fix 1 all kmc ${temp} ${r_type} ${p_type} ${e_type} ${pexp} ${pot} ${seed} ${nevery}


```
## Theory
Please refer to:
Gadea, E. D., Perez Sirkin, Y. A., Molinero, V., & Scherlis, D. A. (2020). Electrochemically generated nanobubbles: invariance of the current with respect to electrode size and potential. The Journal of Physical Chemistry Letters, 11(16), 6573-6579.
[https://doi.org/10.1021/acs.jpclett.0c01404](https://doi.org/10.1021/acs.jpclett.0c01404)
## Nevery usage
for a simulation at 298K with a preexponential factor of 6e-7
the value of nevery*timestep that is safe to use depends on the potential
Here are some reference values for common potentials

| Potential | Nevery*timestep |
|-----------|-----------------|
| 0.1       | 3514            |
| 0.2       | 501             |
| 0.3       | 71              |
| 0.4       | 10              |
| 0.5       | 1.5             |

For higher potentials there might be necessary to decrease the timestep

## Region optional
The reaction can be restricted to a declared region using the region ${regionID} optional argument in the fix

## Density optional
This keyword controls the density of product particles in a region. When the density is greater than ${dmin} every forward reaction produces the back reaction in the region ${regionID}. If the density is greater than ${dmax} every event produces two back reactions.
```python
fix 1 all kmc ${temp} ${r_type} ${p_type} ${e_type} ${pexp} ${pot} ${seed} ${nevery} density ${regionID} ${dmin} ${dmax}
```

## Potential units
The units of the electrode potential have to be compatible with the units used by LAMMPS (energy/charge).

```python
#potential in volts - metal units - si units
variable volts_pot equal 0.2

#real units
variable potential equal ${volts_pot}*6.022e23/2.611e22

#cgs units
variable potential equal ${volts_pot}*1e-7/3.356e-10

#electron units
variable potential equal ${volts_pot}*27.2114

```
