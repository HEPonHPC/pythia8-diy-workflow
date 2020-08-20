# pythia8-diy-workflow

## workflow-sample.py
This script will sample a parameter space and write files to run pyhia8-diy
and so forth into an output directory.

### inputs
You need to specify the basic process card, i.e. the pythia8 commands *EXCEPT* the parameters you want to change.
An example is in the repo, called **process.dat**.

A further required input is the set of parameter ranges, there is an example in the repo called **params.ranges**.

### logic
The dimension of the parameter space (**DIM**) is inferred from all *UNCOMMENTED, non-empty* lines in **params.ranges**.
A hashtag at the beginning of the line is interpreted as comment.

The *minimial* number, **NCOEFF**, of sampled points depends on **DIM** and the requested polynomial orders (CL option **--order**).
**NCOEFF** represents the total number of coefficients the approximations will have.

The *actual* number of sampled points, **NPOINTS**, is further influenced by the oversampling factor, **FOVERSAMPLE**,
which is another CL option (**--fidelity**).
I.e. **NPOINTS = NCOEFF * FOVERSAMPLE**.

### Examples

Use this to see the defaults:

python3 workflow-sample.py -h

python3 workflow-sample.py process.dat params.ranges --sampler sobol --order 3,0

python3 workflow-sample.py process.dat params.ranges --sampler lhs --order 3,0 -o scan2 --fidelity 7.3

### Dependencies
* apprentice
* pyDOE2  for lhs sampling
* sobol for sobol sampling
