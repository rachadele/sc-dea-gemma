## Usage

```
nextflow main.nf -level class
```

Will extract DEA results for the default experiments and factors at the `class` level defined in `nextflow.config`. To provide your own list, create a `params.json` file like so:

```
{
  "experiments": ["GSE1","GSE2"]
  "experiment_factors": [
    "GSE1": ["factor1, factor2"],
    "GSE2": ["factor1"]
  ]
}
```
and run:
```
nextflow main.nf -params-file params.json -level class
```
