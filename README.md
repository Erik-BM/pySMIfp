# pySMIfp

Create SMILES fingerprints. 

Python version of [this](https://pubs.acs.org/doi/10.1021/ci400206h).

## Use
```
>>> from smiles_fingerprints import smiles_fingerprint
>>> smiles = 'C1=CC=C(C=C1)OP(=O)(OC2=CC=CC=C2)OC3=CC=CC=C3' # triphenyl phosphate
>>> smiles_fingerprint(smiles)
>>> [3, 4, 10, 0, 0, 18, 0, 1.0, 1.0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
```
