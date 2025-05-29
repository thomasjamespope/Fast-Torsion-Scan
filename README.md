# Fast_Torsion_Scan
Generate a set of conformers along a torsional degree of freedom in a molecule

Usage:

python torsion-scan.py molecule.xyz i j k l

where:
- molecule.xyz is the moecule in standard xyz format
- i, j, k and l are indexes of the atoms defining the torsion (indices are 0-based)

The easiest way to install Fast_Torsion_Scan is with conda:

```
   git clone https://github.com/thomasjamespope/Fast_Torsion_Scan.git
   cd Fast_Torsion_Scan
   mamba env create -n Fast_Torsion_Scan --file enviroment.yaml
   mamba activate Fast_Torsion_Scan
```
  
