# MDBundle
A collection of independent scripts for MD simulations.


## Contents

### WATER_3to4
Convert a system with a 3-point water model into a new one with a 4-point model.

```bash
python WATER_3to4.py -f PATH/TO/SYSTEM -o PATH/TO/NEW_SYSTEM -wr WATER_RESNAME -nr NEW_RESNAME -nn NEW_O NEW_H1 NEW_H2 NEW_VS -vd DISTANCE_O_VS
```
