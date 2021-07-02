# MESA_ALPs
Edited version of stellar evolution code MESA including energy-loss to ALPs. This code was used in the paper 2102.00379.

This repository contains three files:
  1. bsm.txt - file specifying ALP coupling strength and mass
  2. ALP_Data.txt - file containing pre-processed grid of (F(mu, xi)+G(mu)) values as defined in the paper.
  3. run_star_extras.f - a modified run_star_extras file containing additional routines required for computing energy loss to ALP production by interpolating between the pre-processed grid. These are based in code from 1210.1271.

To use these routines in MESA, add use_other_neu = .true. to your &controls namelist.

This was made using MESA v9691.
