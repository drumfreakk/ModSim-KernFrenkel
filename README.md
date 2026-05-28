# ModSim-KernFrenkel

*Project 3) Crystals of patchy particles*

*Goal: Simulate trivalent colloidal particles and determine their crystaline phases*

-----------------------------------------------------
***Current To Do list*** (**crystaline phases**)
- Decide upon equilibration time or a observable that would indicate EQ time
- Create order parameter analysis code that tells us what phase/unit cell we have, then graph that

-----------------------------------------------------
***Current Achivements***: (**sucesfully simulated**)
- Corretly simulate trihedral or tetrahedral KF patchy particles
- Output coords file that can be used with the CVT to visalise the simulation
- Output energy and snapshot files
- Created MakeFile config that can run multiple set of lambda_cos or temp_density parameters consecutively and output results
-----------------------------------------------------



Plan
- Make a script that analyses radius autocorr, bond order parameter
- Make a script that plot all configurations
- Run simulation with tetrahedral particles and compare to data, then rerun with trihedral
- Potentiall explore the possibilities in the HowTo paper

-----------------------------------------------------
Diamond: 
- each particle is bonded with exactly one other particle

BCC
- each particle has 6 nearest neighbours, 3 with bonds 3 without

*Possibilities:*
- Bond order parameter analysis
- Repurpose radial distribution code
-----------------------------------------------------

Mateos CVT command:
*OUT file:* ./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/coords___.ptc
*Temp_density:* ./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.02_rho_0.3/coords___.ptc
*Lambda_cos:*

-----------------------------------------------------

./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.1_rho_0.1/coords___.ptc

./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.1_rho_0.5/coords___.ptc

./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.1_rho_0.72/coords___.ptc


./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.2_rho_0.1/coords___.ptc

./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.2_rho_0.5/coords___.ptc

./Visualisation/CVT/build-mac-release/cvt ModSim-KernFrenkel/src/out/T_0.2_rho_0.72/coords___.ptc

output_4patch_256_s1e6/*
