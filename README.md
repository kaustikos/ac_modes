# ac_modes
simple MATLAB code for the computation of acoustical normal modes in the ocean

Test cases

Test case 1. Modes for the ASA wedge benchmark

In this example we consider the computation of wavenumbers and mode functions for the well-known ASA wedge benchmark. Sound frequency f = 25 Hz, we have two-layer shallow-water environment with the water column parameters cw = 1500 m/s and ρw = 1 g/cm3. The bottom is a homogeneous liquid sediment with sound speed cb = 1700 m/s, density ρw = 1.5 g/cm3 and attenuation βb = 0.5 dB/λ. Water depth takes the values from 34 m to 380 m (with the step 2 m). Complex wavenumbers (taking attenuation into account) are computed for all values of depth, as well as the values of mode functions at zs = 100 m and zr = 30 m. The computation can be accomplished by running the test script Test_case_1_ASA_wedge_25Hz.m.
	We also compute mode group velocities for each value of water depth. Moreover, the real part of the wavenumber is substituted into the Pekeris dispersion relation, and the resulting error is estimated. 
	Note that 8 modes in total are computed for each water depth. Naturally, when the cut-off depth of certain mode is reached, it turns into the “bottom” mode (approximating some mode of continuous spectrum). 
	In this computation we use 3 grids for Richardson extrapolation with the uniform meshsizes in depths of 0.25, 0.5 and 1 m. We also use false bottom at z = 1000 m, where a Dirichlet boundary condition is set up. 
	Output: text files ‘kj_wedge_att.txt’ (complex-valued wavenumbers), ‘phizr_wedge.txt’ (mode functions at z = zr = 30 m), ‘phizs_wedge.txt’ (mode functions at z = zs = 100 m), ‘vgr.txt’ (mode group velocities – make sense only for trapped modes!), ‘err_pek.txt’ (residual of the Pekeris dispersion relation as computed for waterborne modes).

Test case 2. Shallow-water waveguide with tanh(z) SSP

	In the second example a shallow-water waveguide with non-constant sound speed in the water column is considered. The profile is given by the formula
c(z)=c_0-∆c tanh⁡〖((z-z_0)/σ),〗
where c_0=1490 m/s, z_0=25 m, ∆c=30 m/s, σ=10 m. We assume that h0 = 50 m, and consider the sound frequencies of f = 100, 200 and 400 Hz. We also set cb = 2000 m/s, density ρb = 3 g/cm3. For each frequency we compute wavenumbers and modal functions for the false bottom with the Dirichlet condition at z = 200 m. We use the same grids for doing Richardson extrapolation as in the previous case. In this case 20 modes are computed. In the MATLAB code the step ∆z=0.125 m is used.
	Output: text files ‘kj_f_***Hz.txt’ (complex-valued wavenumbers for the given frequency f = ***), ‘phij_f_***Hz.txt’ (mode functions), ‘vgr.txt’ (mode group velocities – make sense only for trapped modes!).

Test case 3. Deep-water scenario, Munk’s SSP

	In the third example we consider the computation of modes for the canonic Munk’s sound speed profile in the deep ocean. The total depth is 3 km, and the Dirichlet boundary condition is set up both at the bottom and at the surface. The SSP in the water column is described by the formula
c(z)=c_0 (1+ε(z ̅-1+exp⁡〖(-z ̅)〗 )),
where c_0=1500 m/s, ε=0.00737, z ̅=2(z-1300)/1300. First 30 modes are computed for the sound frequencies of f = 100, 200 and 400 Hz.
Output: text files ‘kj_f_***Hz.txt’ (real-valued wavenumbers for the given frequency f = ***), ‘phij_f_***Hz.txt’ (mode functions), ‘vgr.txt’ (mode group velocities – make sense only for trapped modes!).




