import numpy as np
import matplotlib.pyplot as plt
import cahnhilliard as ch

chparams          = ch.CHparams();
chparams.m        = 1.0;
chparams.gam      = pow( 0.01 ,2 );
chparams.b        = 1.0;
chparams.u        = 1.0;
chparams.alpha    = 10.0;
chparams.phi_star = 0.0;
chparams.sigma    = 0.0;
nx          = 128;
dx          = 1./nx;
checkpoint  = 20;
maxsteps    = 200;

ch.run_ch_solver(chparams, nx, dx, checkpoint, maxsteps);
