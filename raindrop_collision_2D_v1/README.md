This is 2D case to simulate collision outcome.

To run this case, type "./external" in terminal window.

External is executable file, which controls the couple between standard solver and external codes

At the preprocessing, the code merge two simulated free falling result together to set up initial conditions

The external file do the following time loops: 1 run interFoam for a short time 2 check the center of raindrop and move properties accordingly
