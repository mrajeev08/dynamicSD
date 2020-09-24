# Modeling rabies dynamics in the Serengeti

## Analysis Plan

1. Profile model & speed up: use native formats as much as possible (lists & vectors); avoid filtering but if you do have to do it do it in data.table context with keys? Separate it out as an R pkg undergoing dev. Follow a POMP style with a reporting model. 

2. Test summary statistics and figure out what can diff between incs & local {spatiotemporal autocorrelation of cases} and high het vs. maintained chains of transmission 

3. Do rejection sampling and interactively explore distance metrics and what is sufficient with simulations

4. Use these as priors for pMCMC (SMC) (to get a smoother posterior to sample from)

5. Simulate the best fit.

What's important:
- space in terms of mixing (endemic)
- space in terms of vacc (SD vacc)
- het in transmission 
- dispersal?
- het in timing?
- density?
- landscape of movement?
 
 6. Vacc and connectivity: a better threshold
 
 7. Vacc and connectivity: interventions
 
 8. Scale & models: approximations 

Don't need POMP. Integrating better with Rebecca's work?
