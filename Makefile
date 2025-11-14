# Makefile for running TriCritical-Dicke phase diagram scan

# Number of Julia workers  
NWORKERS=4

# Parameter ranges
e_min=0.0
e_step=0.025
e_max=1.25

g_min=0.8
g_step=0.025
g_max=1.75

# Additional physical parameters
NSpin=100
NBoson=100
ww0=1.0
h=1.0

JULIA= julia
ENTRY= main_TC-Dicke.jl

all:
	E_MIN=$(e_min) E_STEP=$(e_step) E_MAX=$(e_max) \
	G_MIN=$(g_min) G_STEP=$(g_step) G_MAX=$(g_max) \
	N_SPIN=$(NSpin) N_BOSON=$(NBoson) W0=$(ww0) HFIELD=$(h) \
	NWORKERS=$(NWORKERS) \
	$(JULIA) $(ENTRY)

