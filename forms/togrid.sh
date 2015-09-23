# Create plots
cd hyperelasticity; python ~/Projects/Firedrake/firedrake-bench/forms/togrid_q_nf_mode.py . hyperelasticity; cd ..
cd elasticity; python ~/Projects/Firedrake/firedrake-bench/forms/togrid_q_nf_mode.py . elasticity; cd ..
cd helmholtz; python ~/Projects/Firedrake/firedrake-bench/forms/togrid_q_nf_mode.py . helmholtz; cd ..
cd mass; python ~/Projects/Firedrake/firedrake-bench/forms/togrid_q_nf_mode.py . mass; cd ..
# Copy them
cp hyperelasticity/hyperelasticity.pdf /homes/fl1612/Desktop/Dropbox/MyPapers/OptimalFEMIntegration/TOMS_Paper/perf-results/
cp elasticity/elasticity.pdf /homes/fl1612/Desktop/Dropbox/MyPapers/OptimalFEMIntegration/TOMS_Paper/perf-results/
cp helmholtz/helmholtz.pdf /homes/fl1612/Desktop/Dropbox/MyPapers/OptimalFEMIntegration/TOMS_Paper/perf-results/
cp mass/mass.pdf /homes/fl1612/Desktop/Dropbox/MyPapers/OptimalFEMIntegration/TOMS_Paper/perf-results/
