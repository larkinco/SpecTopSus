################################################
########### 	Initialisations		########
################################################	
EnsembleName='A100_24'
Volume=48*24*24*24
Kappa=0.1632550
Mu=0.010
TwoKappaMu=2*Mu*Kappa
#Kappa=TwoKappaMu/(2*0.0080)
Zp_Zs=0.699
r0_a=5.231
a=0.0863 #fm
Mu_R=1000000000 #MeV#Don't know it, this is shouldn't ever come into the calculation, its been made very large so if it used it will be noticed.
Zp_Zs_error=0.013
r0_a_error=0.038
Zp=0.529
################################################
########### 	Initialisations		########
################################################	
Eig_Number=95
average_mode_number=0.0
mode_numbers=[]
temp_mode_num=0.0
Interval=2
StartOfConfs=501

#M_pi=5.35/24
M_pi=0.22293

r0=r0_a*(a)/197.32
GradFlowInterval=2
Bin_Size=1
