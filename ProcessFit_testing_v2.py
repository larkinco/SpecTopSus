import numpy as np
import top_class as tc
import data_functions as df

import importlib
import Cut_Off as ct
import fit as ft

M_sqr_ev_set=[ct.M_sqr_ev]

DEnsembleNameList=['D15_48','D20_48','D30_48','D45_32']
#DEnsembleNameList=['D15_48','D20_48','D30_48']
BEnsembleNameList=['B25_32','B35_32','B55_32']
AEnsembleNameList=['A30_32','A40_24','A60_24']

Full_EnsembleNameList=[AEnsembleNameList,BEnsembleNameList,DEnsembleNameList]

#EnsembleNameList=['D15_48','D20_48','D30_48','D45_32']
#EnsembleNameList=['D15_48','D20_48','D30_48','D45_32']
#EnsembleNameList=['B55_32']
EnsembleList=[]
AEnsembleList=[]
BEnsembleList=[]
DEnsembleList=[]
Full_Conf_Dir='/Users/conorlarkin/Full_Conf_Dir2/'
#EigProd_List=[]
#Ensemble_Tuple_list=[]

#from collections import namedtuple
#JackData = namedtuple('JackData','Spacing MassPoint JackList')
#data = JackData(2.5, 1.5,[1,2,2.4])

#Ensemble_r0_Mu=[]
AEnsemble_r0_M_pi_sqr=[]
BEnsemble_r0_M_pi_sqr=[]
DEnsemble_r0_M_pi_sqr=[]

for EnsembleName in AEnsembleNameList:
	param_file=EnsembleName+'_params'
	AEnsembleList.append(importlib.import_module(param_file))
for EnsembleName in BEnsembleNameList:
	param_file=EnsembleName+'_params'
	BEnsembleList.append(importlib.import_module(param_file))
for EnsembleName in DEnsembleNameList:
	param_file=EnsembleName+'_params'
	DEnsembleList.append(importlib.import_module(param_file))
for EnsembleName in DEnsembleNameList:
	param_file=EnsembleName+'_params'
	EnsembleList.append(importlib.import_module(param_file))

#####################################################
##########	OPEN AND READ FILE 	#############
#####################################################


EigProd_List=[]
AEnsem_Top_List=[]
BEnsem_Top_List=[]
DEnsem_Top_List=[]
Ensemble_Tuple_list=[]
AEnsemble_SpecSum_Weight=[]
AEnsemble_SpecProj_Weight=[]
BEnsemble_SpecSum_Weight=[]
BEnsemble_SpecProj_Weight=[]
DEnsemble_SpecSum_Weight=[]
DEnsemble_SpecProj_Weight=[]
AM_sqr=[]
BM_sqr=[]
DM_sqr=[]

a_r0_sqr=[]
Spacing=[]

ChiralContFit=ft.Fit("SpecSum Global Fit")
ChiralContFit_SP=ft.Fit("SpecProj Global Fit")

FitA=ft.Fit("SpecSum A Ensemble")
FitB=ft.Fit("SpecSum B Ensemble")
FitD=ft.Fit("SpecSum D Ensemble")
FitA_SP=ft.Fit("SpecProj A Ensemble")
FitB_SP=ft.Fit("SpecProj B Ensemble")
FitD_SP=ft.Fit("SpecProj D Ensemble")

for i in range(0,len(AEnsembleList),1):
	path=Full_Conf_Dir+AEnsembleList[i].EnsembleName
	EigProd_List =df.Data_Loader(path)
	AM_sqr.append(M_sqr_ev_set[0]*(pow(2*AEnsembleList[i].Kappa*AEnsembleList[i].a/197.32,2)))
	#AM_sqr.append(0.0000567*(pow(2*AEnsembleList[i].Kappa,2)))
	#A_M_sqr_ev_set=[0.0000567*pow(AEnsembleList[i].a/197.32,-2)]
	#A_M_sqr_ev_set=[0.0000567]

	AEnsem_Top_List.append(tc.SpectralData(EigProd_List,AEnsembleNameList[i],AM_sqr[i]))
	AEnsem_Top_List[i].print_results(M_sqr_ev_set)
	AEnsem_Top_List[i].Q_plots(AEnsembleList[i].EnsembleName)
	AEnsemble_r0_M_pi_sqr.append(AEnsem_Top_List[i].r0_M_pi_sqr())
	FitA_SP.AddPoint([pow(AEnsembleList[i].r0_a,4)*y for y in AEnsem_Top_List[i].Renorm_SpecProjector_Jack_List()],AEnsembleNameList[i],AEnsem_Top_List[i].r0_M_pi_sqr(),pow(AEnsembleList[i].r0_a,4)*AEnsem_Top_List[i].Renorm_Top_Suscept_SpecProjector(),'A')
	ChiralContFit.AddPoint([pow(AEnsembleList[i].r0_a,4)*y for y in AEnsem_Top_List[i].Renorm_SpecSum_Jack_List()],AEnsembleNameList[i],AEnsem_Top_List[i].r0_M_pi_sqr(),pow(AEnsembleList[i].r0_a,4)*AEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum(),'A')
	ChiralContFit_SP.AddPoint([pow(AEnsembleList[i].r0_a,4)*y for y in AEnsem_Top_List[i].Renorm_SpecProjector_Jack_List()],AEnsembleNameList[i],AEnsem_Top_List[i].r0_M_pi_sqr(),pow(AEnsembleList[i].r0_a,4)*AEnsem_Top_List[i].Renorm_Top_Suscept_SpecProjector(),'A')
	FitA.AddPoint([pow(AEnsembleList[i].r0_a,4)*y for y in AEnsem_Top_List[i].Renorm_SpecSum_Jack_List()],AEnsembleNameList[i],AEnsem_Top_List[i].r0_M_pi_sqr(),pow(AEnsembleList[i].r0_a,4)*AEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum(),'A')
	AEnsemble_SpecSum_Weight.append(pow(pow(AEnsembleList[i].r0_a,4)*AEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum_Error(),-1))
	a_r0_sqr.append(pow(AEnsembleList[i].r0_a,-2))
	Spacing.append('A')
	
	A_system_quad=pow(pow(2*AEnsembleList[i].Zp_Zs_error/pow(AEnsembleList[i].Zp_Zs,-2),2) +pow(4*AEnsembleList[i].r0_a_error/pow(AEnsembleList[i].r0_a,4),2),0.5)
#	FitA.AddSystematicError(AEnsembleList[i].r0_a,A_system_quad)
#	FitA_SP.AddSystematicError(AEnsembleList[i].r0_a,A_system_quad)
	ChiralContFit.AddSystematicError(AEnsembleList[i].r0_a,A_system_quad)
	ChiralContFit_SP.AddSystematicError(AEnsembleList[i].r0_a,A_system_quad)

for i in range(0,len(BEnsembleList),1):
	path=Full_Conf_Dir+BEnsembleList[i].EnsembleName
	EigProd_List =df.Data_Loader(path)
	BM_sqr.append(M_sqr_ev_set[0]*(pow(2*BEnsembleList[i].Kappa*BEnsembleList[i].a/197.32,2)))
	#BM_sqr.append(0.00003555*(pow(2*BEnsembleList[i].Kappa,2)))
	#B_M_sqr_ev_set=[0.00003555*pow(BEnsembleList[i].a/197.32,-2)]
#	B_M_sqr_ev_set=[0.00003555]
	BEnsem_Top_List.append(tc.SpectralData(EigProd_List,BEnsembleNameList[i],BM_sqr[i]))
	BEnsem_Top_List[i].print_results(M_sqr_ev_set)
	BEnsem_Top_List[i].Q_plots(BEnsembleList[i].EnsembleName)
	BEnsemble_r0_M_pi_sqr.append(BEnsem_Top_List[i].r0_M_pi_sqr())
	FitB.AddPoint([pow(BEnsembleList[i].r0_a,4)*y for y in BEnsem_Top_List[i].Renorm_SpecSum_Jack_List()],BEnsembleNameList[i],BEnsem_Top_List[i].r0_M_pi_sqr(),pow(BEnsembleList[i].r0_a,4)*BEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum(),'B')
	ChiralContFit.AddPoint([pow(BEnsembleList[i].r0_a,4)*y for y in BEnsem_Top_List[i].Renorm_SpecSum_Jack_List()],BEnsembleNameList[i],BEnsem_Top_List[i].r0_M_pi_sqr(),pow(BEnsembleList[i].r0_a,4)*BEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum(),'B')
	ChiralContFit_SP.AddPoint([pow(BEnsembleList[i].r0_a,4)*y for y in BEnsem_Top_List[i].Renorm_SpecProjector_Jack_List()],BEnsembleNameList[i],BEnsem_Top_List[i].r0_M_pi_sqr(),pow(BEnsembleList[i].r0_a,4)*BEnsem_Top_List[i].Renorm_Top_Suscept_SpecProjector(),'B')
	FitB_SP.AddPoint([pow(BEnsembleList[i].r0_a,4)*y for y in BEnsem_Top_List[i].Renorm_SpecProjector_Jack_List()],BEnsembleNameList[i],BEnsem_Top_List[i].r0_M_pi_sqr(),pow(BEnsembleList[i].r0_a,4)*BEnsem_Top_List[i].Renorm_Top_Suscept_SpecProjector(),'B')
	BEnsemble_SpecSum_Weight.append(pow(pow(BEnsembleList[i].r0_a,4)*BEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum_Error(),-1))
	a_r0_sqr.append(pow(BEnsembleList[i].r0_a,-2))
	Spacing.append('B')

	B_system_quad=pow(pow(2*BEnsembleList[i].Zp_Zs_error/pow(BEnsembleList[i].Zp_Zs,-2),2) +pow(4*BEnsembleList[i].r0_a_error/pow(BEnsembleList[i].r0_a,4),2),0.5)
#	FitB.AddSystematicError(BEnsembleList[i].r0_a,B_system_quad)
	#FitB_SP.AddSystematicError(BEnsembleList[i].r0_a,B_system_quad)
	ChiralContFit.AddSystematicError(BEnsembleList[i].r0_a,B_system_quad)
	ChiralContFit_SP.AddSystematicError(BEnsembleList[i].r0_a,B_system_quad)
	

for i in range(0,len(DEnsembleList),1):
	path=Full_Conf_Dir+DEnsembleList[i].EnsembleName
	EigProd_List =df.Data_Loader(path)
	DM_sqr.append(M_sqr_ev_set[0]*(pow(2*DEnsembleList[i].Kappa*DEnsembleList[i].a/197.32,2)))
	#DM_sqr.append(0.000021*(pow(2*DEnsembleList[i].Kappa,2)))
#	D_M_sqr_ev_set=[0.000021*pow(DEnsembleList[i].a/197.32,-2)]
	#D_M_sqr_ev_set=[0.000021]
	DEnsem_Top_List.append(tc.SpectralData(EigProd_List,DEnsembleNameList[i],DM_sqr[i]))
	DEnsem_Top_List[i].print_results(M_sqr_ev_set)
	DEnsem_Top_List[i].Q_plots(DEnsembleList[i].EnsembleName)
	DEnsemble_r0_M_pi_sqr.append(DEnsem_Top_List[i].r0_M_pi_sqr())
	FitD.AddPoint([pow(DEnsembleList[i].r0_a,4)*y for y in DEnsem_Top_List[i].Renorm_SpecSum_Jack_List()],DEnsembleNameList[i],DEnsem_Top_List[i].r0_M_pi_sqr(),pow(DEnsembleList[i].r0_a,4)*DEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum(),'D')
	ChiralContFit.AddPoint([pow(DEnsembleList[i].r0_a,4)*y for y in DEnsem_Top_List[i].Renorm_SpecSum_Jack_List()],DEnsembleNameList[i],DEnsem_Top_List[i].r0_M_pi_sqr(),pow(DEnsembleList[i].r0_a,4)*DEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum(),'D')
	ChiralContFit_SP.AddPoint([pow(DEnsembleList[i].r0_a,4)*y for y in DEnsem_Top_List[i].Renorm_SpecProjector_Jack_List()],DEnsembleNameList[i],DEnsem_Top_List[i].r0_M_pi_sqr(),pow(DEnsembleList[i].r0_a,4)*DEnsem_Top_List[i].Renorm_Top_Suscept_SpecProjector(),'D')
	FitD_SP.AddPoint([pow(DEnsembleList[i].r0_a,4)*y for y in DEnsem_Top_List[i].Renorm_SpecProjector_Jack_List()],DEnsembleNameList[i],DEnsem_Top_List[i].r0_M_pi_sqr(),pow(DEnsembleList[i].r0_a,4)*DEnsem_Top_List[i].Renorm_Top_Suscept_SpecProjector(),'D')
	DEnsemble_SpecSum_Weight.append(pow(pow(DEnsembleList[i].r0_a,4)*DEnsem_Top_List[i].Renorm_Top_Suscept_SpecSum_Error(),-1))
	a_r0_sqr.append(pow(DEnsembleList[i].r0_a,-2))
	Spacing.append('D')
	
	D_system_quad=pow(pow(2*DEnsembleList[i].Zp_Zs_error/pow(DEnsembleList[i].Zp_Zs,-2),2) +pow(4*DEnsembleList[i].r0_a_error/pow(DEnsembleList[i].r0_a,4),2),0.5)
#	FitD.AddSystematicError(DEnsembleList[i].r0_a,D_system_quad)
#	FitD_SP.AddSystematicError(DEnsembleList[i].r0_a,D_system_quad)
	ChiralContFit.AddSystematicError(DEnsembleList[i].r0_a,D_system_quad)
	ChiralContFit_SP.AddSystematicError(DEnsembleList[i].r0_a,D_system_quad)
	
#D_Weight=FitD.ComputeEnsembleErrors()
#print([a - b for a, b in zip(Ensemble_SpecSum_Weight, D_Weight)])
#print(Ensemble_SpecSum_Weight-D_Weight)

#print(Jackknife_Mean_List)

#print("The coeff is ")
#print("The coefficients of the SpecSum fit are, slope and intercept respectively,")
#print(SpecSum_jack_coeff_mean)

#print("The standard error of the SpecSum fit is, slope and intercept respectively,  ")
#print(SpecSum_std_err_coeff)
#print("The coefficients of the SpecProj fit are, slope and intercept respectively,")
#print(SpecProj_jack_coeff_mean)

#print("The standard error of the SpecProj fit is, slope and intercept respectively,  ")
#print(SpecProj_std_err_coeff)

#with open("CoeffSpecSumDist.txt", 'w') as f:
#	for coeff in SpecSum_coeff_list:
#		f.write("%s %s \n" % (coeff[0], coeff[1]))

#with open("CoeffSpecProjDist.txt", 'w') as f:
#	for coeff in SpecProj_coeff_list:
#		f.write("%s %s \n" % (coeff[0], coeff[1]))

#Spacing=EnsembleNameList[0][0]
#filepath='Output/FitResults'+Spacing+'.txt'

#f=open(filepath,"w")
#output=Spacing +','+str(SpecSum_jack_coeff_mean[0]) +','+str(SpecSum_jack_coeff_mean[1]) +','+str(SpecSum_std_err_coeff[0])+','+str(SpecSum_std_err_coeff[1])+','+str(SpecProj_jack_coeff_mean[0]) +','+str(SpecProj_jack_coeff_mean[1]) +','+str(SpecProj_std_err_coeff[0])+','+str(SpecProj_std_err_coeff[1])+','+str(Ensem_Top_List[0].a)+','+str(SpecSum_phys_mass_mean)+','+str(SpecProj_phys_mass_mean)+','+str(SpecSum_phys_mass_std_dev)+','+str(SpecProj_phys_mass_std_dev)+','+str(ct.r0_physical_mass)+','+str(pow(Ensem_Top_List[0].r0_a,-1))+'\n'

#f.write(output)

#f.close()


print("fit stuff")

FitA_SP.fit_data()
FitB_SP.fit_data()
FitD_SP.fit_data()

FitA.fit_data()
FitB.fit_data()
FitD.fit_data()

ContFits=[]
#print(FitB.EnsembleAvgList)
for M in ct.MassPoints_sqr:
	MassPointJackLists=[]
#	print(FitB.MassPointJackList(M))
	filepath='Output/FixedMassFitResults'+str(M)+'_SpecSum.txt'
	FitContM=ft.Fit("Cont_FixedMass_fits_SpecSum "+str(M)+" ")	
	FitContM.AddPoint(FitA.MassPointJackList(M),str(M)+'_A',pow(0.1911680,2),FitA.MassPointPredwError(M)[0],str(M))
#	FitContM.AddSystematicError(AEnsembleList[0].r0_a,A_system_quad)
	print("the A fit Mass point and error is")
	print(FitA.MassPointPredwError(M))
	FitContM.AddPoint(FitB.MassPointJackList(M),str(M)+'_B',pow(0.1751313,2),FitB.MassPointPredwError(M)[0],str(M))
#	FitContM.AddSystematicError(BEnsembleList[0].r0_a,B_system_quad)
	print("the B fit Mass point and error is")
	print(FitB.MassPointPredwError(M))
	FitContM.AddPoint(FitD.MassPointJackList(M),str(M)+'_D',pow(0.1326611,2),FitD.MassPointPredwError(M)[0],str(M))
#	FitContM.AddSystematicError(DEnsembleList[0].r0_a,D_system_quad)
	print("the D fit Mass point and error is")
	print(FitD.MassPointPredwError(M))
	FitContM.fit_data()
	FitContM.SaveResults(filepath)

for M in ct.MassPoints_sqr:
	MassPointJackLists=[]
#	print(FitB.MassPointJackList(M))
	filepath='Output/SPFixedMassFitResults'+str(M)+'_SpecSum.txt'
	FitContM_SP=ft.Fit("Cont_FixedMass_fits_SpecProj "+str(M)+" ")	
	FitContM_SP.AddPoint(FitA_SP.MassPointJackList(M),str(M)+'_A',pow(0.1911680,2),FitA_SP.MassPointPredwError(M)[0],str(M))
	#FitContM_SP.AddSystematicError(AEnsembleList[0].r0_a,A_system_quad)
	FitContM_SP.AddPoint(FitB_SP.MassPointJackList(M),str(M)+'_B',pow(0.1751313,2),FitB_SP.MassPointPredwError(M)[0],str(M))
	#FitContM_SP.AddSystematicError(BEnsembleList[0].r0_a,B_system_quad)
	FitContM_SP.AddPoint(FitD_SP.MassPointJackList(M),str(M)+'_D',pow(0.1326611,2),FitD_SP.MassPointPredwError(M)[0],str(M))
	#FitContM_SP.AddSystematicError(DEnsembleList[0].r0_a,D_system_quad)
	FitContM_SP.fit_data()
	FitContM_SP.SaveResults(filepath)

ChiralContFit.FitChiralCont_1Slope(a_r0_sqr,Spacing)
#ChiralContFit.JackknifedSubtraction(a_r0_sqr,'Output/SpecSumContChiralSubtraction1Slope.txt')
ChiralContFit_SP.FitChiralCont_1Slope(a_r0_sqr,Spacing)
#ChiralContFit_SP.JackknifedSubtraction(a_r0_sqr,'Output/SpecProjContChiralSubtraction1Slope.txt')
ChiralContFit.FitChiralCont(a_r0_sqr,Spacing)
#ChiralContFit.JackknifedSubtraction(a_r0_sqr,'Output/SpecSumContChiralSubtraction3slope.txt')
ChiralContFit_SP.FitChiralCont(a_r0_sqr,Spacing)
#ChiralContFit_SP.JackknifedSubtraction(a_r0_sqr,'Output/SpecProjContChiralSubtraction3slope.txt')
#ChiralContFit.AddSystematicError(EnsembleList[i].r0_a,D_system_quad)
ChiralContFit.FitChiralCont_1Slope_Cross_wSubtraction(a_r0_sqr,Spacing,'Output/SubtractionSSResults.txt')
ChiralContFit_SP.FitChiralCont_1Slope_Cross_wSubtraction(a_r0_sqr,Spacing,'Output/SubtractionSPResults.txt')


ChiralContFit.FitChiralCont_1Slope_Cross_wSubtraction_NoIntercept(a_r0_sqr,Spacing,'Output/SubtractionSSResults_NoIntercept.txt')
ChiralContFit_SP.FitChiralCont_1Slope_Cross_wSubtraction_NoIntercept(a_r0_sqr,Spacing,'Output/SubtractionSPResults_NoIntercept.txt')

