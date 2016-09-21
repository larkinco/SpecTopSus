import math
import os
import numpy as np
#import jackknife as jk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import importlib
import Cut_Off as ct

class EigenProd:
	def __init__(self, evls, product):
		self.evls = evls
		self.product = product
class SpectralData:
	def __init__(self, EigenProd_List,EnsembleName,M_sqr_ev):
		param_file=EnsembleName+'_params'
		pm=importlib.import_module(param_file)
		
		self.Q_list=[]
		self.spectral_sum_1_list=[]
		self.spectral_sum_2_list=[]
		self.full_eigenvalue_collection=[] #List of all eigenvalues, used to look at spectral density
		self.TwoKappaMu=pm.TwoKappaMu
		self.EnsembleName=pm.EnsembleName
		self.Lat_Volume=pm.Volume
		self.TwoKappaMu=pm.TwoKappaMu
		self.Kappa=pm.Kappa
		self.Zs_Zp_ratio=pow(pm.Zp_Zs,-1)
		self.r0_a=pm.r0_a
		self.a=pm.a#fm
		self.Mu_R=pm.Mu_R #MeV
		self.Zs_Zp_error=pm.Zp_Zs_error
		self.r0_a_error=pm.r0_a_error
		self.r0=pm.r0
		self.Zp=pm.Zp
		self.M_sqr=pow(self.Zp,2)*M_sqr_ev*pow(2*self.Kappa*self.a/197.32,2)
		print("This is M_sqr")
		print(self.M_sqr)
		self.M_pi=pm.M_pi	
		self.Interval=pm.Interval
		self.StartOfConfs=pm.StartOfConfs	
		self.r0Mu_R=self.r0*self.Mu_R		
		self.ModeList=[]
		self.BinSize=pm.Bin_Size
		product_sum=0.0
		spectral_sum1=0.0
		spectral_sum2=0.0
		self.EigenProd_List=EigenProd_List
		
		for EigenProd in EigenProd_List:
			counter=0
			ModeNumber=0

		#	if __debug__:
		#		print(pow(TwoKappaMu,2)/EigenProd.evls[0])
			if EigenProd.evls[-1]<self.M_sqr:
				print("M_sqr TOO LARGE \n")
			for j in range(0,len(EigenProd.evls),1):
				counter+=1
				spectral_sum1+=pow(self.TwoKappaMu,4)*EigenProd.product[j]*pow(EigenProd.evls[j],-2)####CHANGE
				spectral_sum2+=pow(self.TwoKappaMu,4)*EigenProd.product[j]*pow(EigenProd.evls[j],-2)
				self.full_eigenvalue_collection.append(EigenProd.evls[j]) ##Maybe scale by TwoKappaMu or at least by TwoKappa
				if EigenProd.evls[j] < self.M_sqr:
					#counter=counter+1
					product_sum+=EigenProd.product[j]
					ModeNumber+=1
	#               		cutoff_spectral_sum1+=pow(EigenProd.evls[j],-1)*EigenProd.product[j]
	#               		cutoff_spectral_sum2+=pow(EigenProd.evls[j],-2)*EigenProd.product[j]
######## ADD CHECK that len(EigenProd.evls)==number of eigenvalues expected, 
			if __debug__:
				print(counter)
				print(pow(self.TwoKappaMu,2)/EigenProd.evls[0])	
			if (pow(self.TwoKappaMu,2)/EigenProd.evls[0])<0.90:
				print("A lowest lying eigenvalue is quite large, this may not be a problem but check data is correct")
			
			self.Q_list.append(product_sum)
			self.spectral_sum_1_list.append(spectral_sum1)
			self.spectral_sum_2_list.append(spectral_sum2)	
			
			spectral_sum1=0.0
			spectral_sum2=0.0
			product_sum=0.0
			self.ModeList.append(ModeNumber)
		self.NumOfConfs = len(EigenProd_List)
		self.NumOfBins=math.ceil(self.NumOfConfs/self.BinSize)
		print("The number of configurations is ",len(EigenProd_List))

	def r0_M_pi_sqr(self):
		return pow(self.M_pi*self.r0_a,2)
	
	def Jack_Top_Suscept_SpecSum(self,i):

		Q_arr_a=np.array(self.spectral_sum_1_list)
		Q_arr_b=np.array(self.spectral_sum_2_list)
#		Q_sqr_sum = np.multiply(Q_arr_a,Q_arr_b)			
		Q_sqr = np.multiply(Q_arr_a,Q_arr_b)	
		
		Binned_Q_sqr=self.Binning(Q_sqr)
		
		Binned_Q_sqr_sum=sum(Binned_Q_sqr)

		NumOfBins=self.NumOfBins
		if i >-1:
			Binned_Q_sqr_sum=Binned_Q_sqr_sum - Binned_Q_sqr[i]
			NumOfBins+=-1
		return Binned_Q_sqr_sum/(self.Lat_Volume*NumOfBins)

	def Binning(self,Q_sqr):
		Q_sqr_list=list(Q_sqr)
		Binned_Q_sqr=[]
		i=0
		Q_sqr_Bin=[]
		for Q_sqr in Q_sqr_list:
			Q_sqr_Bin.append(Q_sqr)
			i+=1
			if((i%self.BinSize)==0):
				mean = sum(Q_sqr_Bin)/len(Q_sqr_Bin)
				Binned_Q_sqr.append(mean)
				Q_sqr_Bin=[]
		
		if Q_sqr_Bin:
			i+=1
			mean = sum(Q_sqr_Bin)/len(Q_sqr_Bin)
			Binned_Q_sqr.append(mean)
			Q_sqr_Bin=[]
		self.NumOfBins=i
		self.NumOfBins=len(Binned_Q_sqr)
		return Binned_Q_sqr

	def Set_M_sqr(self,M_sqr):
		self.M_sqr=M_sqr
		a=2 ####CHANGE
		b=2
		spectral_sum_a=0.0
		spectral_sum_b=0.0
		spectral_sum_a_list=[]
		spectral_sum_b_list=[]	
		Q=0.0
		Q_list=[]
		self.ModeList=[]
		
		for EigenProd in self.EigenProd_List:
			ModeNumber=0
			if EigenProd.evls[-1]<M_sqr:
				print("M_sqr TOO LARGE \n")
		#	print(M_sqr)
			for j in range(0,len(EigenProd.evls),1):
				if EigenProd.evls[j] < self.M_sqr:
					Q+=EigenProd.product[j]
					ModeNumber+=1
				if ct.GLOBALMU:
					latMU=ct.MU*self.Zp*self.a/197.32
					spectral_sum_a+=pow(2*self.Kappa*latMU,2.0*a)*EigenProd.product[j]*pow(EigenProd.evls[j]-self.TwoKappaMu*self.TwoKappaMu +4*self.Kappa*self.Kappa*latMU*latMU,-1.0*a)
					spectral_sum_b+=pow(2*self.Kappa*latMU,2.0*b)*EigenProd.product[j]*pow(EigenProd.evls[j]-self.TwoKappaMu*self.TwoKappaMu +4*self.Kappa*self.Kappa*latMU*latMU,-1.0*b)
				else:
					spectral_sum_a+=pow(self.TwoKappaMu,2.0*a)*EigenProd.product[j]*pow(EigenProd.evls[j],-1.0*a)
					spectral_sum_b+=pow(self.TwoKappaMu,2.0*b)*EigenProd.product[j]*pow(EigenProd.evls[j],-1.0*b)
			Q_list.append(Q)
			Q=0.0
			spectral_sum_a_list.append(spectral_sum_a)
			spectral_sum_b_list.append(spectral_sum_b)		
			spectral_sum_a=0.0
			spectral_sum_b=0.0
			self.ModeList.append(ModeNumber)
			
		self.Q_list=Q_list
		self.spectral_sum_1_list=spectral_sum_a_list
		self.spectral_sum_2_list=spectral_sum_b_list
		self.TopChargeSaveFile()
		self.TopSusSaveFile()
		

	def Jack_Top_Suscept_SpecProjector(self,i):
		Q_arr=np.array(self.Q_list)
		Q_sqr = np.multiply(Q_arr,Q_arr)
		Binned_Q_sqr=self.Binning(Q_sqr)	
		Binned_Q_sqr_sum=sum(Binned_Q_sqr)
		NumOfBins=self.NumOfBins

#		NumOfConfs=self.NumOfConfs
		if i >-1:
			Binned_Q_sqr_sum=Binned_Q_sqr_sum - Binned_Q_sqr[i]
			NumOfBins+=-1
		return Binned_Q_sqr_sum/(self.Lat_Volume*NumOfBins)

	def Top_Suscept_SpecSum(self):
		return Jack_Top_Suscept_SpecSum(-1)
	def Top_Suscept_SpecSum_Error(self):
		Spec_array1 = np.array(self.spectral_sum_1_list)
		Spec_array2 = np.array(self.spectral_sum_2_list)
		SpecSum_array = np.multiply(Spec_array1,Spec_array2)
		Binned_Q_sqr=self.Binning(SpecSum_array)
		return self.jackknife_error_mean(Binned_Q_sqr)/self.Lat_Volume

	def Renorm_Top_Suscept_SpecSum(self):
		return (self.Zs_Zp_ratio*self.Zs_Zp_ratio*self.Jack_Top_Suscept_SpecSum(-1))
	def Jack_Renorm_Top_Suscept_SpecProjector(self,i):
		return (self.Zs_Zp_ratio*self.Zs_Zp_ratio*self.Jack_Top_Suscept_SpecProjector(i))
	def Jack_Renorm_Top_Suscept_SpecSum(self,i):
		return (self.Zs_Zp_ratio*self.Zs_Zp_ratio*self.Jack_Top_Suscept_SpecSum(i))
	def Renorm_Top_Suscept_SpecProjector(self):
		return self.Jack_Renorm_Top_Suscept_SpecProjector(-1)
	def Top_Suscept_SpecProjector(self):
		return self.Jack_Top_Suscept_SpecProjector(-1)	
	
	def Renorm_Top_Suscept_Gaussian_Def_ab(self,a,b):
		spectral_sum_a=0.0
		spectral_sum_b=0.0
		spectral_sum_a_list=[]
		spectral_sum_b_list=[]	
		for EigenProd in self.EigenProd_List:
			for j in range(0,len(EigenProd.evls),1):
				if EigenProd.evls[j] < self.M_sqr:
					spectral_sum_a+=pow(self.TwoKappaMu,2.0*a)*EigenProd.product[j]*pow(EigenProd.evls[j],-1.0*a)
					spectral_sum_b+=pow(self.TwoKappaMu,2.0*b)*EigenProd.product[j]*pow(EigenProd.evls[j],-1.0*b)
			spectral_sum_a_list.append(spectral_sum_a)
			spectral_sum_b_list.append(spectral_sum_b)	
			
			spectral_sum_a=0.0
			spectral_sum_b=0.0

		Q_arr_a=np.array(spectral_sum_a_list)
		Q_arr_b=np.array(spectral_sum_b_list)
		Q_sqr_sum = self.Zs_Zp_ratio*self.Zs_Zp_ratio*np.dot(Q_arr_a,Q_arr_b)

		return Q_sqr_sum/(self.Lat_Volume*self.NumOfConfs)


	def Renorm_Top_Suscept_Gaussian_Def01(self):
		Q_arr0=np.array(self.Renorm_Top_Charge_Def0())
		Q_arr1=np.array(self.Renorm_Top_Charge_Def1())
		Q_sqr_sum = np.dot(Q_arr0,Q_arr1)

		return Q_sqr_sum/(self.Lat_Volume*self.NumOfConfs)
	def Renorm_Top_Suscept_Gaussian_Def0(self):
		Q_arr=np.array(self.Renorm_Top_Charge_Def0())
		std_dev = np.std(Q_arr)
		return std_dev*std_dev/(self.Lat_Volume)
	def Renorm_Top_Suscept_Gaussian_Def1(self):
		Q_arr=np.array(self.Renorm_Top_Charge_Def1())
		std_dev = np.std(Q_arr)
		return std_dev*std_dev/(self.Lat_Volume)
	def Renorm_Top_Suscept_Gaussian_Def2(self):
		Q_arr=np.array(self.Renorm_Top_Charge_Def2())
		std_dev = np.std(Q_arr)
		return std_dev*std_dev/(self.Lat_Volume)
	
	def Renorm_Top_Charge_Def2(self):
		return [self.Zs_Zp_ratio*x for x in self.spectral_sum_2_list]
	
#	def Jack_Top_Suscept_SpecProjector(self,i):
#		Q_array = np.array(self.Q_list)
#		Q_sqr_sum = np.dot(Q_array,Q_array)
#		return (Q_sqr_sum/(self.Lat_Volume*self.NumOfConfs))	

	def Top_Suscept_SpecProjector_Error(self):
		Q_arr=np.array(self.Q_list)
		Q_array_sqr=np.multiply(Q_arr,Q_arr)
		Binned_Q_sqr=self.Binning(Q_array_sqr)
		return self.jackknife_error_mean(Binned_Q_sqr)/self.Lat_Volume
	def Renorm_Top_Suscept_SpecSum_Error(self):
		return (self.Zs_Zp_ratio*self.Zs_Zp_ratio*self.Top_Suscept_SpecSum_Error())
	def Renorm_Top_Suscept_SpecProjector_Error(self):
		return (self.Zs_Zp_ratio*self.Zs_Zp_ratio*self.Top_Suscept_SpecProjector_Error())
#	def Number_Of_Confs(self):
#		return self.NumOfBins
	def Renorm_Top_Charge_Def0(self):
		return [self.Zs_Zp_ratio*x for x in self.Q_list]
		#return (self.Zs_Zp_ratio*self.Q_list)
	def Renorm_Top_Charge_Def1(self):
		return [self.Zs_Zp_ratio*x for x in self.spectral_sum_1_list]
	def Full_Eigenvalue_Collection(self):
		return self.full_eigenvalue_collection
	def Save_Full_Eigenvalue_Collection(self,file_link):
		with open(file_link, 'w') as f:
			for i in range(0,len(self.full_eigenvalue_collection),1):
				f.write("%s\n" % (self.full_eigenvalue_collection[i]))	
	
	def Renorm_SpecProjector_Jack_List(self):
		JackList=[]
		for i in range(0,self.NumOfBins,1):
			JackList.append(self.Jack_Renorm_Top_Suscept_SpecProjector(i))
		return JackList
	def Renorm_SpecSum_Jack_List(self):
		JackList=[]
		for i in range(0,self.NumOfBins,1):
			JackList.append(self.Jack_Renorm_Top_Suscept_SpecSum(i))
		return JackList


	def FilePrint(self,file_link):
                with open(file_link, 'w') as f:
 #                       x=StartOfConfs
                         for i in range(0,len(self.Q_list),1):
                                f.write("%s %s %s \n" % (self.Q_list[i],self.spectral_sum_1_list[i],self.spectral_sum_2_list[i]))	
#				x+=Interval

	def TopChargeSaveFile(self):
		with open('Output/Top_Charge_'+self.EnsembleName+'.txt', 'w') as f:
			x=self.StartOfConfs
			f.write("ConfNumber,Q_0,Q_1,Q_2 \n")	
			for i in range(0,len(self.Q_list),1):
				f.write("%d %s %s %s \n" % (x,self.Zs_Zp_ratio*self.Q_list[i],self.Zs_Zp_ratio*self.spectral_sum_1_list[i],self.Zs_Zp_ratio*self.spectral_sum_2_list[i]))	
				x+=self.Interval
	
	def TopSusSaveFile(self):
		with open('Output/Top_Sus_'+self.EnsembleName+'.txt', 'w') as f:
			x=self.StartOfConfs
			f.write("Topological Susceptibility Contribution Per Conf \n")	
			for i in range(0,len(self.Q_list),1):
				xtop_per_conf=pow(self.r0_a,4)*self.Zs_Zp_ratio*self.Zs_Zp_ratio*self.Q_list[i]*self.Q_list[i]/self.Lat_Volume
				f.write("%s \n" %(xtop_per_conf))

	def jackknife_error_mean(self,data_array):
		n = len(data_array)
	#       print(n)
		data_sum_mean = sum(data_array)/float(len(data_array))
		data_jack_omit=[]
		for j in range(0,len(data_array),1):
			data_jack_omit.append((n*data_sum_mean - data_array[j])/(n-1.0))
		jack_mean = sum(data_jack_omit)/len(data_array)
		error_sum=0.0
		for j in range(0,len(data_array),1):
			error_sum+=(jack_mean-data_jack_omit[j])*(jack_mean-data_jack_omit[j])
		return np.sqrt(error_sum*(n-1.0)/n)
	def jackknife_bin_analysis(self,data_array):
		n = len(data_array)
	#       print(n)
		data_sum_mean = sum(data_array)/float(len(data_array))
		data_jack_omit=[]
		for j in range(0,len(data_array),1):
			data_jack_omit.append((n*data_sum_mean - data_array[j])/(n-1.0))
		jack_mean = sum(data_jack_omit)/len(data_array)
		error_sum=0.0
		for j in range(0,len(data_array),1):
			error_sum+=(jack_mean-data_jack_omit[j])*(jack_mean-data_jack_omit[j])
		return np.sqrt(error_sum*(n-1.0)/n)


	def Q_plots(self,EnsembleName):
		fig = plt.figure()

		ax1 = fig.add_subplot(111)

		ax1.set_title("Topological Charge vs. Monte Carlo Time")
		ax1.set_ylabel('Q')
		ax1.set_xlabel('Configurations')

		x=range(self.NumOfConfs)
		#ax1.plot(r0_Mu_R,SpecSum,".", label='the data',linestyle='--')
		ax1.plot(x,self.Q_list,".-", label='the data',color='red')
	#	ax1.plot(x,self.spectral_sum_1_list,".-", label='the data')
		ax1.plot(x,self.spectral_sum_2_list,".-", label='the data',color='blue')

		plt.savefig("Q_vs_montecarlo_t"+EnsembleName + ".pdf", format='pdf')
	


	def print_results_2(self,M_sqr_list_2):
		filepath="Output/Results"+self.EnsembleName+".txt"
		f = open(filepath,"w") 
		conv_const=pow(self.Zp,2)*pow(2*self.Kappa*self.a/(197.32),2)
		conv_const2=1
		M_sqr_list= [x*conv_const2 for x in M_sqr_list_2]

		for M_sqr in M_sqr_list:
			a=0.5
			b=0.5
			self.Set_M_sqr(M_sqr)
			print(M_sqr)
			spec_proj=self.Renorm_Top_Suscept_SpecProjector()
			#spec_proj=x.Jack_Renorm_Top_Suscept_SpecProjector(-1)
			spec_proj=self.Renorm_Top_Suscept_SpecProjector()
			spec_sum=self.Renorm_Top_Suscept_SpecSum()
			#spec_proj_error=1.0
			spec_proj_error=self.Renorm_Top_Suscept_SpecProjector_Error()
			spec_sum_error=self.Renorm_Top_Suscept_SpecSum_Error()
			QDef1_top_sus=self.Renorm_Top_Suscept_Gaussian_Def1()
			QDef_ab_top_sus=self.Renorm_Top_Suscept_Gaussian_Def_ab(a,b)
			spec_sum_quad=pow( spec_sum_error/spec_sum,2)
			spec_proj_quad=pow( spec_proj_error/spec_proj,2)
			renorm_quad=pow(2*self.Zs_Zp_error/pow(self.Zs_Zp_ratio,2),2)
			#print(renorm_quad)
			r0_a_quad=pow(4*self.r0_a_error/pow(self.r0_a,4),2)
			#print(r0_a_quad)
			spec_sum_total_error=pow(self.r0_a,4)*spec_sum*pow(spec_sum_quad+renorm_quad+r0_a_quad,0.5)
			spec_proj_total_error=pow(self.r0_a,4)*spec_proj*pow(spec_proj_quad+renorm_quad+r0_a_quad,0.5)
			#### CAREFUL WITH THIS
			
			self.Save_Full_Eigenvalue_Collection("Output/Eigenvalue_Dist"+self.EnsembleName+".txt")
			Renorm_M_sqr_ev=M_sqr/conv_const

			print(self.EnsembleName)
			print ("Spectral Values")
			print ('The average mode number is %.5e' % (sum(self.ModeList)/len(self.ModeList)))
			print ('The spectral projector a^4 \chi is: %.5e' % (spec_proj))
			print ('The spectral projector r0^4 \chi is: %.5e' % (pow(self.r0_a,4)*spec_proj))
			print ('The spectral sum a^4 \chi is: %.5e' % (spec_sum))
			print ('The spectral sum r0^4 \chi is: %.5e' % (pow(self.r0_a,4)*spec_sum))
			print ("These are the Jackknife std. deviations")
			print ('The spectral projector r_0^4 \chi error is: %.5E' % (pow(self.r0_a,4)*(spec_proj_error)))
			print ('The spectral projector a^4 \chi error is: %.5E' % ((spec_proj_error)))
			print ('The spectral sum r_0^4 \chi error is: %.5E' %  (pow(self.r0_a,4)*spec_sum_error))
			print ('The total spectral projector r_0^4 \chi error is: %.5E' % (spec_proj_total_error))
			print ('The total spectral sum r_0^4 \chi error is: %.5E' %  (spec_sum_total_error))
			print ('The QDef1 r_0^4 topological suceptibility is: %.5E' % (pow(self.r0_a,4)*QDef1_top_sus))
			print ('The QDef_ab r_0^4 topological suceptibility is: %.5E' % (pow(self.r0_a,4)*QDef_ab_top_sus))
			output=self.EnsembleName+','+str(pow(self.r0_a,4)*spec_sum)+','+str(pow(self.r0_a,4)*spec_sum_error)+','+str(pow(self.r0_a,4)*spec_proj)+','+str(pow(self.r0_a,4)*spec_proj_error)+','+str(self.r0Mu_R) +','+ self.EnsembleName[:1] +','+ str(pow(self.r0_a,4)*QDef1_top_sus)+',' + str(pow(self.r0_a,4)*QDef_ab_top_sus)+','+str(spec_sum_total_error)+','+str(spec_proj_total_error)+','+str(Renorm_M_sqr_ev)+','+str(pow(self.r0_a,-1))+','+ str(self.r0_a_error)+',' + str(pow(self.r0_a*self.M_pi,2)) +'\n'
			f.write(output)
		f.close()
	def print_results(self,M_sqr_ev_list):
		filepath="Output/Results"+self.EnsembleName+".txt"
		f = open(filepath,"w") 
		conv_const=pow(self.Zp,2)*pow(2*self.Kappa*self.a/(197.32),2)
		M_sqr_list= [x*conv_const for x in M_sqr_ev_list]

		for M_sqr in M_sqr_list:
			a=0.5
			b=0.5
			self.Set_M_sqr(M_sqr)
			print(M_sqr)
			spec_proj=self.Renorm_Top_Suscept_SpecProjector()
			#spec_proj=x.Jack_Renorm_Top_Suscept_SpecProjector(-1)
			spec_proj=self.Renorm_Top_Suscept_SpecProjector()
			spec_sum=self.Renorm_Top_Suscept_SpecSum()
			#spec_proj_error=1.0
			spec_proj_error=self.Renorm_Top_Suscept_SpecProjector_Error()
			spec_sum_error=self.Renorm_Top_Suscept_SpecSum_Error()
			QDef1_top_sus=self.Renorm_Top_Suscept_Gaussian_Def1()
			QDef_ab_top_sus=self.Renorm_Top_Suscept_Gaussian_Def_ab(a,b)
			spec_sum_quad=pow( spec_sum_error/spec_sum,2)
			spec_proj_quad=pow( spec_proj_error/spec_proj,2)
			renorm_quad=pow(2*self.Zs_Zp_error/pow(self.Zs_Zp_ratio,2),2)
			#print(renorm_quad)
			r0_a_quad=pow(4*self.r0_a_error/pow(self.r0_a,4),2)
			#print(r0_a_quad)
			spec_sum_total_error=pow(self.r0_a,4)*spec_sum*pow(spec_sum_quad+renorm_quad+r0_a_quad,0.5)
			spec_proj_total_error=pow(self.r0_a,4)*spec_proj*pow(spec_proj_quad+renorm_quad+r0_a_quad,0.5)
			#### CAREFUL WITH THIS
			
			self.Save_Full_Eigenvalue_Collection("Output/Eigenvalue_Dist"+self.EnsembleName+".txt")
			Renorm_M_sqr_ev=M_sqr/conv_const

			print(self.EnsembleName)
			print ("Spectral Values")
			print ('The average mode number is %.5e' % (sum(self.ModeList)/len(self.ModeList)))
			print ('The spectral projector a^4 \chi is: %.5e' % (spec_proj))
			print ('The spectral projector r0^4 \chi is: %.5e' % (pow(self.r0_a,4)*spec_proj))
			print ('The spectral sum a^4 \chi is: %.5e' % (spec_sum))
			print ('The spectral sum r0^4 \chi is: %.5e' % (pow(self.r0_a,4)*spec_sum))
			print ("These are the Jackknife std. deviations")
			print ('The spectral projector r_0^4 \chi error is: %.5E' % (pow(self.r0_a,4)*(spec_proj_error)))
			print ('The spectral projector a^4 \chi error is: %.5E' % ((spec_proj_error)))
			print ('The spectral sum r_0^4 \chi error is: %.5E' %  (pow(self.r0_a,4)*spec_sum_error))
			print ('The total spectral projector r_0^4 \chi error is: %.5E' % (spec_proj_total_error))
			print ('The total spectral sum r_0^4 \chi error is: %.5E' %  (spec_sum_total_error))
			print ('The QDef1 r_0^4 topological suceptibility is: %.5E' % (pow(self.r0_a,4)*QDef1_top_sus))
			print ('The QDef_ab r_0^4 topological suceptibility is: %.5E' % (pow(self.r0_a,4)*QDef_ab_top_sus))
			output=self.EnsembleName+','+str(pow(self.r0_a,4)*spec_sum)+','+str(pow(self.r0_a,4)*spec_sum_error)+','+str(pow(self.r0_a,4)*spec_proj)+','+str(pow(self.r0_a,4)*spec_proj_error)+','+str(self.r0Mu_R) +','+ self.EnsembleName[:1] +','+ str(pow(self.r0_a,4)*QDef1_top_sus)+',' + str(pow(self.r0_a,4)*QDef_ab_top_sus)+','+str(spec_sum_total_error)+','+str(spec_proj_total_error)+','+str(Renorm_M_sqr_ev)+','+str(pow(self.r0_a,-1))+','+ str(self.r0_a_error)+',' + str(pow(self.r0_a*self.M_pi,2)) +'\n'
			f.write(output)
		f.close()
