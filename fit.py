import numpy as np

class Fit:
	def __init__(self,Name):
		self.EnsembleJackList=[]
		self.EnsembleJAvgList=[]
		self.EnsembleNameList=[]
		self.Indep_Var_List=[]
		self.PointIdList=[]
		self.PointValueList=[]
		self.Name=Name
		self.EnsembleError=[]
		self.EnsembleWeight=[]
		#self.Ensemble_a_r0=[]

	def AddPoint(self,JackList,PointName,Indep_Var,Avg,point_id):
		self.EnsembleJackList.append(JackList)
		avg=sum(JackList)/len(JackList)
		#self.Ensemble_a_r0.append(a_r0)
		self.EnsembleJAvgList.append(avg)
		self.EnsembleNameList.append(PointName)
		self.Indep_Var_List.append(Indep_Var)
		self.PointIdList.append(point_id)
		self.PointValueList.append(Avg)
		self.CoeffList=[]
		self.ComputeEnsembleErrors()
		self.InterceptError=None
		self.SlopeError=None
		self.Intercept=None
		self.Slope=None
#		self.SystematicError=0

	def AddSystematicError(self,r0_a,quad):
		x=self.EnsembleError[-1]
		y=self.PointValueList[-1]
#		print("INside systematic error")
#		if(value!=y):
#			exit()
#			print("THIS IS BAD")
		error=y*pow(pow(quad,2)+pow(x,2)*pow(y,-2),0.5)
		self.EnsembleError[-1]=error
		self.EnsembleWeight[-1]=pow(error,-1)

	def MassPointJackList(self,MassPoint):
		MassPoint_Eval_Jack=[]
		for coeff in self.CoeffList:
			MassPoint_Eval_Jack.append(np.polyval(coeff,MassPoint))
		return MassPoint_Eval_Jack

	def MassPointPredwError(self,MassPoint):
		Pred=np.polyval(self.FitAvg,MassPoint)###### 
		MassPointStdErr=self.ComputeStdError(self.MassPointJackList(MassPoint))
		PredWError=(Pred,MassPointStdErr)
		return PredWError

	def ComputeEnsembleErrors(self):
		#self.EnsembleError=[]
		#self.EnsembleWeight=[]
		x=self.EnsembleJackList[-1]
	#	for x in self.EnsembleJackList:
		#	n=len(x)
		#	avg=sum(x)/len(x)
		#	error=np.sqrt(((n-1.0)/n)*sum([pow(y - avg,2) for y in x]))
		error=self.ComputeStdError(x)
		self.EnsembleError.append(error)
		self.EnsembleWeight.append(pow(error,-1))
		return self.EnsembleWeight

	def fit_data(self):
		Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.CoeffList=[]
		for j in range(0,len(self.EnsembleNameList),1):
			for x in (self.EnsembleJackList[j]):
				Temp_Jack_Mean_List[j]=x
				Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1,w=self.EnsembleWeight)
				#Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1)
				self.CoeffList.append(Coeff)
			Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.StdErr=self.ComputeStdError(self.CoeffList)
		JackAvg=sum(self.CoeffList)/len(self.CoeffList)
		self.FitAvg=np.polyfit(self.Indep_Var_List,self.PointValueList,1,w=self.EnsembleWeight)
		#print("The following are the fit average and jack average")
		n=len(self.CoeffList)
	#	print(n*self.FitAvg-(n-1)*JackAvg)
		self.InterceptError=self.StdErr[1]
		self.SlopeError=self.StdErr[0]
		self.FitCoeff=self.FitAvg# - (len(self.CoeffList)-1)*(JackAvg-self.FitAvg)	##Fix this to be unbiased	
		self.Intercept=self.FitCoeff[1]
		self.Slope=self.FitCoeff[0] 
		self.PrintResults()

	def ComputeStdError(self,Coeff):
		n=len(Coeff)
		avg= sum(Coeff)/n
		error=np.sqrt(((n-1.0)/n)*sum([pow(y-avg,2) for y in Coeff]))
		return error

	def FitChiralCont_1Slope(self,a_r0_sqr,Categorical_Var):
		import rpy2.robjects as robjects
		from rpy2.robjects import FloatVector
		from rpy2.robjects import FactorVector
		from rpy2.robjects.packages import importr
		
		self.NumSlope=1

		stats = importr('stats')
		base = importr('base')

		Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.CoeffList=[]
		robjects.globalenv["a_r0_sqr"] = FloatVector(a_r0_sqr)
		robjects.globalenv["Spacing"] = FactorVector(Categorical_Var)
		robjects.globalenv["M_pi_sqr"] = FloatVector(self.Indep_Var_List)
		r_weight=[pow(x,2) for x in self.EnsembleWeight]
		robjects.globalenv["PointValueList"] = FloatVector(self.PointValueList)
		for j in range(0,len(self.EnsembleNameList),1):
			for x in (self.EnsembleJackList[j]):
				Temp_Jack_Mean_List[j]=x
				robjects.globalenv["TempJackMean"] = FloatVector(Temp_Jack_Mean_List)
				#Coeff=(stats.lm("TempJackMean ~ Spacing:M_pi_sqr+a_r0_sqr","weights=EnsembleWeight"))[0]
				Coeff=np.asarray((stats.lm("TempJackMean ~ M_pi_sqr+a_r0_sqr",weights=FloatVector(r_weight)))[0])
				#Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1)
#				vector=numpy.asarray(vector_R)
				self.CoeffList.append(Coeff)
			Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.StdErr=self.ComputeStdError(self.CoeffList)
		JackAvg=sum(self.CoeffList)/len(self.CoeffList)
		self.FitAvg=np.asarray((stats.lm("PointValueList ~ M_pi_sqr+a_r0_sqr",weights=FloatVector(r_weight)))[0])
                #print(self.StdErr)
		#CoeffList
		self.InterceptError=self.StdErr[1]
		self.SlopeError=self.StdErr[0]
		self.FitCoeff=self.FitAvg	####Fix	
		#self.FitCoeff=self.FitAvg - (len(self.CoeffList)-1)*(JackAvg-self.FitAvg)	##Fix this to be unbiased	
		self.Intercept=self.FitCoeff[1]
		self.Slope=self.FitCoeff[0]
		self.PrintResults()
	

	def FitChiralCont(self,a_r0_sqr,Categorical_Var):
		import rpy2.robjects as robjects
		from rpy2.robjects import FloatVector
		from rpy2.robjects import FactorVector
		from rpy2.robjects.packages import importr

		stats = importr('stats')
		base = importr('base')
		
		self.NumSlope=3
		
		Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.CoeffList=[]
		robjects.globalenv["a_r0_sqr"] = FloatVector(a_r0_sqr)
		robjects.globalenv["Spacing"] = FactorVector(Categorical_Var)
		robjects.globalenv["M_pi_sqr"] = FloatVector(self.Indep_Var_List)
		r_weight=[pow(x,2) for x in self.EnsembleWeight]
		robjects.globalenv["PointValueList"] = FloatVector(self.PointValueList)
		for j in range(0,len(self.EnsembleNameList),1):
			for x in (self.EnsembleJackList[j]):
				Temp_Jack_Mean_List[j]=x
				robjects.globalenv["TempJackMean"] = FloatVector(Temp_Jack_Mean_List)
				#Coeff=(stats.lm("TempJackMean ~ Spacing:M_pi_sqr+a_r0_sqr","weights=EnsembleWeight"))[0]
				Coeff=np.asarray((stats.lm("TempJackMean ~ factor(Spacing):M_pi_sqr+a_r0_sqr",weights=FloatVector(r_weight)))[0])
				#Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1)
#				vector=numpy.asarray(vector_R)
				self.CoeffList.append(Coeff)
			Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.StdErr=self.ComputeStdError(self.CoeffList)
		JackAvg=sum(self.CoeffList)/len(self.CoeffList)
		self.FitAvg=np.asarray((stats.lm("PointValueList ~ Spacing:M_pi_sqr+a_r0_sqr",weights=FloatVector(r_weight)))[0])
                #print(self.StdErr)
		#CoeffList
		self.InterceptError=self.StdErr[1]
		self.SlopeError=self.StdErr[0]
		self.FitCoeff=self.FitAvg	####Fix	
		#self.FitCoeff=self.FitAvg - (len(self.CoeffList)-1)*(JackAvg-self.FitAvg)	##Fix this to be unbiased	
		self.Intercept=self.FitCoeff[1]
		self.Slope=self.FitCoeff[0]
		self.PrintResults()	

	def FitChiralCont_1Slope_Cross(self,a_r0_sqr,Categorical_Var):
		import rpy2.robjects as robjects
		from rpy2.robjects import FloatVector
		from rpy2.robjects import FactorVector
		from rpy2.robjects.packages import importr
		
		self.NumSlope=1

		stats = importr('stats')
		base = importr('base')

		Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.CoeffList=[]
		robjects.globalenv["a_r0_sqr"] = FloatVector(a_r0_sqr)
		robjects.globalenv["Spacing"] = FactorVector(Categorical_Var)
		robjects.globalenv["M_pi_sqr"] = FloatVector(self.Indep_Var_List)
		r_weight=[pow(x,2) for x in self.EnsembleWeight]
		robjects.globalenv["PointValueList"] = FloatVector(self.PointValueList)
		for j in range(0,len(self.EnsembleNameList),1):
			for x in (self.EnsembleJackList[j]):
				Temp_Jack_Mean_List[j]=x
				robjects.globalenv["TempJackMean"] = FloatVector(Temp_Jack_Mean_List)
				#Coeff=(stats.lm("TempJackMean ~ Spacing:M_pi_sqr+a_r0_sqr","weights=EnsembleWeight"))[0]
				Coeff=np.asarray((stats.lm("TempJackMean ~ M_pi_sqr*a_r0_sqr",weights=FloatVector(r_weight)))[0])
				#Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1)
#				vector=numpy.asarray(vector_R)
				self.CoeffList.append(Coeff)
			Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.StdErr=self.ComputeStdError(self.CoeffList)
		JackAvg=sum(self.CoeffList)/len(self.CoeffList)
		self.FitAvg=np.asarray((stats.lm("PointValueList ~ M_pi_sqr*a_r0_sqr",weights=FloatVector(r_weight)))[0])
                #print(self.StdErr)
		#CoeffList
		self.InterceptError=self.StdErr[1]
		self.SlopeError=self.StdErr[0]
		self.FitCoeff=self.FitAvg	####Fix	
		#self.FitCoeff=self.FitAvg - (len(self.CoeffList)-1)*(JackAvg-self.FitAvg)	##Fix this to be unbiased	
		self.Intercept=self.FitCoeff[1]
		self.Slope=self.FitCoeff[0]
		self.PrintResults()


	def FitChiralCont_1Slope_Cross_wSubtraction(self,a_r0_sqr,Categorical_Var,filepath):
		import rpy2.robjects as robjects
		from rpy2.robjects import FloatVector
		from rpy2.robjects import FactorVector
		from rpy2.robjects.packages import importr
		
		self.NumSlope=1

		stats = importr('stats')
		base = importr('base')

		Temp_Jack_Mean_List=self.EnsembleJAvgList
		SubtractionJackList=[]
		self.CoeffList=[]
		robjects.globalenv["a_r0_sqr"] = FloatVector(a_r0_sqr)
		robjects.globalenv["Spacing"] = FactorVector(Categorical_Var)
		robjects.globalenv["M_pi_sqr"] = FloatVector(self.Indep_Var_List)
		r_weight=[pow(x,2) for x in self.EnsembleWeight]
		robjects.globalenv["PointValueList"] = FloatVector(self.PointValueList)
		for j in range(0,len(self.EnsembleNameList),1):
			for x in (self.EnsembleJackList[j]):
				Temp_Jack_Mean_List[j]=x
				robjects.globalenv["TempJackMean"] = FloatVector(Temp_Jack_Mean_List)
				#Coeff=(stats.lm("TempJackMean ~ Spacing:M_pi_sqr+a_r0_sqr","weights=EnsembleWeight"))[0]
				Coeff=np.asarray((stats.lm("TempJackMean ~ M_pi_sqr*a_r0_sqr",weights=FloatVector(r_weight)))[0])
				#Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1)
#				vector=numpy.asarray(vector_R)
				self.CoeffList.append(Coeff)
				SubtractionList=[]
				for k in range(0,len(a_r0_sqr),1):
					SubtractionList.append(Temp_Jack_Mean_List[k]-a_r0_sqr[k]*Coeff[2] - a_r0_sqr[k]*self.Indep_Var_List[k]*Coeff[3])
		#			SubtractionList.append(Temp_Jack_Mean_List[k])
				SubtractionJackList.append(SubtractionList)
			Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.StdErr=self.ComputeStdError(self.CoeffList)
		JackAvg=sum(self.CoeffList)/len(self.CoeffList)
		self.FitAvg=np.asarray((stats.lm("PointValueList ~ M_pi_sqr*a_r0_sqr",weights=FloatVector(r_weight)))[0])

		Rev_SubtractionJackList=[]
#		print(len(SubtractionJackList[]))
		SubtractJackAvg=[]
		SubtractJackErr=[]
		for i in range(0,len(a_r0_sqr),1):
			Temp=[]
			for j in range(0,len(SubtractionJackList),1):
				Temp.append(SubtractionJackList[j][i])
			Rev_SubtractionJackList.append(Temp)
			#SubtractJackAvg.append(sum(Rev_SubtractionJackList[i])/len(Rev_SubtractionJackList))
			#SubtractJackErr.append(self.ComputeStdError(Rev_SubtractionJackList[i]))
		print(len(Rev_SubtractionJackList[1]))
		#print(self.StdErr)
		#CoeffList
		
		for x in Rev_SubtractionJackList:
			#SubtractJackAvg.append(sum(x)/len(x))
#			print(len(x))
#                        SubtractJackAvg.append()
			SubtractJackErr.append(self.ComputeStdError(x))
		#SubtractJackAvg=self.PointValueList
		for k in range(0,len(a_r0_sqr),1):
			SubtractJackAvg.append(self.PointValueList[k]-a_r0_sqr[k]*self.FitAvg[2] - a_r0_sqr[k]*self.Indep_Var_List[k]*self.FitAvg[3])
			
#		print(self.FitAvg)	

		self.InterceptError=self.StdErr[1]
		self.SlopeError=self.StdErr[0]
		self.FitCoeff=self.FitAvg	####Fix	
		#self.FitCoeff=self.FitAvg - (len(self.CoeffList)-1)*(JackAvg-self.FitAvg)	##Fix this to be unbiased	
		self.Intercept=self.FitCoeff[1]
		self.Slope=self.FitCoeff[0]
		self.PrintResults()
	
		f=open(filepath,"w")
		output='PointName,TopSus,TopSus_Error,r0_Mpi_sqr,PlotName,PointId\n'
		f.write(output)
		for i in range(0,len(self.EnsembleNameList),1):
			output=self.EnsembleNameList[i]+','+str(SubtractJackAvg[i])+','+str(SubtractJackErr[i])+','+str(self.Indep_Var_List[i])+','+str(self.Name)+','+str(self.PointIdList[i])+'\n'
			f.write(output)
		f.close()

	def FitChiralCont_1Slope_Cross_wSubtraction_NoIntercept(self,a_r0_sqr,Categorical_Var,filepath):
		import rpy2.robjects as robjects
		from rpy2.robjects import FloatVector
		from rpy2.robjects import FactorVector
		from rpy2.robjects.packages import importr
		
		self.NumSlope=1

		stats = importr('stats')
		base = importr('base')

		Temp_Jack_Mean_List=self.EnsembleJAvgList
		SubtractionJackList=[]
		self.CoeffList=[]
		robjects.globalenv["a_r0_sqr"] = FloatVector(a_r0_sqr)
		robjects.globalenv["Spacing"] = FactorVector(Categorical_Var)
		robjects.globalenv["M_pi_sqr"] = FloatVector(self.Indep_Var_List)
		r_weight=[pow(x,2) for x in self.EnsembleWeight]
		robjects.globalenv["PointValueList"] = FloatVector(self.PointValueList)
		for j in range(0,len(self.EnsembleNameList),1):
			for x in (self.EnsembleJackList[j]):
				Temp_Jack_Mean_List[j]=x
				robjects.globalenv["TempJackMean"] = FloatVector(Temp_Jack_Mean_List)
				#Coeff=(stats.lm("TempJackMean ~ Spacing:M_pi_sqr+a_r0_sqr","weights=EnsembleWeight"))[0]
				Coeff=np.asarray((stats.lm("TempJackMean ~0.0+ M_pi_sqr*a_r0_sqr",weights=FloatVector(r_weight)))[0])
				#Coeff=np.polyfit(self.Indep_Var_List,Temp_Jack_Mean_List,1)
#				vector=numpy.asarray(vector_R)
				self.CoeffList.append(Coeff)
				SubtractionList=[]
				for k in range(0,len(a_r0_sqr),1):###### FIX
					SubtractionList.append(Temp_Jack_Mean_List[k]-a_r0_sqr[k]*Coeff[1] - a_r0_sqr[k]*self.Indep_Var_List[k]*Coeff[2])
		#			SubtractionList.append(Temp_Jack_Mean_List[k])
				SubtractionJackList.append(SubtractionList)
			Temp_Jack_Mean_List=self.EnsembleJAvgList
		self.StdErr=self.ComputeStdError(self.CoeffList)
		JackAvg=sum(self.CoeffList)/len(self.CoeffList)
		self.FitAvg=np.asarray((stats.lm("PointValueList ~ 0+ M_pi_sqr*a_r0_sqr",weights=FloatVector(r_weight)))[0])

		Rev_SubtractionJackList=[]
#		print(len(SubtractionJackList[]))
		SubtractJackAvg=[]
		SubtractJackErr=[]
		for i in range(0,len(a_r0_sqr),1):
			Temp=[]
			for j in range(0,len(SubtractionJackList),1):
				Temp.append(SubtractionJackList[j][i])
			Rev_SubtractionJackList.append(Temp)
			#SubtractJackAvg.append(sum(Rev_SubtractionJackList[i])/len(Rev_SubtractionJackList))
			#SubtractJackErr.append(self.ComputeStdError(Rev_SubtractionJackList[i]))
		print(len(Rev_SubtractionJackList[1]))
		#print(self.StdErr)
		#CoeffList
		
		for x in Rev_SubtractionJackList:
			#SubtractJackAvg.append(sum(x)/len(x))
#			print(len(x))
#                        SubtractJackAvg.append()
			SubtractJackErr.append(self.ComputeStdError(x))
		#SubtractJackAvg=self.PointValueList
		for k in range(0,len(a_r0_sqr),1):###FIX
			SubtractJackAvg.append(self.PointValueList[k]-a_r0_sqr[k]*self.FitAvg[1] - a_r0_sqr[k]*self.Indep_Var_List[k]*self.FitAvg[2])
			
#		print(self.FitAvg)	

		self.InterceptError=self.StdErr[1]
		self.SlopeError=self.StdErr[0]
		self.FitCoeff=self.FitAvg	####Fix	
		#self.FitCoeff=self.FitAvg - (len(self.CoeffList)-1)*(JackAvg-self.FitAvg)	##Fix this to be unbiased	
		self.Intercept=self.FitCoeff[1]
		self.Slope=self.FitCoeff[0]
		self.PrintResults()
	
		f=open(filepath,"w")
		output='PointName,TopSus,TopSus_Error,r0_Mpi_sqr,PlotName,PointId\n'
		f.write(output)
		for i in range(0,len(self.EnsembleNameList),1):
			output=self.EnsembleNameList[i]+','+str(SubtractJackAvg[i])+','+str(SubtractJackErr[i])+','+str(self.Indep_Var_List[i])+','+str(self.Name)+','+str(self.PointIdList[i])+'\n'
			f.write(output)
		f.close()


	def JackknifedSubtraction(self,Data_a_r0_sqr,filepath):
		DataCoeff=[]
		DataStdDev=[]
		SubtractedData=[]
		Data=self.PointValueList
		a_index=0
		if(self.NumSlope==1):
			a_index=1
		if(self.NumSlope==3):
			a_index=2
		if(a_index==0):
			exit()

		for i in range(0,len(Data),1):
			for Coeff in self.CoeffList:
				DataCoeff.append(Data[i]-Coeff[a_index]*Data_a_r0_sqr[i])
			DataStdDev.append(self.ComputeStdError(DataCoeff))
			SubtractedData.append(Data[i]-self.FitCoeff[a_index]*Data_a_r0_sqr[i])
			DataCoeff=[]

		f=open(filepath,"w")
		output='PointName,TopSus,TopSus_Error,r0_Mpi_sqr,PlotName,PointId\n'
		f.write(output)
		for i in range(0,len(self.EnsembleNameList),1):
			output=self.EnsembleNameList[i]+','+str(SubtractedData[i])+','+str(DataStdDev[i])+','+str(self.Indep_Var_List[i])+','+str(self.Name)+','+str(self.PointIdList[i])+'\n'
			f.write(output)
		f.close()
			
#		DataTup = (SubtractedData,SubtractedStdDev);
#		return DataTup
	    

	def plot():
	    temp=0

	def PrintResults(self):
		print("The coefficients of the "+self.Name+  " fit are, slope and intercept respectively,")
		print(self.FitCoeff)
		print("The standard error of the "+self.Name+  " fit are, slope and intercept respectively,")
		print(self.StdErr)

	def SaveResults(self,filepath):
		f=open(filepath,"w")
		#print("Ensemble Error")
		#print(self.EnsembleError)
		for i in range(0,len(self.EnsembleNameList),1):
			output=self.EnsembleNameList[i]+','+str(self.PointValueList[i])+','+str(self.EnsembleError[i])+','+str(self.Indep_Var_List[i])+','+str(self.Intercept) +','+str(self.InterceptError)+','+str(self.Name)+','+str(self.PointIdList[i])+'\n'
			f.write(output)
		f.close()
        
