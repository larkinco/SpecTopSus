import top_class as tc
import os
import sys
import numpy as np
#def Data_Loader_GradFlow(root_string,TwoKappaMu,Interval,Bin_Size):
#        Eigen_list=[]
#        qbfile = open(root_string+'_topological_charge.dat',"r")
#       print(os.path.join(subdir, file))
#        for data_row in qbfile:
#            values = data_row.split()
#            evl=[pow(TwoKappaMu,2)]
#            charge=[float(values[1])]
#            Eigen_list.append(tc.EigenProd(evl,charge))
#        qbfile.close()
#        return Eigen_list

#def Data_Loader_GradFlow(root_string,TwoKappaMu,Interval,Bin_Size):
#        Eigen_list=[]
#        qbfile = open(root_string+'_topological_charge.dat',"r")
#        i=0
#        diff=Bin_Size/Interval
#        charge_sqr_bin=[]
#        for data_row in qbfile:
 #           values = data_row.split()
#            evl=[pow(TwoKappaMu,2)]
#            charge_sqr_bin.append(pow(float(values[1]),2))
#            i+=1
#            if(i%(Bin_Size/Interval)==0):
#                mean=sum(charge_sqr_bin)/len(charge_sqr_bin)
#                charge=[np.sqrt(mean)]
#                Eigen_list.append(tc.EigenProd(evl,charge))
#                charge_sqr_bin=[]
#        qbfile.close()
#        return Eigen_list
	
def Data_Loader_GradFlow(root_string,TwoKappaMu,Interval,Bin_Size,Zp_Zs):
        Eigen_list=[]
        qbfile = open(root_string+'_topological_charge.dat',"r")
        i=0
#        diff=Bin_Size/Interval
        charge_sqr_bin=[]
        for data_row in qbfile:
            values = data_row.split()
            evl=[pow(TwoKappaMu,2)]
            charge_sqr_bin.append(pow(float(values[1])*Zp_Zs,2))
            i+=1
            if(i%(Bin_Size/Interval)==0):
                mean=sum(charge_sqr_bin)/len(charge_sqr_bin)
                charge=[np.sqrt(mean)]
#                charge=[mean]
#                charge=[float(values[1])]
                Eigen_list.append(tc.EigenProd(evl,charge))
                charge_sqr_bin=[]
        qbfile.close()
        return Eigen_list
	

def Data_Loader(root_string):
        Eigen_list=[]
        cutoff_list=[]
        evls=[]
        product=[]
#        print(root_string)
        for subdir, dirs, files in os.walk(root_string):
 #               print(subdir)
                for file in files:
  #                      print(os.path.join(subdir, file))
                        if(file=='v_g5_product.txt'):
                                qbfile = open(os.path.join(subdir, file),"r")
   #                             print(os.path.join(subdir, file))
                                for data_row in qbfile:
                                        values = data_row.split()
                                        if(values[0]=='V_G5_Product:'):
                                                evls.append(float(values[5]))
                                                product.append(float(values[2]))
                                        else:
                                                evls.append(float(values[4]))
                                                product.append(float(values[1]))
                                        if evls[-1]>1.0:
                                                print("EIGENVALUE IS BIG, maybe arpackproblem, use debug to avoid crashing")
                                                if not __debug__:
                                                       exit()
                                        
                                qbfile.close()
                                if __debug__:
                                        print(file,subdir)
                                        print(len(evls))
                                if len(evls)==0:
                                        print("The v_g5_product_file is empty for")
                                        print(file,subdir)
                                Eigen_list.append(tc.EigenProd(evls,product))
                                evls=[]
                                product=[]	
        return Eigen_list

def CutOff_Compute(Eig_Number,List_EigProd):
#	M_sqr= sum(cutoff_evls)/float(len(list_EigProd))	
	cut_off_sum=0.0
	for eig_prod in List_EigProd:
		cut_off_sum+=eig_prod.evls[Eig_Number]
	return cut_off_sum/float(len(List_EigProd))

def print_results(x,r0_a,EnsembleName,r0Mu_R,Zs_Zp,Zs_Zp_error,r0_a_error):
	a=0.25
	b=0.25
	spec_proj=x.Renorm_Top_Suscept_SpecProjector()
	#spec_proj=x.Jack_Renorm_Top_Suscept_SpecProjector(-1)
	spec_proj=x.Renorm_Top_Suscept_SpecProjector()
	spec_sum=x.Renorm_Top_Suscept_SpecSum()
	#spec_proj_error=1.0
	spec_proj_error=x.Renorm_Top_Suscept_SpecProjector_Error()
	spec_sum_error=x.Renorm_Top_Suscept_SpecSum_Error()
	QDef1_top_sus=x.Renorm_Top_Suscept_Gaussian_Def1()
	QDef_ab_top_sus=x.Renorm_Top_Suscept_Gaussian_Def_ab(a,b)
	spec_sum_quad=pow( spec_sum_error/spec_sum,2)
	spec_proj_quad=pow( spec_proj_error/spec_proj,2)
	renorm_quad=pow(2*Zs_Zp_error/pow(Zs_Zp,2),2)
	#print(renorm_quad)
	r0_a_quad=pow(4*r0_a_error/pow(r0_a,4),2)
	#print(r0_a_quad)
	spec_sum_total_error=pow(r0_a,4)*spec_sum*pow(spec_sum_quad+renorm_quad+r0_a_quad,0.5)
	spec_proj_total_error=pow(r0_a,4)*spec_proj*pow(spec_proj_quad+renorm_quad+r0_a_quad,0.5)
	M_sqr=x.M_sqr
	x.Save_Full_Eigenvalue_Collection("Output/Eigenvalue_Dist"+EnsembleName+".txt")

	print ("Spectral Values")
	print ('The spectral projector a^4 \chi is: %.5e' % (spec_proj))
	print ('The spectral projector r0^4 \chi is: %.5e' % (pow(r0_a,4)*spec_proj))
	print ('The spectral sum a^4 \chi is: %.5e' % (spec_sum))
	print ('The spectral sum r0^4 \chi is: %.5e' % (pow(r0_a,4)*spec_sum))
	print ("These are the Jackknife std. deviations")
	print ('The spectral projector r_0^4 \chi error is: %.5E' % (pow(r0_a,4)*(spec_proj_error)))
	print ('The spectral sum r_0^4 \chi error is: %.5E' %  (pow(r0_a,4)*spec_sum_error))
	print ('The total spectral projector r_0^4 \chi error is: %.5E' % (spec_proj_total_error))
	print ('The total spectral sum r_0^4 \chi error is: %.5E' %  (spec_sum_total_error))	
	print ('The QDef1 r_0^4 topological suceptibility is: %.5E' % (pow(r0_a,4)*QDef1_top_sus))
	print ('The QDef_ab r_0^4 topological suceptibility is: %.5E' % (pow(r0_a,4)*QDef_ab_top_sus))
	output=EnsembleName+','+str(pow(r0_a,4)*spec_sum)+','+str(pow(r0_a,4)*spec_sum_error)+','+str(pow(r0_a,4)*spec_proj)+','+str(pow(r0_a,4)*spec_proj_error)+','+str(r0Mu_R) +','+ EnsembleName[:1] +','+ str(pow(r0_a,4)*QDef1_top_sus)+',' + str(pow(r0_a,4)*QDef_ab_top_sus)+','+str(spec_sum_total_error)+','+str(spec_proj_total_error)+','+str(M_sqr)+'\n'
	filepath="Output/Results"+EnsembleName+".txt"
	f = open(filepath,"w")
	f.write(output)
	f.close()
	
