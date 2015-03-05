#−*−coding: utf−8−*− 
"""
Created on Fri Mar 21 21:14:21 2014 SPCATPLOT for SPCAT by PICKETT
This will create Lorenztian lineshapes and simulate the output from a .cat file
@author : Zachary Glassman

INSTRUCTIONS
Take a .cat file and place it in same directory as SPCATPLOT and then
run the program. You will be prompted for the name of the file , the range (can say
all), and the gamma parameter.
The lines will be extrapolated to a gaussian, lorentzian or (maybe voigt) profile with
the following function call make_plot(filename ,window min,window max,gamma,profile type ,step)
filename: the name of the file as a string ex. ’mypredicts.cat’
window min,window max, floats setting boundaries that will be profiled. Note
that the entire spectra will still have a stick plot.
Set to a reasonable value for speed purposes gamma: parameter specificying full width at half maximum
sigma : only used for Voigt profile , controls the amount of each type convolved . 
Note: necessary for all calls , however , does not do anything for
Lorenztian or Gaussian profile profile ￼ type : Gauss , Lorentz , Voigt
NOTE: must be typed exactly as above.
step : step to be taken . Note this is heavily depednant on type of dataset .
in general , sufficiently smaller than resolving width or gamma.

V1: added Gaussian and Lorentzian Profiles
V1.1: added numpy routines for speed improvement. 
V1.2: added Voigt profile

future possibilities
− Add noise
− Multiple batch style .CAT file input −
− GUI
"""
#################################################### 
#Main Code− Don’t touch this 
################################################# 
#import proper packages
import numpy as np
import matplotlib.pyplot as plt 
from scipy import stats

from scipy import special
#fmt : This function prduces the proper number format for the plots 
#inputs x, a number
#output: a string with 5 decimal places

def fmt(x):
	return "%10.5f" %x

#get_data : gets data from .cat file
#input: file: a string with the name of the .cat file 
#output , a string of input data

def get_data(last_test): 
	data = []
	f = open(last_test, 'r')
	for line in f:
		data.append(line) 
	f.close()
	return(data)

#parse_level : parse an individual line into its consitutent parts 
	#inputs: a string of one line of the .cat file
	#outputs: [transition frequency,intensity, quantum numbers] as a list

def parse_level(line):
	freq = float(line[:13]) 
	intens = float(line[21:29]) 
	qn = line[55:].rstrip() 
	return([freq, intens ,qn])


#parse_data : parse all data and return a list of parsed data 
	#inputs: full data set from get data
	#out: list of parsed data
def parse_data(data_in): 
	parsed = []
	for i in data_in:
		parsed.append(parse_level(i)) 
	return(parsed)
	

#Lorentz : the lorenztan function
	#inputs : center , line center ; x, distance from line center;
		#gamma, profile parameter
	#outputs: values of Lorentzian at that point

def Lorentz(x,intens ,center ,gamma,sigma):
	L = np.pi * gamma * intens * stats.cauchy.pdf(x,center,gamma)
	return(L)
	

#Gauss : the gaussian function
	#inputs : center , line center ; x, array of ;
		#gamma, profile parameter
	#outputs: values of Gaussian at that point
	
def Gauss(x, intens,center,gamma,sigma):
	G = intens * np.sqrt(np.pi * 2) * gamma * stats.norm.pdf(x,center ,gamma)
	return(G)
	

#Voigt profile : convolutions of Lorenztian and Gausssian 
	#computed as V(x,sigma ,gamma)= Re(w(z))/(sigma * 2 pi) where
		#z = (x+ i * gamma)/(sigma sqrt(2)) 
	
def Voigt(x,intens ,center, gamma, sigma):
	tmp = np.real(1/special.wofz(np.zeros((len(x))) + 1j* np.sqrt(np.log(2.0))* sigma))
	V = np.real (tmp * intens * special.wofz(2 * np.sqrt(np.log (2.0) )* (x-center)/gamma + 1j * np.sqrt(np.log(2.0))*sigma))
	return(V)
	

#line_profile: add profile profile to specific line for some range 
	#inputs : line , line to be profiled , line range , range of profile ,
		#gamma, profile parameter ; profile type , type of profile 
	#return: list of values

def line_profile(line ,line_range ,gamma,sigma,profile_type): 
	freq = line [0]
	intens = pow(10,line[1])
	y = profile_type(line_range , intens , freq ,gamma, sigma) 
	return(y)
	

#add_lines : sum profiles together
	#inputs : profile type , y: a matrix of individual line profiles

def add_lines(profile_type ,y): 
	if profile_type == 'Gauss ':
		b = np.sum(y * y,axis = 0) 
		total = np.sqrt(b)
		return(total)
	else :
		total = np.sum(y,axis = 0) 
		return(total)


#total_profile:find total profile
	#inputs : parsed data , plot min , plot max ,gamma, profile type , step : step size 
	#outputs , vector to be plotted

def total_profile(parsed_data,plot_min,plot_max,gamma,sigma,profile_type ,step): 
	x = np.arange(plot_min,plot_max,step)
	data_length = len(parsed_data)
	k=0
	y = np.zeros((data_length ,len(x))) 
	for i in parsed_data :
		y[k,:] = line_profile(i,x,gamma,sigma,profile_type) 
		print('Parsed',k + 1,'profiles out of',data_length,'total') 
		k=k+1
	total = add_lines(profile_type ,y) 
	return(x,total)
	

#plot_sticks : makes stick plot
	#input : parsed array , start
	#output: stick plot of frequences which can be offset by some value

def plot_sticks(parsed,start): 
	for i in parsed :
		intens = pow(10,i[1]) 
		plt.vlines(i[0],start ,start + intens)
		

#makes overall plot, calls all other functions

def make_plot(filename,window_min,window_max,gamma,sigma,profile_type ,step):
	data = get_data(filename)
	parsed = parse_data(data)
	x, total = total_profile(parsed,window_min,window_max,gamma,sigma,profile_type,step)
	plot_sticks (parsed ,0) 
	plt.plot (x, total)


############################################## 
######################PLOT OPTIONS######################## 
fig = plt.figure()
ax = fig.add_subplot(111)

############################################## 
#Main Routine 
######################################### 
#example calls
#make plot('oP2YbF173.cat ',18104,18105,.0007,1,Gauss,.0001) 
#make plot('oP2YbF173.cat ',18104,18105,.0007,1,Lorentz ,.0001) 
#make plot('oP2YbF173.cat’,18104,18105,.0007,.0002,Voigt,.0001)

make_plot('last_test.cat',339235,339300,.0007,2,Voigt,.0001)

plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)

ax.fmt_xdata = fmt
ax.set_ylabel('Intensity')
ax.set_xlabel('Frequency') 
plt.show()
