#!/home/ocean/anaconda3/bin/python3
from numpy import cos, arccos, sin, arctan, tan, pi, sqrt; from numpy import array as ary; import numpy as np; tau = 2*pi
from matplotlib import pyplot as plt
#Program for calculating the coefficient of heat transfer in the given reactors in the exercise at the bottom of slide set number 4.

class Reactor:
	def __init__(self,coolant):
		if coolant=="water" or coolant=="H2O":
			self.k = 0.5
			self.rho=704
			self.mu =8.69E-5
			self.cp =6270
			self.Pr = self.cp*self.mu/self.k
			self.P  = 12.6E-3
			self.D  = 9.5E-3
		elif coolant=="sodium" or coolant=="Na":
			self.k = 62.6
			self.rho=817.7
			self.mu =2.28E-4
			self.cp =1254
			self.Pr = self.cp*self.mu/self.k
			self.P  = 9.8E-3
			self.D  = 8.5E-3
def FindInLocals(DICT,varnames):
	for name in varnames:
		if not name in DICT:
			if not name in DICT["kwargs"]:
				print("error: can't find the variable", name)
				exit()

def getNu(method,**kwargs):
	default_method="Dittus-Boelter"

	if method in ["Kazimi","square array","triangular array"]:#For calculating array case Nu
		FindInLocals(locals(),["P","D"])
		if method == "Kazimi":
			FindInLocals(locals(),["Pe"])
			Nu = 4.0+ 0.33*(P/D)**(3.8)*(Pe/100)**(0.86)+0.16*(P/D)**5.0
		
		elif method=="square array":
			phi= 0.9217+0.1478*(P/D)-0.1130*np.exp(-7*((P/D)-1))
			if default_method=="Dittus-Boelter": FindInLocals(locals(),["Re","Pr"])
			Nu = phi*getNu(default_method,Pr=Pr, Re=Re,Pe=Pe)
		
		elif method=="triangular array":
			phi= 0.9090+0.0783*(P/D)-0.1283*np.exp(-2.4*( (P/D)-1) )
			if default_method=="Dittus-Boelter": FindInLocals(locals(),["Re","Pr"])
			Nu = phi*getNu(default_method,Pr=Pr, Re=Re,Pe=Pe)
	
	elif method in ["Lyon", "Seban", "Dittus-Boelter"]:#for calculating the cylindrical case Nu
		if method=="Lyon":
			FindInLocals(locals(),["Pe"])
			Nu = 7+0.025*Pe**0.8
		elif method=="Seban":
			FindInLocals(locals(),["Pe"])
			Nu = 5.0+0.025*Pe**0.8
		elif method=="Dittus-Boelter":
			FindInLocals(locals(),["Re","Pr"])
			Nu = 0.023*Re**0.8*Pr**0.4
	else:
		print("method invalid!")
		exit()
	return Nu

def getRe(rho,u,D,mu):#Reynolds number
	return rho*u*D/mu

def getPe(Re,Pr):#Peclet number
	return Re*Pr

def getPr(cp,mu,k):#Prandtl's number
	return cp*mu/k

def get_h(Nu,k,De):#get convective coefficient given the numerical value of Nu and the array type
	return Nu*k/De
def getDe(array_type,P,D):
	#Determining the factor to put inside the expression,
	#obtained after evaluating the equation
	#De = Pw/(Af/4)
	if array_type=="square":
		factor=4
	elif array_type=="triangular":
		factor=2*sqrt(3)
	De = D*( factor/pi*(P/D)**2 -1)
	return De
def getAttr(reactor):
	assert type(reactor)==Reactor, "expected reactor input"
	return reactor.k ,reactor.rho ,reactor.mu ,reactor.cp ,reactor.Pr ,reactor.P ,reactor.D 
if __name__=="__main__":
	PWR=Reactor("water")
	SFBR=Reactor("sodium")
	Rxr=input("Reactor type?")
	v  =float(input("flow velocity?"))
	if Rxr=="PWR":
		k,rho,mu,cp,Pr,P,D = getAttr(PWR)
		De=getDe("square",P=P,D=D)
		Re=getRe(rho=rho,u=v,D=De,mu=mu)
		Pe=getPe(Re,Pr)
		Nu=getNu("square array",P=P,D=D,Re=Re,Pr=Pr,Pe=Pe)
	elif Rxr=="SFBR":
		k,rho,mu,cp,Pr,P,D = getAttr(SFBR)
		De=getDe("triangular",P=P,D=D)
		Re=getRe(rho=rho,u=v,D=De,mu=mu)
		Pe=getPe(Re,Pr)
		Nu=getNu("Kazimi",P=P,D=D,Re=Re,Pr=Pr,Pe=Pe)
	print("De=",De)
	print("Re=",Re)
	print("Pe=",Pe)
	print("Nu=",Nu)
		#use "square array" for PWR
		#use "Kazimi" for SFBR
	h =get_h(Nu,k,De)
	print(h/1000)