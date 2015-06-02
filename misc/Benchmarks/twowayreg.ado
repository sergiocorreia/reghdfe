capture program drop twowayset
capture mata mata drop sparse()
capture mata mata drop proddiag()
capture mata mata drop diagprod()
capture mata mata drop diagminus()
capture mata mata drop projDummies()
capture mata mata drop saveMat()
capture mata mata drop readMat()
//Mata programs:

mata:
real matrix sparse(real matrix x)
 {
  real matrix y
  real scalar k
 
  y = J(colmax(x[,1]),colmax(x[,2]),0)
  for (k=1; k<=rows(x); k++) {
    y[x[k,1],x[k,2]] = y[x[k,1],x[k,2]] + x[k,3]
  }
 
  return(y)
 }
 
 //sparse matrix function ends
 
 
 // multiplying a diagonal matrix represented by a vector times a matrix.
 // Diag*A multiplies each rows.
 real matrix diagprod(real colvector x, real matrix A)
 {
  real matrix y
  real scalar k
  if(rows(x)<cols(x)) x = x'
 
  y = J(rows(A),cols(A),0)
  for (k=1; k<=rows(x); k++) {
    y[k,] = A[k,] * x[k,1]
  }
 
  return(y)
 }
 
 real matrix readMat(string s,string n)
 {

fh = fopen(s+"_"+n, "r")
X = fgetmatrix(fh)
fclose(fh)
return(X)
 }
 
 void saveMat(string s,string n,real matrix X)
 {

fh = fopen(s + "_" + n, "rw")
fputmatrix(fh, X)
fclose(fh)
 }
 
  
 
  real matrix proddiag(real matrix A,real colvector x)
 {
  real matrix y
  real scalar k
  if(rows(x)<cols(x)) x = x'
 
  y = J(rows(A),cols(A),0)
  for (k=1; k<=rows(x); k++) {
    y[,k] = A[,k] * x[k,1]
  }
 
  return(y)
 }
 
   real matrix diagminus(real colvector x,real matrix A)
 {
  //real matrix y
  real scalar k
  if(rows(x)<cols(x)) x = x'
 
  //y = -A
  for (k=1; k<=rows(x); k++) {
    A[k,k] = A[k,k] - x[k,1]
  }
 
  return(-A)
 }

void projDummies()
{
real matrix D, DH1, DH, CinvHHDH, AinvDDDH, A, B, C
real colvector DD, HH, invDD, invHH
real scalar N, T
string scalar id, t, w,sampleVarName
D=.
//printf("Hola Paulo, todo functiona hasta aqui.")

id = st_local("twoway_id")
t = st_local("twoway_t")
w = st_local("twoway_w")
root =st_local("root")
sampleVarName = st_local("twoway_sample")
if (w==""){
D = st_data(.,(id,t),sampleVarName)
D = (D,J(rows(D),1,1))
}
else {
D = st_data(.,(id,t,w),sampleVarName)
}
//printf(sampleVarName)
//printf("Incluso aca\n")
//D[1..10,]
//printf("y aca")

DH1=sparse(D)
//printf("Wohoo")
DD=quadrowsum(DH1)
HH=quadcolsum(DH1)'
HH=HH[1..cols(DH1)-1]



DH=DH1[.,1..cols(DH1)-1]

 
invDD=DD:^-1 
invHH=HH:^-1

N=colmax(D)[.,1]
T=colmax(D)[.,2]
saveMat(root,"twoWayN1", N)
saveMat(root,"twoWayN2", T)
saveMat(root,"twoWayinvDD", invDD)
saveMat(root,"twoWayinvHH", invHH)
//st_matrix("twoWayD", D...)
 if (N<T)
		{
        
        CinvHHDH=diagprod(invHH,DH')
		A=qrinv(diagminus(DD,CinvHHDH'*DH'))
		//st_matrix("CinvHHDH",CinvHHDH)
        B=-A*CinvHHDH'
		saveMat(root,"twoWayCinvHHDH", CinvHHDH)
		saveMat(root,"twoWayA", A)
		saveMat(root,"twoWayB", B)
		
		
		}
    else
	{
        AinvDDDH=diagprod(invDD,DH)
		C=qrinv(diagminus(HH,AinvDDDH'*DH))
		//st_matrix("AinvDDDH",AinvDDDH)
        B=-AinvDDDH*C
		saveMat(root,"twoWayAinvDDDH", AinvDDDH)
		saveMat(root,"twoWayC", C)
		saveMat(root,"twoWayB", B)
		
    }
 }
 
 end



program define twowayset, rclass
version 11
syntax varlist(min=2 max=3) [if] [in], [Root(name)]
//summ `varlist'
// I need to make it robust to non 1,2,3... ids.
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux
if ("`root'" == "") {
	local root="last"
	}

//di in gr "`twoway_id'"
//di in gr "`twoway_t'"

tempvar twoway_sample
mark `twoway_sample' `if' `in'
markout `twoway_sample' `varlist'
mata projDummies()
//di in gr "Checkpoint 1"
//ret li
//di in gr "Checkpoint 2"
scalar twoWayid="`twoway_id'"
scalar twoWayt="`twoway_t'"
scalar twoWayw="`twoway_w'"
scalar twoWayif="`if'"
scalar twoWayin="`in'"
//return post r(B), esample(`twoway_sample') 
//obs(`nobs') dof(`dof')

end
 

capture program drop projvar
capture mata mata drop projVar()

mata
void projVar()
{
	real matrix V, varIn, D,aux,delta,tau,varOut,A,B,CinvHHDH,AinvDDDH,C
	real colvector invHH,invDD,Dy,Ty
	real scalar N,T
	string scalar id, t, currvar,newvar,sampleVarName,w
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	id=st_strscalar("twoWayid")
	root =st_local("root")
	N=readMat(root,"twoWayN1")
	T=readMat(root,"twoWayN2")
	//D=readMat(root,"twoWayD")
	w=st_strscalar("twoWayw")
	t=st_strscalar("twoWayt")
	sampleVarName = st_local("twoway_sample")
	V = st_data(.,(id,t,currvar),sampleVarName)
	varIn=V[.,3]
	
	if (w==""){
	D = st_data(.,(id,t),sampleVarName)
	D = (D,J(rows(D),1,1))
	}
	else {
	D = st_data(.,(id,t,w),sampleVarName)
	}
	
	V[.,3]=V[.,3]:*D[.,3]
	aux=sparse(V)
	//printf("3")
	Dy=rowsum(aux)
	Dy=Dy
	Ty=colsum(aux)
	Ty=Ty[1,1..cols(aux)-1]'
	B=readMat(root,"twoWayB")
	
	//rows(Ty)
    //cols(Ty)
	//rows(Dy)
	//cols(Dy)
			

	 if (N<T)
			{
			
			A=readMat(root,"twoWayA")
			invHH=readMat(root,"twoWayinvHH")
			CinvHHDH=readMat(root,"twoWayCinvHHDH")
			//printf("b")
			delta=A*Dy+B*Ty
			tau=B'*(Dy-CinvHHDH'*Ty)+(invHH:*Ty) \0
			}
		else
		{
			//printf("1")
			C=readMat(root,"twoWayC")
			invDD=readMat(root,"twoWayinvDD")
			AinvDDDH=readMat(root,"twoWayAinvDDDH")
			delta=(invDD:*Dy)+B*(Ty-AinvDDDH'*Dy)
			tau=B'*Dy+C*Ty \0 
			//printf("c")
		}

	//how to index
	//varout=(var-delta(struc.hhid)-tau(struc.tid')).*sqrt(struc.w);
	varOut=(varIn-delta[V[.,1]]-tau[V[.,2]]):*sqrt(D[.,3])
	//printf("4")
	//st_matrix("DD2",B)
	st_store(., newvar, varOut)
	//printf("5")
}
end


program define projvar, nclass
version 11
syntax varlist, [Prefix(name)] [Root(name)] [REPLACE]
tempvar twoway_sample
loc tif=twoWayif
loc tin=twoWayin
mark `twoway_sample' `tif' `tin'
markout `twoway_sample' `varlist'
//mata mata describe
//summ `varlist'
//summ `twoway_sample'
// I need to make it robust to non 1,2,3... ids.
if ("`prefix'" == "") {
	local prefix="proj_"
	}
if ("`root'" == "") {
	local root="last"
	}

foreach currvar of varlist `varlist' {
	local newvar="`prefix'`currvar'"
	if ("`replace'" != "") {
	local newvar="`currvar'"
	}
	else {
	gen `newvar'=.
	}
	//di "`currvar'"
	//di "`newvar'"
	mata projVar()
	/*
	mata
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	printf(".")
	V = st_data(.,(id,t,currvar),sampleVarName)
	varIn=V[.,3]
	V[.,3]=V[.,3]:*D[.,3]
	aux=sparse(V)
	printf(".")
	Dy=rowsum(aux)
	Ty=colsum(aux)
	Ty=Ty[1,1..cols(aux)-1]'

	 if (N<T)
			{
			delta=A*Dy+B*Ty
			tau=B'*(Dy-CinvHHDH'*Ty)+invHH*Ty \0
			}
		else
		{
			delta=(invDD:*Dy)+B*(Ty-AinvDDDH'*Dy)
			tau=B'*Dy+C*Ty \0 
			
		}

	//how to index
	//varout=(var-delta(struc.hhid)-tau(struc.tid')).*sqrt(struc.w);
	varOut=(varIn-delta[V[.,1]]-tau[V[.,2]]):*sqrt(D[.,3])
	printf(".")
	//st_matrix("DD2",B)
	st_store(., newvar, varOut)
	printf(".")
	end
	*/
}


//gettoken twoway_id aux: varlist
//gettoken twoway_t twoway_w: aux

//di in gr "`twoway_id'"
//di in gr "`twoway_t'"

//tempvar twoway_sample
//mark `twoway_sample' `if' `in'
//markout `twoway_sample' `varlist'
//mata projDummies()
//di in gr "Checkpoint 1"
//ret li
//di in gr "Checkpoint 2"
//return add
//return post r(B), esample(`twoway_sample') 
//obs(`nobs') dof(`dof')

end
 
