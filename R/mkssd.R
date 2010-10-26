#m=number of factors
#n=number of runs
#q=number of levels
#D is the supersaturated design

genvec=function(m,n,q,k)
#The function created initial k-circulant q-level SSD randomly.
{
	t=n/q	#checks condition 1
	if (n%%q !=0) stop("n should be divisible by q")
	if (m!=(n-1)*k) stop("k-circulant design not possible") #checks condition 1	
	gvec=matrix(0,1,m)
	gmat=matrix(0,n-1,k)
	for (i in 1:(q-1))
	{
		for (j in 1:t)
		{
			gmat[j+(i-1)*t,1]=i
		}
	}
	for (i in (t*(q-1)+1):(n-1))
	{
			gmat[i,1]=q
	}
	for (i in 1:k)
	{
		gmat[,i]=sample(gmat[,1],(n-1))
	}	#condition 2
	count=0
	for (j in 1:(n-1))
	{
		for (i in 1:k)
		{
			count=count+1
			gvec[,count]=gmat[j,i]
		}
	}	
	gvec
}

D_from_G=function(genvec,m,n,q,k)
{
#the function creates the design obtained by k-circulating the generator vector
	D=matrix(0,n,m)
	D[1,]=genvec
	#k-circulating the first row to next (n-2) rows
	for (rown in 1:(n-2))
	{
		for (coln in 1:m)
		{
			coln_plusk=coln+k
			if ((coln+k)>m) coln_plusk= (coln+k) %% m
			D[rown+1,coln_plusk]=D[rown,coln]
		}
	}	
	#Allocation of level q to last row of the design
	for (i in 1:m) D[n,i]=q
	D
}


ave_chisq=function(D,m,n,q)
{
	#The function calculates average chisquare of the design
	aliased=0	
	count=0
	sumchi=0	
	max_chisq=0	
	for (i in 1:(m-1))
	{
		for (j in (i+1):m)
		{
			freq_levels=matrix(-n/(q*q),q,q) #freq_levels is the matrix to contain freq of level pairs in the columns i and j
			for (i1 in 1:n)
			{
				freq_levels[D[i1,i],D[i1,j]]=freq_levels[D[i1,i],D[i1,j]]+1
			}
			for (a in  1:q)
			{
				for (b in 1:q)
				{
					freq_levels[a,b]=(freq_levels[a,b])**2
				}
			}
			temp1=sum(freq_levels)/(n/(q*q))
			if (temp1==n*(q-1)) aliased=aliased + 1	#checks the aliasness of two columns i and j and counts number of such pairs					
			sumchi=sumchi+temp1			
			if (temp1>max_chisq) max_chisq=temp1						
		}
	}	
	ave_chisq=sumchi*2/(m*(m-1))	
	output1=list(ave_chisq,max_chisq,aliased)
}

Design_swap=function(D,m,n,k,i,j)
{
#The function generates the design after swaping the column i and j of the generator vector
	t_i=D[1,i]
	t_j=D[1,j]
		
	for (row in 1:(n-1))
	{
		coli_plusk=i+(row-1)*k
		colj_plusk=j+(row-1)*k
		if ((coli_plusk)>m) coli_plusk= (coli_plusk) %% m
		if ((colj_plusk)>m) colj_plusk= (colj_plusk) %% m
		D[row,coli_plusk]=t_j
		D[row,colj_plusk]=t_i		
	}
D
}	

mkssd=function(m,n,q,k,mef)
{
#The function tries to obtain efficient multilevel supersaturated design that has efficiency more than 'mef' by interchanging factors of the generator vector. Here, first position of
#the generator vector is interchanges with most important factor, then second position with most important factor and so on.
#After full iteration, design may not have chisquare efficiency 1. In that case recall the function. mef is the minimum efficiency required, should be between 0 to 1.
	#time of execution
	stime=proc.time()
	L_chisq=(q-1)*n*((q-1)*m-n+1)/((n-1)*(m-1))
	aliased=100000
	trial=0
	while (aliased>0 && trial<=100)  #aliased=1 means two columns in the design are fully aliased
	{
		trial=trial+1
		aliased=0
		gvec=genvec(m,n,q,k)
		D=D_from_G(gvec,m,n,q,k)		
		outersuccess=1
		final_success=0 #final_success=1 implies a required design was found
		while (outersuccess!=0)  #outersuccess=1 implies that after a full round of iteration the criterion value decreased and a better design was found
		{					
			outersuccess=0
			outer_avechi=ave_chisq(D,m,n,q)[[1]]
			loop=0			
			total=k*(n-1)*(n-2)/2
			pb <- txtProgressBar(min = 0, max = total, style=3)	#progress bar variable
			i=1				
			while (i<=(m-1) && (final_success==0) )
			{
				avechi=ave_chisq(D,m,n,q)[[1]]
				if (mef*avechi<=L_chisq) final_success=1		
				D1=D
				min_avechi=avechi				
				j=i+k			
				while ((j<=m) && (final_success==0))
				{
					loop=loop+1 #loop is used for creating progress bar
					Sys.sleep(0.1)
   					setTxtProgressBar(pb, loop) #creates progress bar
					if (D[1,i]!=D[1,j])
					{
						temp_gvec=D[1,]					
						temp_gvec[i]=D[1,j]
						temp_gvec[j]=D[1,i]						
						D_temp=Design_swap(D,m,n,k,i,j)
						avechi_temp=ave_chisq(D_temp,m,n,q)[[1]]						
						if (avechi_temp<min_avechi) 
						{
							min_avechi=avechi_temp
							D1=D_temp
							success=1							
						}
						if (mef*min_avechi<=L_chisq) final_success=1 #for terminating the outer while loop
					}
					j=j+k				
				}
				D=D1
				i=i+1			
			}
			close(pb) #closes progress bar
			if ((min_avechi<outer_avechi) && (final_success==0))
			{
				outersuccess=1
				outer_avechi=min_avechi
			}			 		
		}
		aliased=ave_chisq(D,m,n,q)[[3]]
	}	
	Deff=L_chisq/min_avechi	
	max_chisq=ave_chisq(D,m,n,q)[[2]]	
	t_taken=proc.time()-stime  #t_taken is time taken to generate the design
	genv=D[1,]
	if (q==2) 
	{
		LB=n*n*(m-n+1)/((n-1)*(m-1))
		if (n*k %% 4==0) LB=n*n*(k-1)/(n*k-k-1)
		h=(k*(k-1)*(n-1)*n*n+2*n*(n-2))/(k*(n-1)*(k*n-k-1))		
		if ((n %% 4==2) && ( k %% 2== 1)) LB = 4+64*ceiling(m*(m-1)*(h-4)/64)/(m*(m-1))
		Es2=n*min_avechi
		Deff=LB/Es2
		for (i in 1:m)
		{
			if (genv[i]==1) genv[i]=-1
			if (genv[i]==2) genv[i]=1
		}
		result=list(m,n,q,k,genv,Deff,t_taken, aliased)	
	}
	else
	{
		result=list(m,n,q,k,genv,Deff,max_chisq,t_taken, aliased)
	}
	return(result)	
}