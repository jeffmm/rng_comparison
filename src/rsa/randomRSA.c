#include "randomRSA.h"

//private functions
static int vrandomRSA0();
static uint64_t getMicrotime();
static int primeq(uint64_t n);
static int safeprimeq(uint64_t n);
static uint64_t nextprime(uint64_t n, long offset);
static uint64_t nextsafeprime(uint64_t n, long offset);
static uint32_t safeprime(long i);
static uint64_t powermod(uint64_t x, uint64_t e, uint64_t n);
static inline uint32_t powermod32(uint32_t x, uint32_t e, uint32_t n);;
static uint64_t gcd(unsigned long u, unsigned long v);
static uint64_t modinverse(uint64_t u, uint64_t v);

//validate the operation of the vectorized pseudorandom number generator
static int randomRSA_validate(int numvectorcalls);

//private data

//well-tested primitive roots mod q, where q=2^63-25 is largest prime less than 2^63 
//(#0 and #10 are repeated because of trusted source (P. L’Ecuyer. F.-O. Blouin, 
// and R. Couture) and need to have prime number (11)). All others are Sezgin and Sezgin
static uint32_t two63minus25primitiveroots[11]={2307085864,3157107955,3163786287,3200261722,3211103532,3338736601,3423977237,3465965455,3474009732,3512424704,2307085864};

//safe primes less than 2^32: safeprime[0] is largest, safeprime[1] is 100000th largest etc.
static uint32_t safeprimes[31]={4294967087, 4222725707, 4150883063, 4079008547, 4006857719, 3935091059, 3863505827, 3792054359, 3720601559, 3649365899, 3578380259, 3507339959, 3436690847, 3366364187, 3295757987, 3225327167, 3155044799, 3085122659, 3015377447, 2945577263, 2876115563, 2806385207, 2736992087, 2667700019, 2598759287, 2530175939, 2461716743, 2393244167, 2324812307, 2256486719, 2188727483};

//time stamp used in initializing generator
static uint64_t timemicrosecondssince1970;

//prime modulus for skip generator, q=2^63-25, and default values of oter parameters
static uint64_t q=9223372036854775783UL;
static uint32_t a=2307085864;   //restricted primitive root mod , a<sqrt(q)
static uint32_t q1=2915295535;   //=q/a;  used in calculation of (a s) mod q
static uint32_t q2=669447238;   //=q%a;  used in calculation of (a s) mod q

static uint32_t p1=4294967087; //larger of two primes selected from [2^31.5 .. 2^32]
static uint32_t p2=2147483783; //smaller of two primes in [2^31 .. 2^31.5] near q/p1
static uint64_t n=9223372167851250121UL; //RSA composite modulus n=p1*p2 close to q
static uint32_t e=DEFAULTEXPONENT,exponentsaved=DEFAULTEXPONENT;  //RSA exponent coprime to (p1-1)*(p2-1)
static uint32_t p2inversemodp1=932518950;   //=p2^-1 mod p1 used in Garner's formula in Chinese Remainder Theorem calculations

static uint64_t s[RSAVECTORSIZE]; //local vector of skips
static uint32_t m1[RSAVECTORSIZE]; //local vector of messages mod p1 for use in Chinese Remainder Theorem calculation
static uint32_t m2[RSAVECTORSIZE]; //local vector of messages mod p2 for use in Chinese Remainder Theorem calculation
static double randomRSAvector[RSAVECTORSIZE]; //local vector of double precision pseudorandom numbers

//used to force initialization of state if generator not initialized by user
//prior to first call to randomRSA or vrandomRSA
static unsigned char firsttime=1;  

//index used to determine when above local arrays need to be filled or refilled
uint32_t vectorindex=RSAVECTORSIZE; 

//total number of rands returned by generator 
long randsreturned=0;

static uint32_t N_a=11;   //number of well-tested primitive roots mod q

//number of safe primes in [2^31..2^32]
static uint32_t N_p=3060793;  

// number of safe primes > 2^31.5 that we will use, full number is 1768946
static uint32_t N_p1=1768747;  // number of safe primes > 2^31.5 that we will use, full number is 1768946


// avoid safe primes close to 2^32 and 2^31.5 to ensure p1>2^31.5>p2>2^31
static uint32_t N_p1_offset=100; 

//max number of safe primes near q/p1 we will use, ensures |n/q|-1 < 10^-5 
static uint32_t N_p2=59; 

//return one double precision pseudorandom number 
double randomRSA()
{
	if (vectorindex==RSAVECTORSIZE) 
	{
		//fill or refill local array of pseudorandom numbers with length RSAVECTORSIZE
		vrandomRSA0();
	}
	vectorindex++; 
	randsreturned++;
	return randomRSAvector[vectorindex-1];	
}


//return a vector of nrandomRSA double precision pseudorandom numbers 
int vrandomRSA(double * randomRSA, long nrandomRSA)
{
	//#pragma omp parallel for schedule(static)
	for(int i=0;i<nrandomRSA;i++)
	{
		if (vectorindex==RSAVECTORSIZE) 
		{
			//fill or refill local array of pseudorandom numbers with length RSAVECTORSIZE
			vrandomRSA0();
		}
		vectorindex++;
		randomRSA[i]=randomRSAvector[vectorindex-1];
	}
	randsreturned+=nrandomRSA;
	return 0;
}


//update local vectors of messages, skips, and double precision pseudorandom numbers 
//based on non-cryptographic RSA 
int vrandomRSA0()
{
	if (firsttime) 
	{
		firsttime=0;
		randomRSA_init();
	}
	
	//reset index to beginning of the vector
	vectorindex=0;
	
	//update vectors s,m,c,r 
	//s=(a s) mod q
	//m=(m+s) mod n
	//c=m^e mod n
	//r=(double) c / (double) n

//vectorized calculation of skips, messages, ciphertexts, and randoms 
//using Chinese Remainder Theorem.
//this is the vectorized inner most loop 
#ifdef OMP	
	#pragma omp parallel for schedule(static)
#endif
	for (int i=0;i<RSAVECTORSIZE;i++)
	{
		uint64_t s1,s2;
		uint32_t c1,c2,y11,y22,h;
		uint64_t c;

		//calculate s = a s mod q using restricted multiplier a	
		s1=a*(s[i]%q1);
		s2=q2*(s[i]/q1);
		if (s2>s1) s[i]=s1+(q-s2); else s[i]=s1-s2;
		
		// calculate m=(m+s) mod n using Chinese Remainder Theorem
		c1=m1[i]=(m1[i]+s[i])%p1;
		c2=m2[i]=(m2[i]+s[i])%p2;
		
		y11 = (c1 * (uint64_t) c1)%p1;
		y22 = (c2 * (uint64_t) c2)%p2;
		uint32_t e0=e>>1;

		//calculate m^e mod n using Chinese Remainder Theorem
		//assumes e odd exponent greater than or equal to 3 		
    	while(e0>1)
    	{
        	if (e0 & 1) 
        	{
        		c1 = ((uint64_t) c1*y11)%p1;
        		c2 = ((uint64_t) c2*y22)%p2;
        	}
        	y11=(y11*(uint64_t) y11)%p1;
        	y22=(y22*(uint64_t) y22)%p2;
        	e0 >>= 1;
    	}
        c1 = ((uint64_t) c1*y11)%p1;
        c2 = ((uint64_t) c2*y22)%p2;
        
  		//calculate ciphertext c from c1 and c2 using Garner's formula
        if (c2>c1) h=c1+(p1-c2); else h=c1-c2;
        h=((uint64_t) h*p2inversemodp1)%p1;
        c=(uint64_t)h*p2+c2;

		//calculate uniform random double on (0,1], ensuring 1.0 is not returned
        randomRSAvector[i]=(double) c / (double) n;
        if (randomRSAvector[i]==1.0) randomRSAvector[i]=0.999999999999999889;
      
 	}	

	return 0;

}





#define TWO23 8388608
#define MILLION 1000000


//initialize state using time in microseconds since 1/1/1970, mpirank, and mpisize
int randomRSA_init()
{
		timemicrosecondssince1970=getMicrotime();
		int mpirank=0,mpisize=1;
#ifdef MPI
		int mpiinitialized=0;
		MPI_Initialized(&mpiinitialized);
		if (mpiinitialized)
		{
#ifdef DEBUGPRINT
		printf("initializing using MPI\n");
#endif		
			MPI_Comm_size(MPI_COMM_WORLD, &mpisize); 
			MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
		}
		else
		{
#ifdef DEBUGPRINT
		printf("initializing without using MPI\n");
#endif			
		}
#endif		
		randomRSA_init_seed_MPI(timemicrosecondssince1970,mpirank,mpisize);
		return 0;
}

//initialize all instances using assigned value of seed, with mpirank=0 and mpisize=1
int randomRSA_init_seed(uint64_t seed)
{
		randomRSA_init_seed_MPI(seed,0,1);
		return 0;
}


//initialize state using assigned values of seed, mpirank, and mpisize
int randomRSA_init_seed_MPI(uint64_t seed, int mpirank, int mpisize)
{
		//mark generator as initialized
		firsttime=0;
		

		timemicrosecondssince1970=seed;
#ifdef DEBUGPRINT
		printf("initializing randomRSA\n");
		printf("seed=%llu\n",seed);
#endif
		if (mpirank<0) mpirank=0;
		if (mpisize<1) mpisize=1;
		
		//set exponent to small odd number which ensures e coprime to (p1-1)(p2-2) for safe primes p1 and p2	
		e=exponentsaved;  //recommend e=5,9,17

#ifdef DEBUGPRINT		
		printf("mpirank,mpisize=%d %d\n",mpirank,mpisize);
		printf("e=%u\n",e);
#endif		

		//prime used in pseudorandom skip 
		//q=2^63-25 largest prime less than 2^63
		q=9223372036854775783UL; 
		//a=2307085864; //ECuyer,Blounin,Couture restricted primitive root
		
		//use fewer values of p2 close to q/p for smaller number of processes so that
		//n=p1*p2 is closer to q
		//with N_p2=7, |n-q|/q - 1 < 10^-6
		N_p2=59;
		if (mpisize <= 32768) N_p2=47;
		if (mpisize <= 4096) N_p2=23;
		if (mpisize <= 256) N_p2=7;
	

		int sigma=(N_p1 * N_p2 * N_a)/(mpisize+N_p2+N_a);
		while (gcd(sigma,N_p1 * N_p2 * N_a)!=1) sigma--;
		uint64_t beta = powermod((timemicrosecondssince1970/(MILLION/N_p2)+mpirank*sigma)%(N_p1*N_p2),257,N_p1*N_p2*N_a);
		//beta = powermod(beta,17,(N_p1*N_p2));
		
		
		//choose 32-bit restricted primitive root
		a = two63minus25primitiveroots[beta%N_a];
		
		//32-bit parameters so pseudorandom skip calculation can be done
		//using unsigned 64-bit arithmetic
		q1=q/a;
		q2=q%a;
		
#ifdef DEBUGPRINT			
		printf("a,q1,q2=%u %u %u\n",a,q1,q2);	
		printf("time,rank,beta,sigma= %llu %d %d %llu %d\n",timemicrosecondssince1970,mpirank,mpisize,beta,sigma);
#endif		

		//choose 32-bit safe primes p1 and p2 and calculate public key composite modulus n=p1*p2

		//choose safeprime p1 between 2^31.5 and 2^3; p1 unique for each process 
		long safeprime1index=(beta%N_p1+N_p1_offset);
		p1=safeprime(safeprime1index);
		
		//choose safeprime p2 near q/p1
		int safeprime2offset=beta%N_p2 - (N_p2/2);
		if (safeprime2offset>=0) safeprime2offset++; 
		p2=nextsafeprime(q/p1,safeprime2offset);
		
		//64-bit composite public key 
		n=(uint64_t)p1*(uint64_t)p2;
		
		//parameter used in Chinese Remainder Theorem calculation
		//using unsigned 64-bit arithmetic
		p2inversemodp1=modinverse(p2,p1);

#ifdef DEBUGPRINT			
		printf("safeprime1index=%ld\n",safeprime1index);
		printf("p1=%u\n",p1);
		printf("safeprime2offset=%d\n",safeprime2offset);
		printf("p2=%u\n",p2);			
		printf("p1,p2,p2inv,n=%u %u %u %llu\n",p1,p2,p2inversemodp1,n);
#endif
		
		//initialize skip array to prevent skip overlap during any run shorter than about 2^63 random numbers 
		uint64_t delta=((__uint128_t)timemicrosecondssince1970 * TWO23 + (q/mpisize)*mpirank)%(q-1);
		s[0]=powermod(a,delta,q);
		uint64_t delta1=(((q-1)/RSAVECTORSIZE))%(q-1);
		uint64_t b=powermod(a,delta1,q);
#ifdef DEBUGPRINT	
		printf("delta=%llu   s[0]=%llu\n",delta,s[0]);
		printf("delta1=%llu   s1=%llu\n",delta1,b);
#endif
		uint64_t m;
		m = ((__uint128_t)timemicrosecondssince1970)%n;
		m=powermod(m,17,n);
		m1[0]=m%p1;
		m2[0]=m%p2;
#ifdef DEBUGPRINT			
       	printf("0,s[0],m1[0],m2[0],m= %ld %llu %u %u %llu\n",0L,s[0],m1[0],m2[0],m);
#endif
		for (long i=1;i<RSAVECTORSIZE;i++) 
		{	
			s[i]=((__uint128_t) s[i-1] * (__uint128_t) b)%q;
			m = ((__uint128_t)timemicrosecondssince1970*(i+1))%n;
			m=powermod(m,17,n);
			m1[i]=m%p1;
			m2[i]=m%p2;

#ifdef DEBUGPRINT				
       //debug
        if (i<5)
        {
 			uint64_t m,h;
			if (m2[i]>m1[i]) h=m1[i]+(p1-m2[i]); else h=m1[i]-m2[i];
        	h=((uint64_t) h*p2inversemodp1)%p1;
        	m=h*p2+m2[i];
       		printf("i,s[i],m1[i],m2[i],m= %ld %llu %u %u %llu\n",i,s[i],m1[i],m2[i],m);
        }
#endif
		}

#ifdef DEBUGPRINT
		printf("finished initializing randomRSA\n");		
#endif

		randomRSA_validate(10);
		return 0;
}




static int randomRSA_validate(int numvectorcalls)
{
	
	uint64_t *skips=malloc(RSAVECTORSIZE*sizeof(uint64_t));
	uint64_t *messages=malloc(RSAVECTORSIZE*sizeof(uint64_t));
	double *randoms=malloc(RSAVECTORSIZE*sizeof(uint64_t));
	
	int ok=1;
	
	vrandomRSA0();
	for (int iter=0;iter<numvectorcalls;iter++)
	{
		//printf("testing: m1[0],m2[0],s[0],r[0]=%u %u %llu %f\n",m1[0],m2[0],s[0],randomRSAvector[0]);
	for(int i=0;i<RSAVECTORSIZE;i++)
	{
		uint32_t h;
		uint64_t c;

		//set current values of skips and messages
		skips[i]=s[i];
		//use Garner's formula to extract messages from CRT representation
		if (m2[i]>m1[i]) h=m1[i]+(p1-m2[i]); else h=m1[i]-m2[i];
        h=((uint64_t) h*p2inversemodp1)%p1;
        messages[i]=(uint64_t)h*p2+m2[i];

		//calculate the next values using 128-bit arithmetic
		skips[i]=((__uint128_t) a * (__uint128_t) skips[i])%q;
		messages[i]= ((__uint128_t) messages[i] + (__uint128_t) skips[i])%n;
        c=powermod(messages[i],e,n);
        randoms[i]=(double) c / (double) n;
        if(randoms[i]==1.0) randoms[i]=0.999999999999999889;
	}
	//compare to values using 64-bit arithmetic
	vrandomRSA0();
	for(int i=0;i<RSAVECTORSIZE;i++)
	{
		uint32_t h;
		uint64_t message;
		

		if(skips[i] != s[i])
		{
			ok=0;
			//printf("error with skips:i,s[i],skips[i]=%d %llu %llu\n",i,s[i],skips[i]);
		}
		//use Garner's formula to extract messages from CRT representation
		if (m2[i]>m1[i]) h=m1[i]+(p1-m2[i]); else h=m1[i]-m2[i];
        h=((uint64_t) h*p2inversemodp1)%p1;
        message=(uint64_t)h*p2+m2[i];
		if(messages[i] != message)
		{
			ok=0;
			//printf("error with messages:i,m[i],messages[i]=%d %llu %llu\n",i,message,messages[i]);
		}
		if(randoms[i] != randomRSAvector[i])
		{
			ok=0;
			//printf("error with randoms:i,randomRSAvector[i],randoms[i]=%d %f %f\n",i,randomRSAvector[i],randoms[i]);
		}		
	}	
	}
	
	free(skips);
	free(messages);
	free(randoms);

	if (!ok)
	{
		fprintf(stderr,"Error: randomRSA validation failed\n");
		return 1;
	}


#ifdef DEBUGPRINT
		printf("randomRSA validation complete with no errors after %d calls to vrandomRSA0()\n",numvectorcalls);
#endif	
	return 0;
}


//get the current value of the exponent
int randomRSA_get_exponent(uint32_t *exponent)
{
	*exponent=e;
	return 0;
}



int randomRSA_set_exponent(uint32_t exponent)
{
	if (firsttime) 
	{
		firsttime=0;
		randomRSA_init(); //initialize state on first call using time, mpirank and mpisize
	}
	if (exponent%2==1 && exponent >= 3 && exponent <= 257)
	{
		e=exponentsaved=exponent;
#ifdef DEBUGPRINT
		printf("randomRSA exponent reset to %u\n",e);
#endif
		vrandomRSA0(); //reset vector
		return 0;
	}
	fprintf(stderr,"Error: attempt to reset randomRSA exponent to illegal value %u: no changes made.\n",exponent);
	return 1;
}

//get the current values of the safe primes, as well as composite n=p1*p2
int randomRSA_get_primes(uint32_t *prime1, uint32_t *prime2, uint64_t *composite)
{
	*prime1=p1;
	*prime2=p2;
	*composite=n;
}


//set safe primes 2^32>prime1>prime2>2^31, and composite n=p1*p2. If illegal, make no changes.
int randomRSA_set_primes(uint32_t prime1, uint32_t prime2)
{
	if (firsttime) 
	{
		firsttime=0;
		randomRSA_init(); //initialize state on first call using time, mpirank and mpisize
	}
	if (prime1>prime2 && prime2 > 2147483648 && safeprimeq(prime1) && safeprimeq(prime2))
	{
		p1=prime1;
		p2=prime2;
		n=(uint64_t) p1*(uint64_t) p2;
		p2inversemodp1=modinverse(p2,p1);
#ifdef DEBUGPRINT
		printf("randomRSA safe primes reset to p1=%u and p2=%u,\n",p1,p2);
		printf("and n=%llu, and p2inversemodp1=%u\n",n,p2inversemodp1);
#endif
		for (int i=0;i<RSAVECTORSIZE;i++)
		{
			m1[i]=m1[i]%p1;
			m2[i]=m2[i]%p2;
		}
		vrandomRSA0(); //reset vector
		return 0;
	}
	fprintf(stderr,"Error: attempt to set randomRSA safe primes to illegal values: no changes made.\n");	
}


//get q=2^63-25, prime used in skip generator. q can not be changed)
int randomRSA_get_skipprime(uint64_t *skipprime)
{
	*skipprime=q;
	return 0;
}

//return the current value of the primitive root (mod q) (q=2^63-25 is prime used in skip generator and can not be changed)
int randomRSA_get_primitiveroot(uint32_t *primitiveroot)
{
	*primitiveroot=a;
	return 0;
}


//get vector size
int randomRSA_vectorsize()
{
	return RSAVECTORSIZE;
}


//get total number of random numbers returned
long randomRSA_randsreturned()
{
	return randsreturned;
}




//return state of generator to user
int randomRSA_get_state(uint32_t *prime1, uint32_t *prime2, uint32_t *exponent, uint32_t *primitiveroot, uint32_t *index, uint32_t *hash, uint64_t *messages,uint64_t *skips)
{
	if (firsttime) 
	{
		firsttime=0;
		randomRSA_init(); //initialize state on first call using time, mpirank and mpisize
	}
	*hash=0;
	*hash = (*hash + (uint64_t)p1*(uint64_t)p2)%4294967291;
	*hash = (*hash + (uint64_t)e*(uint64_t)a*(uint64_t)(vectorindex+1))%4294967291;

	*prime1=p1;
	*prime2=p2;
	*exponent=e;
	*primitiveroot=a;
	*index=vectorindex;
	for (int i=0;i<RSAVECTORSIZE;i++)
	{
		uint64_t m;
		uint32_t h;
		if (m2[i]>m1[i]) h=m1[i]+(p1-m2[i]); else h=m1[i]-m2[i];
        h=((uint64_t) h*p2inversemodp1)%p1;
        m=(uint64_t)h*p2+m2[i];
		messages[i]=m;
		skips[i]=s[i];
		
		*hash = (*hash+m)%4294967291;
		*hash = (*hash+s[i])%4294967291;
	}
#ifdef DEBUGPRINT
	printf("hash=%llu\n",*hash);
#endif
	return 0;
}



int randomRSA_set_state(uint32_t prime1, uint32_t prime2, uint32_t exponent, uint32_t primitiveroot, uint32_t index, uint32_t hash, uint64_t *messages, uint64_t *skips)
{

	//check to see state is legal in every respect
	unsigned char legalstate=1; 
	legalstate = legalstate && safeprimeq(p1);
	legalstate = legalstate && safeprimeq(p2);
	legalstate = legalstate && (e%2==1);
	legalstate = legalstate && (e<=257);
	unsigned char primitiverootok=0; 
	for (int i=0;i<N_a;i++) if (primitiveroot == two63minus25primitiveroots[i]) primitiverootok=1;
	legalstate = legalstate && primitiverootok;
	legalstate = legalstate && (index<=RSAVECTORSIZE);
	for (int i=0;i<RSAVECTORSIZE;i++) legalstate && (messages[i] < p1*p2);
	for (int i=0;i<RSAVECTORSIZE;i++) legalstate && ((skips[i] > 0) && (skips[i] < q));
	if (!legalstate)
	{
		fprintf(stderr,"Error: attempt to set randomRSA with illegal state. No changes made.\n");
		return 1;
	}
	
	



	//check to see newly calculated value hash0 value agrees with saved hash
	uint64_t hash0=0;
	hash0 = (hash0 + (uint64_t)prime1*(uint64_t)prime2)%4294967291;
	hash0 = (hash0 + (uint64_t)exponent*(uint64_t)primitiveroot*(uint64_t)(index+1))%4294967291;
	for (int i=0;i<RSAVECTORSIZE;i++)
	{
		hash0 = (hash0+messages[i])%4294967291;
		hash0 = (hash0+skips[i])%4294967291;
	}		



	if (hash0!=hash)
	{
		fprintf(stderr,"Error: attempt to set randomRSA state that does not agree with hash value of saved state. No changes made.\n");
		return 1;
	}
	
	
	p1=prime1;
	p2=prime2;
	e=exponentsaved=exponent;
	a=primitiveroot;
	vectorindex=index;
	
	//calculate derived quantities
	n=(uint64_t)p1*(uint64_t)p2;
	p2inversemodp1=modinverse(p2,p1);
	q1=q/a;
	q2=q%a;

	for (int i=0;i<RSAVECTORSIZE;i++)
	{
		uint64_t c;
		m1[i]=messages[i]%p1;
		m2[i]=messages[i]%p2;
		s[i]=skips[i];	
		c=powermod(messages[i],e,n);
        randomRSAvector[i]=(double) c / (double) n;
        if (randomRSAvector[i]==1.0) randomRSAvector[i]=0.999999999999999889;		
	}	
	
#ifdef DEBUGPRINT	
	printf("resetting state\n");
	printf("prime1=%u\n",prime1);
	printf("prime2=%u\n",prime2);
	printf("exponent=%u\n",exponent);
	printf("primitiveroot=%u\n",primitiveroot);
	printf("index=%u\n",index);
	printf("hash=%llu\n",hash);
	printf("messages[0-2]=%llu %llu %llu\n",messages[0],messages[1],messages[2]);
	printf("skips[0-2]=%llu %llu %llu\n",skips[0],skips[1],skips[2]);
#endif

	return 0;

}

//*****************************************************************

//private functions

//Returns the current time in microseconds since 0:00:00 1/01/1970
static uint64_t getMicrotime()
{
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (uint32_t)1e6 + currentTime.tv_usec;
}


// primality test for n<2^32
static int primeq(uint64_t n)
{

    uint64_t prime[12]={2,3,5,7,11,13,17,19,23,29,31,37};


    int ntests=12;//for n<2^64 use ntests=12
    if (n<4294967296) ntests=5;

    for (int i=0;i<12;i++)
    {
        if (n==prime[i]) return 1;
        if (n%prime[i]==0) return 0;
    }


    uint64_t n1,q,k,x,y,s;

    n1=n-1;
    q=n-1;
    k=0;
    while(q%2==0)
    {
        k++;
        q /= 2;
    }

    s=1;
    for (int repeat=0;repeat<ntests;repeat++)
    {
        //s=multiplymod(a,s,p);
        x=prime[repeat];  //  first ntests primes
        y=powermod32(x,q,n);
        int ok=0;
        //cout <<"n,k,q,x,y="<<n<<" "<<k<<" "<<q<<" "<<x<<" "<<y<<endl;
        for (int j=0;j<k;j++)
        {
            //cout <<"n,k,q,j,x,y="<<n<<" "<<k<<" "<<q<<" "<<j<<" "<<x<<" "<<y<<endl;
            if ( (j==0 && y==1) || (y==n1))
            {
                ok = 1;
                break;
            }
            if (j>0 && y==1) return 0;  //not prime
            y=(y*y)%n;
            //y=multiplymod(y,y,n);
        }
        if (ok==0) return 0; // not prime
    }
    return 1; //prime
}

// tests to see if n is a strong prime for n<2^32
static int safeprimeq(uint64_t n)
{
    uint64_t  n12;

    if(n%2==0) return 0;

    if (n < 84)
    {
        if (n==5) return 1;
        if (n==7) return 1;
        if (n==11) return 1;
        if (n==23) return 1;
        if (n==47) return 1;
        if (n==59) return 1;
        if (n==83) return 1;
        return 0;
    }

	if(n%12!=11) return 0;
	
    n12=(n-1)/2;

   	uint64_t prime[12]={2,3,5,7,11,13,17,19,23,29,31,37};

    for (int i=0;i<12;i++) if (n%prime[i]==0 || n12%prime[i]==0) return 0;

    if (primeq(n)==0) return 0;

    if (primeq(n12)==0) return 0;

    return 1;

}

// selects the next prime above n if offset = +1, jth prime above n if offset = +j
// selects the next prime below n if offset = -1, jth prime below n if offset = -j
static uint64_t nextprime(uint64_t n, long offset)
{
    int updown,nprimes,pq;
    if (offset == 0) return 0;
    if (offset > 0) updown=+1; else updown=-1;

    offset=labs(offset);
    nprimes=0;
    do
    {
        n += updown;
        if (n<2) return 0;
        pq=primeq(n);
        if (pq) nprimes++;
    }
    while (!pq || nprimes < offset);

    return n;
}

// selects the next strong prime above n if offset = +1, jth strong prime above n if offset = +j
// selects the next strong prime below n if offset = -1, jth stong prime below n if offset = -j
static uint64_t nextsafeprime(uint64_t n, long offset)
{
    int updown,nprimes,pq;
    if (offset == 0) return 0;
    if (offset > 0) updown=+1; else updown=-1;

    offset=labs(offset);
    nprimes=0;
    do
    {
        n += updown;
        if ( n < 2) return 0;
        pq=safeprimeq(n);
        if (pq) nprimes++;
    }
    while (!pq || nprimes < offset);

    return n;
}


//returns one of 3060794 strongprimes between 2^31 and 2^32 starting near 2^32
// with 0 <= i < 3060794
// safeprime(0) = 4294967087  is the largest strong prime  below 2^32
// safeprime(3060793) = 2147483783  is the smallest strong prime above 2^31
// safeprime(1768946)=3037000943 is smallest safeprime greater than floor(2^31 sqrt(2))
// for use in choosing p1 in RSA 
static uint32_t safeprime(long i)
{
    uint64_t safeprimes[31]={4294967087, 4222725707, 4150883063, 4079008547, 4006857719, 3935091059, 3863505827, 3792054359, 3720601559, 3649365899, 3578380259, 3507339959, 3436690847, 3366364187, 3295757987, 3225327167, 3155044799, 3085122659, 3015377447, 2945577263, 2876115563, 2806385207, 2736992087, 2667700019, 2598759287, 2530175939, 2461716743, 2393244167, 2324812307, 2256486719, 2188727483};

    if (i<0) return 0;
    if (i>=3060794) return 0;


    uint64_t n=safeprimes[i/100000];
    int j=i%100000;
    if (j > 0) n=nextsafeprime(n, -j);
    return n;
}


// x^e mod n
// assumes n<2^64
static uint64_t powermod(uint64_t x, uint64_t e, uint64_t n)
{

   uint64_t y=x%n,z=1;

    if (n < 4294967296UL)
    {
    	while(e>0)
    	{
        	if (e & 1) z = (z*y)%n;
        	if (e == 1) return z;
        	y=(y*y)%n;
        	e >>= 1;
    	}
    }
    else
    {
     	while(e>0)
    	{
        	if (e & 1) z=((__uint128_t) z * (__uint128_t) y) % n;  //z = (z*y)%n;
        	if (e == 1) return z;
        	y=((__uint128_t) y * (__uint128_t) y) % n;  //y=(y*y)%n;
        	e >>= 1;
    	}
    }

    return z;
}

// returns x^e mod n for 32-bit numbers
static inline uint32_t powermod32(uint32_t x, uint32_t e, uint32_t n)
{

    uint64_t y,z;

    y=x%n;
    z=1;

     while(e>0)
    {
        if (e & 1) z = (z*y)%n;
        if (e == 1) return z;
        y=(y*y)%n;
        e >>= 1;
    }


    return z;
}





// greatest common denominator of u and v bv Euclid's algorithm (check for u==0 or v==0)??
static uint64_t gcd(unsigned long u, unsigned long v)
{
 	while (v != 0) 
 	{
 		uint64_t temp = v;
 		v = u % v;
 		u = temp;
    }
    return u;
}

static uint64_t modinverse(uint64_t u, uint64_t v)
{
    uint64_t inv, u1, u3, v1, v3, t1, t3, q;
    int iter;
    // Step X1. Initialise 
    u1 = 1;
    u3 = u;
    v1 = 0;
    v3 = v;
    // Remember odd/even iterations 
    iter = 1;
    /* Step X2. Loop while v3 != 0 */
    while (v3 != 0)
    {
        // Step X3. Divide and "Subtract" 
        q = u3 / v3;
        t3 = u3 % v3;
        t1 = u1 + q * v1;
        // Swap 
        u1 = v1; v1 = t1; u3 = v3; v3 = t3;
        iter = -iter;
    }
    /* Make sure u3 = gcd(u,v) == 1 */
    if (u3 != 1)
        return 0;   /* Error: No inverse exists */
    /* Ensure a positive result */
    if (iter < 0)
        inv = v - u1;
    else
        inv = u1;
    return inv;
}

