   /*******************************/
   /*
   /*  Project Name: lab2
   /*  Description: read in number determine factors/prime
   /*  File names: if multiple .c files exist lab2.c
   /*  Date: Project due date
   /*  Lab Section: 0701
   /*  Programmer:Reinhart, Tony 
   /*  email: aereinha@purdue.edu
   /*
   /*******************************/


#include <stdio.h>



main()
{
int t = 0;     //t is a counter that finds if the number
		// is prime or what factors  until the given number is reached. 
int k = 0;     //k counts how many facotrs there are.
int a = 0;     //shows what number is a factor
int N = 0;	//the user entered number
int b = 0;     //whole number entered
int c = 0;	//number entered

printf("Please enter a number between 1 and 100:");

N = getchar() - '0';
       while(N != EOF - '0' && N != '\n'-'0')
	 {
		if( N >= 0 && N <=9)
		{
			if(t >= 1)
			c = b*10 + N;

		  b = N;
       		   t = t++;
		if(t == 1)
			c = b;
		}
	N = getchar() - '0';
	} 
	if(c < 1 || c > 100) 
	{
		printf("Invalid Input"); 
		return(0);
	}
	for(t=1;t < c; t++)
	{
		a = c % t;
		if(a == 0)
		{
			printf("%d ", t);
			k++;
		}
	}

printf("%d\n", t);

	if(c == t && k==1)
		printf("\nprime\n");
return(0);
}
