#include <iostream>
#include <cmath>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{ 
	int inneriteration = 30;
	int outeriteration = 300;
	double startime, finishtime;
	int N = 64;					// Number of intervals in one dimension .			
	startime = MPI_Wtime();
	double h = 1.0/N; 			// mesh size
	double error;				//to calculate Euclidean norm
	double errorsquare;			//to calculate Euclidean norm
	double sumoferrorsquare = 0.0;//to calculate Euclidean norm
	double cmidpoint;			//correct midpoint
	double midpointerror;		//error in midpoint
	int messagerows = 2;			//Every message sent has two rows of values
	int messagesize;				//no. of elements in message
	MPI::Status status;			// status to be used in mpi  messages
	
	MPI::Init(argc, argv);		//Initializing mpi
	int P = MPI::COMM_WORLD.Get_size(); 	//get the total no. of processors
	int rank = MPI::COMM_WORLD.Get_rank();	//get the rank of the processor in which this code is running

	int *borders = new int[P];		//the index of the array corresponds to the processor and the value corresponds to the no. of borders the processor has to 										  deal with
	borders[0] = 1;			
	borders[P-1] = 1;
	for (int i = 1 ; i < P-1; i++)
		borders[i] = 2;
		
/*------------------------------------------------------------------------------CODE FOR 1st PROCESSOR------------------------------------------------------------------------------------------------------*/
	if(rank == 0)					
	{
		int noderows = (N/P) + 1 + borders[rank];		//No. of  rows the processor has to deal with
		int nodearraysize = noderows*(N+1);			//No . of elements = row*(N+1) 	
		int *globaltolocal = new int[N+1];				//An array that is used to map the local array index to the global array index (the final combined result). 													  The index of the array corresponds to the global address of the rows that the processor contains and 														  the value is the local row number in processor's local array
		int dfgh = 0;
		for(int akak = 0; akak < noderows; akak++)		//initializing global array for processor 0
		{
			globaltolocal[akak] = dfgh;
			dfgh++;
		}

		double *cu = new double[(N+1)*(N+1)];		//cu is array with the correct solution to compare with our calculated one
		double *u   = new double[nodearraysize];		//u is the array used for iterations
		double *calcu = new double[(N+1)*(N+1)];		// calcu is the array where the values from other processors are combined to form final result
		
		messagesize = messagerows*(N+1);			//Size of each messages to be sent or received
		double *sendarray = new double[messagesize];		//used for comm
		double *receivearray = new double[messagesize];
		double *globalreceivearray = new double[(N/P +1)*(N+1)];		//array where the final results is received using mpirecv
	
		cmidpoint = 2*((floor((N+1)/2))*(floor((N+1)/2)))*h*h;	//correct midpoint

		for(int j = 0; j < N+1 ; j++)						//initalizing the values for the array with correct solution
		{
			for ( int i = 0; i < N+1; i++)
			{
				cu[j*(N+1)+i] = i*i*h*h + j*j*h*h ;
			}
		}
	
		for (int i = 0; i < nodearraysize; i++)		//initializing the array elements to zero	
			u[i] = 0.0;
	
		for(int i = 0; i < N+1; i++)					//calculating boundary values
			u[i] = i*i*h*h;
		
		for ( int i = 1; i < noderows ; i++)			//calculating boundary values
		{
			u[ (i+1)*N + i ] = i*i*h*h+1.0;
			u[ i*N + i ] = i*i*h*h;
		}
		
		int fl;

		for(int outer = 0; outer < outeriteration; outer++)			//outer Gauss Seidel iteration
		{
			for(int inner = 0; inner < inneriteration; inner++)		//inner iteration: after each inner iteration, the processors exchange the boundary values
			{
				for(int j = 1; j < noderows-1 ; j++)		//to traverse the rows
				{
					for ( int i = 1; i < N; i++)			//to traverse the elements within a row
					{
						u[j*(N+1) + i] = 0.25* (-h*h*4.0 + u[j*(N+1) +i-1] + u[j*(N+1) +i+1] + u[(j-1)*(N+1) +i] + u[(j+1)*(N+1)+i]);
					}
				}
			}
			
			
			for( int jj=0; jj < messagerows; jj++)		//loading the send array with the boundary values. Will be sent to 2nd processor.
			{
				for( int ii=0; ii < N+1; ii++)
				{
					sendarray[ii+jj*(N+1)] = u[ii + ((noderows-1)-messagerows+jj)*(N+1)];
				}
			}
			MPI::COMM_WORLD.Send(sendarray, messagesize, MPI_DOUBLE, rank+1,1);
					
			MPI::COMM_WORLD.Recv(receivearray, messagesize, MPI_DOUBLE, rank+1,1);	//receiving values from the second processor.
			for( int ii=0; ii < N+1; ii++)													//calculating average for the shared grid point
			{
				u[ii + ((noderows-1)-1+0)*(N+1)] = (receivearray[ii+0*(N+1)]+u[ii + ((noderows-1)-1+0)*(N+1)])/2;
			}
			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + ((noderows-1)-1+1)*(N+1)] = receivearray[ii+1*(N+1)];
			}
			
		}
		
		/*------------------------------------Now starts the gathering of results from all other processors ---------------------------------------*/
		int globaladdress;
		int localaddress;
		int globalsendsize = ((N/P)+1)*(N+1);
		
		for(int iakak = 0; iakak < (N/P+1); iakak++)			//copying the values from iteration array to the calcu where final values are available
		{
			globaladdress = N/P*rank+iakak;	
			localaddress = globaltolocal[globaladdress];
			for(int jakak = 0; jakak < N+1; jakak++)
			{
				calcu[jakak+iakak*(N+1)] = u[jakak + localaddress*(N+1)];
			}
		}
		
		int somecounter = 0;
		for(int rrr = 1; rrr < P; rrr++)				// loop to gather the calculated values from each processor. 
		{
			somecounter = 0;
			MPI::COMM_WORLD.Recv(globalreceivearray, globalsendsize, MPI_DOUBLE, rrr, 2);		//diff message tag of 2 has been used for this gather
			for(int j = (N/P)*rrr; j <= (N/P)*(rrr+1); j++)
			{
				for (int i = 0; i < N+1 ; i++)
				{
					calcu[i + j*(N+1)] = globalreceivearray[somecounter];
					somecounter++;
				}
			}
		}
			
		cout <<"\n";
		cout <<"\n";
		//cout<<"\nThe caluclated values are "<<"\n";
		
		finishtime=MPI_Wtime();
		
		cout<<"Elapsed time is "<<finishtime-startime<<"\n";
	
		for(int j = 0; j < N+1 ; j++)
		{
			for ( int i = 0; i < N+1; i++)
			{
				//cout<<calcu[j*(N+1)+i]<<"	";
				error = cu[j*(N+1)+i] - calcu[j*(N+1)+i];
				errorsquare = error*error;
				sumoferrorsquare = errorsquare + sumoferrorsquare;
			}
			//cout<<"\n";
		}
		
		cout <<"\n";
		cout <<"\n";
		cout<<"\nThe euclidean norm is  "<<sqrt(sumoferrorsquare)<<"\n";
		
		fl = floor((N+1)/2);
		cout<<"The midpoint error is "<<cu[(N+1)*fl + fl] - cmidpoint<< " for a calc midpoint of  "<< cu[(N+1)*fl + fl] <<"\n";
		
		delete[] u;							//free up all the memories allocated in the heap
		delete[] cu;
		delete[] sendarray;
		delete[] receivearray;
		delete[] globaltolocal;
		delete[] globalreceivearray;
		delete[] calcu;
	}
	
/*---------------------------------------------------------------------------CODE FOR LAST PROCESSOR----------------------------------------------------------------------------------------------*/	
	if(rank == P-1 && P>1)
	{
		int noderows = (N/P) + 1 + borders[rank];	//No. of  rows the processor has to deal with
		int nodearraysize = noderows*(N+1);		//No . of elements = row*(N+1) 	

		double *u   = new double[nodearraysize];    //Array where calculations are carried out
		int *globaltolocal = new int[N+1];			//An array that is used to map the local array index to the global array index (the final combined result). 												The index of the array corresponds to the global address of the rows that the processor contains and 													the value is the local row number in processor's local array
		
		int dfgh = 0;
		for(int akak = ((N/P)*rank-1); akak < (((N/P)*rank-1)+noderows); akak++)			//initializing global array for processor 0
		{
			globaltolocal[akak] = dfgh;
			dfgh++;
		}

		messagesize = messagerows*(N+1);								//Size of each messages to be sent or received
		double *sendarray = new double[messagesize];					//used for comm
		double *receivearray = new double[messagesize];
		double *globalsend = new double[(N/P + 1)*(N+1)];					//array to send the calculated values after outer iteration to the 1st processor
	
		for (int i = 0; i < nodearraysize; i++)								//initializing the array elements to zero	
			u[i] = 0.0;
	
		for(int i = 0; i < N+1; i++)											//calculating boundary values
			u[(noderows-1)*(N+1)+i] = i*i*h*h+1.0;
			
		int bakbak = 0;
		for ( int i = (((N/P)*rank)-1); i < noderows+(((N/P)*rank)-1) ; i++)		//calculating boundary values
		{
			u[ (bakbak+1)*N + bakbak ] = i*i*h*h+1.0;
			u[ bakbak*N + bakbak ] = i*i*h*h;
			bakbak++;
		}
		
		int fl;

		for(int outer = 0; outer < outeriteration; outer++)							//outer Gauss Seidel iteration
		{
			for(int inner = 0; inner < inneriteration; inner++)		//inner iteration: after each inner iteration, the processors exchange the boundary values
			{
				for(int j = 1; j < noderows-1 ; j++)						//to traverse the rows
				{
					for ( int i = 1; i < N; i++)							//to traverse the elements within a row
					{
						u[j*(N+1) + i] = 0.25* (-h*h*4.0 + u[j*(N+1) +i-1] + u[j*(N+1) +i+1] + u[(j-1)*(N+1) +i] + u[(j+1)*(N+1)+i]);
					}
				}
			}
			
			
			for( int jj=0; jj < messagerows; jj++)
			{
				for( int ii=0; ii < N+1; ii++)
				{
					sendarray[ii+jj*(N+1)] = u[ii + (1+jj)*(N+1)];
				}
			}
			MPI::COMM_WORLD.Send(sendarray, messagesize, MPI_DOUBLE, rank-1,1);
			
			
			MPI::COMM_WORLD.Recv(receivearray, messagesize, MPI_DOUBLE, rank-1,1);
			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + (1)*(N+1)] = (receivearray[ii+1*(N+1)]+u[ii + (1)*(N+1)])/2;
			}
			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + (0)*(N+1)] = receivearray[ii+0*(N+1)];
			}

		}
		
		int globaladdress;
		int localaddress;
		int globalsendsize = (N/P+1)*(N+1);
		
		/*--------------------------------------// sending the calculated values to 1st processor  after  all Gauss Seidel iterations//---------------------------------*/
		for(int iakak = 0; iakak < (N/P+1); iakak++)				
		{
			globaladdress = N/P*rank+iakak;
			localaddress = globaltolocal[globaladdress];
			for(int jakak = 0; jakak < N+1; jakak++)
			{
				globalsend[jakak+iakak*(N+1)] = u[jakak + localaddress*(N+1)];
			}
		}
		
		MPI::COMM_WORLD.Send(globalsend, globalsendsize, MPI_DOUBLE, 0,2);
		
		delete[] u;
		delete[] sendarray;
		delete[] receivearray;
		delete[] globaltolocal;
		delete[] globalsend;
	}


/*---------------------------------------------------------------------------CODE FOR THE MIDDLE PROCESSORS-------------------------------------------------------------------------*/

	if(0 < rank && (P-1) > rank && P>2)
	{
		int noderows = (N/P) + 1 + borders[rank];			//No. of  rows the processor has to deal with
		int nodearraysize = noderows*(N+1);				//No . of elements = row*(N+1) 

		double *u   = new double[nodearraysize];
		int *globaltolocal = new int[N+1];			//An array that is used to map the local array index to the global array index (the final combined result). 											The index of the array corresponds to the global address of the rows that the processor contains and 											  the value is the local row number in processor's local array
		
		int dfgh = 0;
		for(int akak = ((N/P)*rank-1); akak < (((N/P)*rank-1)+noderows); akak++)
		{
			globaltolocal[akak] = dfgh;
			dfgh++;
		}

		messagesize = messagerows*(N+1);
		double *sendarray = new double[messagesize];			// array to send messages to prev processor
		double *receivearray = new double[messagesize];			// array to receive messages from prev processor
		double *rsendarray = new double[messagesize]; 			// array to send messages to next processor
		double *rreceivearray = new double[messagesize];			// array to receive messages from next processor
		double *globalsend = new double[(N/P + 1)*(N+1)];			// array to send calculated values after all iterations to the 1st processor

		for (int i = 0; i < nodearraysize; i++)
			u[i] = 0.0;
		
		int bakbak=0;
		for ( int i = (((N/P)*rank)-1); i < noderows+(((N/P)*rank)-1) ; i++)		//boundary conditions
		{
			u[ (bakbak+1)*N + bakbak ] = i*i*h*h+1.0;
			u[ bakbak*N + bakbak ] = i*i*h*h;
			bakbak++;
		}

		int fl;
		for(int outer = 0; outer < outeriteration; outer++)			//outer Gauss Seidel iteration
		{
			for(int inner = 0; inner < inneriteration; inner++)		//inner iteration: after each inner iteration, the processors exchange the boundary values
			{
				for(int j = 1; j < noderows-1 ; j++)		//to traverse the rows
				{
					for ( int i = 1; i < N; i++)			//to traverse the elements within a row
					{
						u[j*(N+1) + i] = 0.25* (-h*h*4.0 + u[j*(N+1) +i-1] + u[j*(N+1) +i+1] + u[(j-1)*(N+1) +i] + u[(j+1)*(N+1)+i]);
					}
				}
			}
			
			/*---------------------------------------part to send messages to neighboring processor after each inner iteration loop exit--------------------------*/
				
			for( int jj=0; jj < messagerows; jj++)
			{
				for( int ii=0; ii < N+1; ii++)
				{
					sendarray[ii+jj*(N+1)] = u[ii + (1+jj)*(N+1)];
				}
			}
			MPI::COMM_WORLD.Send(sendarray, messagesize, MPI_DOUBLE, rank-1,1);

		
			for( int jj=0; jj < messagerows; jj++)
			{
				for( int ii=0; ii < N+1; ii++)
				{
					rsendarray[ii+jj*(N+1)] = u[ii + ((noderows-1)-messagerows+jj)*(N+1)];
				}
			}
			MPI::COMM_WORLD.Send(rsendarray, messagesize, MPI_DOUBLE, rank+1,1);
			
			/*---------------------------------------part to receive messages to neighboring processor after each inner iteration loop exit--------------------------*/			
			
			MPI::COMM_WORLD.Recv(receivearray, messagesize, MPI_DOUBLE, rank-1,1);		
			
			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + (1)*(N+1)] = (receivearray[ii+1*(N+1)]+u[ii + (1)*(N+1)])/2;
			}

			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + (0)*(N+1)] = receivearray[ii+0*(N+1)];
			}
			
			MPI::COMM_WORLD.Recv(rreceivearray, messagesize, MPI_DOUBLE, rank+1,1);		
			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + ((noderows-1)-1+0)*(N+1)] = (rreceivearray[ii+0*(N+1)]+u[ii + ((noderows-1)-1+0)*(N+1)])/2;
			}
			for( int ii=0; ii < N+1; ii++)
			{
				u[ii + ((noderows-1)-1+1)*(N+1)] = rreceivearray[ii+1*(N+1)];
			}
		}
		
		int globaladdress;
		int localaddress;
		int globalsendsize = (N/P+1)*(N+1);
		
					/*---------------------------------------part to send calculated values to 1st processor after exit of outer iteration loop--------------------------*/
		
		for(int iakak = 0; iakak < (N/P+1); iakak++)			
		{
			globaladdress = N/P*rank+iakak;
			localaddress = globaltolocal[globaladdress];
			for(int jakak = 0; jakak < N+1; jakak++)
			{
				globalsend[jakak+iakak*(N+1)] = u[jakak + localaddress*(N+1)];
			}
		}	
		MPI::COMM_WORLD.Send(globalsend, globalsendsize, MPI_DOUBLE, 0,2);
		
		delete[] u;
		delete[] sendarray;
		delete[] receivearray;
		delete[] rsendarray;
		delete[] rreceivearray;
		delete[] globaltolocal;
		delete[] globalsend;
	}
	
	delete[] borders;
	MPI::Finalize();
	return 0;
}
	
	
