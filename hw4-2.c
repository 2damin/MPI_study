#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#define MaxNum 10000

int main(int argc, char** argv)
{
    int rank, size;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int i,j;
    int sum =0, TotalSum =0;
    int start = MaxNum*rank/size+1;
    int end = MaxNum*(rank+1)/size+1;

    for(i=start;i<end;i++)
    {
        sum += i;
    }

    MPI_Reduce(&sum,&TotalSum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    if(rank==0)
    {
        printf("Rank%d -- Total sum : %d\n",rank,TotalSum);
    }
    else if(rank != 0)
    {
        printf("Rank%d -- Total sum : %d\n",rank,TotalSum);
    }

    MPI_Finalize();

}
