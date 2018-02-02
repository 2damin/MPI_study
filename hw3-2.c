#include <stdio.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char** argv)
{
    int rank,size;
    MPI_Init(NULL,NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int TotalSum = 0;
    int sum = 0;
    int MaxNum = 10000;
    int i,j;
    clock_t start, end;

    start = clock();

    if(rank == size-1)
    {
        if(rank != 0)
        {
        MPI_Recv(&sum,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        for(i =(MaxNum*rank)/size+1; i<MaxNum*(rank+1)/size+1; i++)
        {
            sum += i;
        }
        end = clock();

        TotalSum = sum;
        printf("rank %d sum: %d\n\n",rank,sum);
        printf("---------------------------------------------------\n");
        printf("Total sum is %d\n\n", TotalSum);
        printf("Total running time is %lf \n", (end-start)/(double)1000);
    }
    else if(rank != size-1)
    {
        if(rank != 0)
        {
        MPI_Recv(&sum,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }

        for(i=(MaxNum*rank)/size+1; i<MaxNum*(rank+1)/size+1; i++)
        {
            sum += i;
        }
        printf("rank %d sum: %d\n\n",rank,sum);

        MPI_Send(&sum,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);

    }

    MPI_Finalize();
    return 0;
}
