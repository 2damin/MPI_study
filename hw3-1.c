#include <stdio.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char** argv)
{
    int rank, size;
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int MaxNum = 10000;
    int sum = 0;
    int TotalSum = 0;
    int i,j;
    clock_t start, end;
    start = clock();

    if(rank != 0)
    {
        for(i=MaxNum*(rank)/size+1; i < MaxNum*(rank+1)/size+1; i++)
        {
            sum += i;
        }
        MPI_Send(&sum,1,MPI_INT,0,0,MPI_COMM_WORLD);
        sum = 0;
    }
    else if(rank == 0)
    {
        for(i=MaxNum*(rank)/size+1; i < MaxNum*(rank+1)/size+1; i++)
        {
            sum += i;
        }
        printf("Rank0 calculated sum '%d' by itself \n\n",sum);

        TotalSum = sum;

        for(j=1; j<size; j++)
        {
            MPI_Recv(&sum,1,MPI_INT,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            TotalSum += sum;

            printf("Rank0 received sum '%d' from rank%d \n\n",sum, j);
        }

        end = clock();
        printf("--------------------------------------------------------------------- \n");
        printf("Total Sum is %d\n",TotalSum);
        printf("Total running time is %lf seconds \n",(end-start)/(double)1000);
    }
    MPI_Finalize();
    return 0;
}
