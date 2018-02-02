#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#define MaxNum 10000

int main(int argc, char** argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i,j;
    int sum = 0, TotalSum = 0;
    int sendsize = MaxNum/size, recvsize = MaxNum/size;
    int senddata[MaxNum];
    int *sum_array = (int*)malloc(sizeof(int)*size);
    int *recvdata = (int*)malloc(sizeof(int)*recvsize);
    if(rank == 0)
    {
        for(i=0;i<MaxNum;i++)
        {
            senddata[i]=i+1;
        }
    }

    MPI_Scatter(senddata,sendsize,MPI_INT,recvdata,
            recvsize,MPI_INT,0,MPI_COMM_WORLD);

    for(j=0;j<recvsize;j++)
    {
        sum += recvdata[j];
    }
    printf("partial sum in rank%d: %d\n\n",rank,sum);

    MPI_Gather(&sum,1,MPI_INT,sum_array,1,MPI_INT,0,
            MPI_COMM_WORLD);

    if(rank == 0)
    {
        for(i=0;i<size;i++)
        {
            TotalSum += sum_array[i];
        }
        printf("--------------------------------\n");
        printf("\ntotal sum : %d \n",TotalSum);
    }

    free(recvdata);
    free(sum_array);

    MPI_Finalize();
}
