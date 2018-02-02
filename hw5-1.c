#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int main(int argc, char** argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i,j,N;
    if(rank == 0)
    {
        printf("what is N?\n");
        scanf("%d",&N);
        printf("\nN is %d\n",N);
    }
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);

    clock_t starttime, endtime;
    double totaltime;
    double sum = 0.0, TotalSum = 0.0;
    int sendsize = N/size, recvsize = N/size;
    double *senddata = (double*)malloc(sizeof(double)*N);
    double *sum_array = (double*)malloc(sizeof(double)*size);
    double *recvdata = (double*)malloc(sizeof(double)*recvsize);
    if(rank == 0)
    {
        starttime=clock();

        for(i=0;i<N;i++)
        {
            senddata[i]=sqrt((double)(i+1))/(double)N;
        }
    }

    MPI_Scatter(senddata,sendsize,MPI_DOUBLE,recvdata,
            recvsize,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for(j=0;j<recvsize;j++)
    {
        sum += recvdata[j];
    }
    printf("partial sum in rank%d: %f\n\n",rank,sum);

    MPI_Gather(&sum,1,MPI_DOUBLE,sum_array,1,MPI_DOUBLE,0,
            MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0)
    {
        for(i=0;i<size;i++)
        {
            TotalSum += sum_array[i];
        }
        endtime =clock();
        totaltime = ((double)(endtime-starttime))/CLOCKS_PER_SEC;

        printf("--------------------------------\n");
        printf("\ntotal sum : %f \n",TotalSum);
        printf("running time: %f\n",totaltime);
    }

    free(senddata);
    free(recvdata);
    free(sum_array);

    MPI_Finalize();
}
