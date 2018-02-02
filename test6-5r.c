#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include<time.h>
#include<stdlib.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    int size,rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int i,j,N,partialnum,myrowsize,realrowsize;
    clock_t start, end;
    double totaltime =0.0;

    if(rank==0)
    {
        printf("What is N? (N*N matrix)\n");
        scanf("%d",&N);
        printf("\nN is %d\n",N);
    }
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);

    int *A=NULL;                                          /* A matrix(N*N) */
    int *B=(int*)malloc(sizeof(int)*N);                   /* B matrix(N)   */
    int *C=NULL;                                          /* C matrix(N)   */
    int *recvcounts=(int*)malloc(sizeof(int)*size);       /* the number of receving data of each rank (gatherv) */
    int *senddisp=(int*)malloc(sizeof(int)*size);         /* the displacement of data size at each rank (scatterv) */
    int *recvdisp=(int*)malloc(sizeof(int)*size);         /* the displacement of data size at each rank (gatherv) */
    int *sendcounts=(int*)malloc(sizeof(int)*size);       /* the number of sending data of each rank (scatterv) */

    myrowsize=(int)(N/size);

    senddisp[0]=0;
    recvdisp[0]=0;
    for(i=0;i<size;i++)
    {
        if(i == size-1)
        {
            sendcounts[i]=N*N-(sendcounts[0]*(size-1));
            recvcounts[i]=N-myrowsize*(size-1);
        }
        else if(i != size-1)
        {
            sendcounts[i]=N*myrowsize;
            recvcounts[i]=myrowsize;
            senddisp[i+1]=senddisp[i]+sendcounts[i];
            recvdisp[i+1]=recvdisp[i]+myrowsize;
        }
    }

    partialnum=sendcounts[rank];                         /*the number of data at each rank*/
    realrowsize=partialnum/N;
    int *Arecv=(int*)malloc(sizeof(int)*partialnum);
    int *result=(int*)malloc(sizeof(int)*realrowsize);

    if(rank==0)
    {
        A=(int*)malloc(sizeof(int)*N*N);
        C=(int*)malloc(sizeof(int)*N);

        //printf("A matrix------------------B matrix =\n \n");
        for(j=0;j<N;j++)
        {
            for(i=0;i<N;i++)
            {
                A[i+j*N]=rand()%100;
               // printf("%d ",A[i+j*N]);
            }
            B[j]= rand()%100;
           // printf("-----%d\n",B[j]);
        }
        start = clock();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(B,N,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Scatterv(A,sendcounts,senddisp,MPI_INT,Arecv,partialnum,MPI_INT,0,MPI_COMM_WORLD);

    int k=-1;
    for(i=0;i<realrowsize;i++)
    {
        result[i]=0;
    }

    for(i=0;i<partialnum;i++)
    {
        if(i%N==0)
        {
            k=k+1;
        }
        result[k]= result[k] + Arecv[i]*B[i-N*k];
    }

    MPI_Gatherv(result,realrowsize,MPI_INT,C,recvcounts,recvdisp,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

   // free(result);
   // free(Arecv);
   // free(Brecv);

    if(rank==0)
    {
        end = clock();
       // printf("\n\nC MATRIX :\n\n");
       // for(i=0;i<N;i++)
       // {
       //     printf("%d\n",C[i]);
       // }

        totaltime = ((double)(end-start))/CLOCKS_PER_SEC;
        printf("---------------------------\n\n");
        printf("total time : %f \n",totaltime);
        free(A);
        free(C);
    }

    free(sendcounts);
    free(recvcounts);
    free(senddisp);
    free(recvdisp);
    free(B);
    free(result);
    free(Arecv);
    

    MPI_Finalize();
    return 0;
}
