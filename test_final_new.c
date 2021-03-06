#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include<time.h>
#include<stdlib.h>
#include<malloc.h>

#pragma warning(disable:4996)

#define pi 3.14159265358979323846264338327950288
#define np 50                   /*node number*/
#define del_x 2.0/(np-1.0)
#define del_y 2.0/(np-1.0)


double exact(int i, int j)                                             /*exact solution*/
{
	return cos(pi*i*del_x)*sin(pi + j*del_y);
}

double f(int i, int j)                                                    /*Poisson eq*/
{
	return cos(pi*i*del_x)*sin(j*del_y) + pi*pi*cos(pi*i*del_x)*sin(j*del_y);
}


/*---------------------------------------------------------------------------------------------------------------------*/



int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    int size,rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    int i,j,partialnum,myrowsize,realrowsize;
    clock_t start, end;
    double totaltime =0.0;
    double error=0.0;

   /*
    if(rank==0)
    {
        printf("What is np?\n");
        scanf("%d",&np);
        printf("\nnp is %d\n",np);
    }
    MPI_Bcast(&np,1,MPI_INT,0,MPI_COMM_WORLD);
    */

    int *recvcounts=(int*)malloc(sizeof(int)*size);       /* the number of receving data of each rank (gatherv) */
    int *senddisp=(int*)malloc(sizeof(int)*size);         /* the displacement of data size at each rank (scatterv) */
    int *recvdisp=(int*)malloc(sizeof(int)*size);         /* the displacement of data size at each rank (gatherv) */
    int *sendcounts=(int*)malloc(sizeof(int)*size);       /* the number of sending data of each rank (scatterv) */
	double *u = NULL;
    double *u_m= NULL;
    double *X,*Y =NULL;
    double L1error=0.0;
    myrowsize=(int)(np/size)+2;


    senddisp[0]=0;
    recvdisp[0]=0;
    for(i=0;i<size;i++)
    {
        if(i == size-1)
        {
            sendcounts[i]=np*np-(np*myrowsize*(size-1)-2*np*(size-2))+2*np;
            recvcounts[i]=sendcounts[i]-np;
        }
        else if(i != size-1)
        {
            sendcounts[i]=np*myrowsize;
            recvcounts[i]=sendcounts[i]-2*np;
            recvcounts[0]=sendcounts[0]-np;
            senddisp[i+1]=senddisp[i]+sendcounts[i]-2*np;
            recvdisp[i+1]=recvdisp[i]+recvcounts[i];
        }
    }

    partialnum=sendcounts[rank];                         /*the number of data at each rank*/
    realrowsize=partialnum/np;
    if(size==1)
    {
        recvcounts[0]=sendcounts[0];
    }

    double *u_k=(double*)malloc(sizeof(double)*recvcounts[rank]);
    double *urecv=(double*)malloc(sizeof(double)*partialnum);

    if(rank==0)
    {
        u=(double*)malloc(sizeof(double)*np*np);
        X = (double*)malloc(sizeof(double)*np*np);
	    Y = (double*)malloc(sizeof(double)*np*np);

        for(j=0;j<np;j++)
        {
            for(i=0;i<np;i++)
            {
                u[i+j*np]=0.0;
                X[i+j*np]=del_x*i;
                Y[i+j*np]=del_y*j;
            }
        }
        start = clock();
    }

    int k=0;
    do{
        if(rank==0)
        {
            //dirichlet boundary condition//
		    for (j = 0; j < np; j++)
		    {
			    u[j*np] = exact(0, j);
			    u[j*np + np - 1] = exact(np-1, j);
		    }

		    for (i = 0; i < np; i++)
		    {
			    u[i] = exact(i, 0);
			    u[(np - 1)*np + i] = exact(i, np-1);
            }
        }

        L1error=0.0;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatterv(u,sendcounts,senddisp,MPI_DOUBLE,urecv,partialnum,MPI_DOUBLE,0,MPI_COMM_WORLD);

        error=0.0;
        if(rank == 0)
        {
            for(j=1;j<realrowsize-1;j++)
            {
                for(i=1;i<np-1;i++)
                {
                    u_k[j*np + i] = -(f(i,j) - (urecv[j*np + i - 1] + urecv[j*np + i + 1])/pow(del_x, 2) - (urecv[(j - 1)*np + i] + urecv[(j + 1)*np + i])/pow(del_y, 2))/(2 * (1/pow(del_x, 2) + 1/pow(del_y, 2)));

                    error = error + fabs(u_k[j*np+i]-urecv[j*np+i]);
                }
            }

        }
        else if(rank != 0)
        {
            for(j=1;j<realrowsize-1;j++)
            {
                for(i=1;i<np-1;i++)
                {
                    u_k[(j-1)*np + i] = -(f(i,j+rank*(realrowsize-2)) - (urecv[j*np + i - 1] + urecv[j*np + i + 1])/pow(del_x, 2) - (urecv[(j - 1)*np + i] + urecv[(j + 1)*np + i])/pow(del_y, 2))/(2 * (1/pow(del_x, 2) + 1/pow(del_y, 2)));

                    error = error + fabs(u_k[(j-1)*np+i]-urecv[j*np+i]);
                }
            }
        }

        //boundary condition
        for(j=1;j<realrowsize-1;j++)
        {
        u_k[(j-1)*np+(np-1)] = exact(np-1,j+rank*(realrowsize-2));
        u_k[(j-1)*np] = exact(0,j+rank*(realrowsize-2));
        }

        for(i=0;i<np;i++)
        {
            if(rank==size-1)
            {
                u_k[(realrowsize-2)*np+i]=exact(i,np-1);
            }
            else if(rank==0)
            {
                u_k[i]=exact(i,0);
            }
        }

    MPI_Gatherv(u_k,recvcounts[rank],MPI_DOUBLE,u,recvcounts,recvdisp,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Reduce(&error,&L1error,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank ==0)
    {
        printf("L1error: %f\n",L1error);
        k=k+1;

        if(k==200)
        {
        FILE*fp;
        char name[200];
        sprintf(name,"test.dat");
        fp = fopen(name,"w");
        fprintf(fp, "variables= X Y u \n");
        fprintf(fp,"zone i=%d j=%d\n",np,np);
        for(j=0;j<np;j++)
        {
            for(i=0;i<np;i++)
            {
                fprintf(fp,"%f   %f   %f\n",X[j*np+i],Y[j*np+i],u[j*np+i]);
            }
        }
        fclose(fp);
        }
    }

    MPI_Bcast(&L1error,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    }while(L1error>0.0001);




    if(rank==0)
    {
        end = clock();

        totaltime = ((double)(end-start))/CLOCKS_PER_SEC;
        printf("---------------------------\n\n");
        printf("total time : %f \n",totaltime);
        
        FILE*fp;
        char name[200];
        sprintf(name,"jacobi_parallel%d.dat",np);
        fp = fopen(name,"w");
        fprintf(fp, "variables= X Y u \n");
        fprintf(fp,"zone i=%d j=%d\n",np,np);
        for(j=0;j<np;j++)
        {
            for(i=0;i<np;i++)
            {
                fprintf(fp,"%f   %f   %f\n",X[j*np+i],Y[j*np+i],u[j*np+i]);
            }
        }
        fclose(fp);
        

        free(u);
        free(X);
        free(Y);
    }

    free(sendcounts);
    free(recvcounts);
    free(senddisp);
    free(recvdisp);
    free(u_k);
    free(urecv);

    MPI_Finalize();
    return 0;
}
