#include<stdio.h>
#include<stdlib.h>

struct Vertex
{
    int vertex_no;
    struct Vertex *Next;
    struct Adj_list *Neighbours;
    int neighbours_count;
};

struct Adj_list
{
    int vertex_no;
    struct Adj_list *Next;
    float jc,ks,ct,sat_ct_score;
    int flag,sat_ct;
};

struct Vertex* Create_graph(char []);
void add_in_Adjlist(struct Vertex*, struct Adj_list*);
void add_edge(int,int,struct Vertex**);
struct Vertex* complement_graph(struct Vertex*);
void delete_edge(struct Vertex**,int,int);
float*** Hitting_time_prob(int,struct Vertex*,int);
float Jaccard(int,int,struct Vertex*);
float** katz(struct Vertex*,int);
void jc_sort(struct Vertex*,int,int);
void ks_sort(struct Vertex*,int,int);
void ct_sort(struct Vertex*,int,int);
void sat_ct_sort(struct Vertex*,int,int);



struct Vertex* Create_graph(char file_path[])
{
    FILE *fp;
	struct Vertex *header;
	header = NULL;
	int i,j,w;
	fp = fopen(file_path,"r");
	while(fscanf(fp,"%d %d %d",&i,&j,&w) == 3)
    {
        add_edge(i,j,&header);
        add_edge(j,i,&header);
	}
	fclose(fp);

	return header;
}

// add given edge in Adjacency list of given vertex in increasing order of vertex no. appearing in Adjacency list of given vertex
void add_in_Adjlist(struct Vertex *vptr, struct Adj_list *eptr)
{
    struct Adj_list *temp_e;   //temporary edge pointer to traverse the Adjacency list of given vertex
    temp_e = vptr->Neighbours;
    if (eptr->vertex_no < temp_e->vertex_no)
    {
        eptr->Next = temp_e;
        vptr->Neighbours = eptr;
    }
    else
    {
        while(temp_e->Next != NULL && (temp_e->Next)->vertex_no < eptr->vertex_no) {temp_e = temp_e->Next;}
        if (temp_e->Next != NULL) {temp_e->Next = eptr; temp_e = temp_e->Next; eptr->Next = temp_e;}
        else {temp_e->Next = eptr;}
    }
    vptr->neighbours_count++;
    return;
}

// add edge in graph adjacency list in sorted way
void add_edge(int i, int j, struct Vertex **header)
{
    struct Vertex *vptr,*temp_v,*prev_temp_v;
    struct Adj_list *eptr,*temp_e;
    eptr = (struct Adj_list*)malloc(sizeof(struct Adj_list));
    eptr->vertex_no = j;
    eptr->Next = NULL;
//    printf("%d %d\n",i,j);

    if ((*header) == NULL)  // first edge case
    {
        vptr = (struct Vertex*)malloc(sizeof(struct Vertex));
        vptr->vertex_no = i;
        vptr->Neighbours = eptr;
        vptr->neighbours_count = 1;
        vptr->Next = NULL;
        *header = vptr;
        return;
    }

    temp_v = *header;
    while(temp_v->Next!=NULL) {temp_v = temp_v->Next;}

    if (temp_v->vertex_no < i) // when new edge has first vertex no. more than highest vertex no. in graph
    {
        vptr = (struct Vertex*)malloc(sizeof(struct Vertex));
        vptr->vertex_no = i;
        vptr->Neighbours = eptr;
        vptr->neighbours_count = 1;
        vptr->Next = NULL;
        temp_v->Next = vptr;
        return;
    }

    else
    {
        temp_v = *header;
        while(temp_v->vertex_no < i) {prev_temp_v = temp_v; temp_v = temp_v->Next;}
        if (temp_v->vertex_no == i) {add_in_Adjlist(temp_v, eptr); return;} // when new edge has first vertex no. less than highest and vertex no exists in graph
        else  // when vertex no doesn't exist in graph
        {
            vptr = (struct Vertex*)malloc(sizeof(struct Vertex));
            vptr->vertex_no = i;
            vptr->Neighbours = eptr;
            vptr->neighbours_count = 1;
            vptr->Next = temp_v;
            prev_temp_v->Next = vptr;
            return;
        }
    }
    return;
}

// find graph with non-existent edges of given graph
struct Vertex* complement_graph(struct Vertex *header)
{
    int vertex_count;
    struct Vertex *vptr,*c_header=NULL;
    struct Adj_list *eptr,*ehead=NULL;
    FILE *fptr;

    // find vertex count
    for(vptr=header;vptr->Next!=NULL;vptr=vptr->Next) {continue;}
    vertex_count = vptr->vertex_no;


    // create full graph from all possible edges except self-edges
    for(int i=vertex_count;i>0;i--)
    {
        vptr = (struct Vertex*)malloc(sizeof(struct Vertex));
        vptr->vertex_no = i;
        vptr->Next = c_header;
        c_header=vptr;
        ehead=NULL;

        for(int j=vertex_count;j>0;j--)
        {
            if(j==i) {continue;} // self edge excluded
            eptr = (struct Adj_list*)malloc(sizeof(struct Adj_list));
            eptr->vertex_no = j;
            eptr->Next = ehead;
            ehead=eptr;
        }
        vptr->Neighbours = ehead;
        vptr->neighbours_count = vertex_count-1;
    }


    // delete edges of given graph from full graph
    for(vptr=header; vptr!=NULL; vptr=vptr->Next)
    {
        for(eptr=vptr->Neighbours; eptr!=NULL; eptr=eptr->Next)
        {
            delete_edge(&c_header, vptr->vertex_no, eptr->vertex_no);
        }
    }
    return c_header;
}



void delete_edge(struct Vertex **header, int i, int j)
{
    struct Vertex *vptr,*prev_vptr;
    struct Adj_list *eptr,*prev_eptr;

    vptr = *header;
    while(vptr->vertex_no!=i) {prev_vptr=vptr; vptr=vptr->Next;}
    eptr=vptr->Neighbours;
    while(eptr->vertex_no!=j) {prev_eptr=eptr; eptr=eptr->Next;}

    if(eptr==vptr->Neighbours) {vptr->Neighbours=eptr->Next; free(eptr);}
    else {prev_eptr->Next=eptr->Next; free(eptr);}
    (vptr->neighbours_count)--;
    if((vptr->neighbours_count)==0)
    {
        if(vptr==*header) {*header=vptr->Next; free(vptr);}
        else {prev_vptr->Next=vptr->Next; free(vptr);}
    }
}



float*** Hitting_time_prob(int k, struct Vertex *header, int vertex_count)
{
    float ***P; // prob. matrix which includes prob. transition matrix and hitting time prob. upto k step
    struct Vertex *vptr;
    struct Adj_list *eptr;
    float sum;


    P = (float ***)malloc(vertex_count*sizeof(float **));
    for(int i=0;i<vertex_count;i++)
    {
        P[i] = (float **)malloc(vertex_count*sizeof(float *));
        for(int j=0;j<vertex_count;j++)
        {
            P[i][j] = (float *)calloc(k,sizeof(float));
        }
    }


    // prob. transition matrix
    for(vptr=header;vptr!=NULL;vptr=vptr->Next)
    {
        for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
        {
            P[(vptr->vertex_no)-1][(eptr->vertex_no)-1][0] = 1/(float)(vptr->neighbours_count);
        }
    }

    // computation of prob of different hitting time steps
    for(int m=1;m<k;m++)
    {
        for(int i=0;i<vertex_count;i++)
        {
            for(int j=0;j<vertex_count;j++)
            {
                for(int l=0;l<vertex_count;l++)
                {
                    if(l!=j) {P[i][j][m] += P[i][l][0]*P[l][j][m-1];}
                }
            }
        }
    }


    return P;
}


float Jaccard(int i, int j, struct Vertex *header)
{
    struct Vertex *vptr1, *vptr2;
    struct Adj_list *eptr1, *eptr2;
    float union_count=0, intersection_count=0,jc;

    vptr1 = header;
    vptr2 = header;

    while(vptr1->vertex_no!=i) {vptr1=vptr1->Next;}
    while(vptr2->vertex_no!=j) {vptr2=vptr2->Next;}

    eptr1=vptr1->Neighbours;
    eptr2=vptr2->Neighbours;


    while(eptr1!=NULL && eptr2!=NULL)
    {
        if(eptr1->vertex_no==eptr2->vertex_no)
        {
            intersection_count++;
            eptr1=eptr1->Next;
            eptr2=eptr2->Next;
        }
        else if(eptr1->vertex_no>eptr2->vertex_no) {eptr2=eptr2->Next;}

        else {eptr1=eptr1->Next;}
    }

    union_count = vptr1->neighbours_count + vptr2->neighbours_count - intersection_count;
    jc = intersection_count / union_count;

    return jc;
}


float** katz(struct Vertex *header, int vertex_count)
{
    struct Vertex *vptr;
    struct Adj_list *eptr;
    int*** adj_matrix; // adjacency matrix and matrices indicating no. of 2 to 6 length paths
    float** katz_matrix;
    float b;



	adj_matrix = (int ***)malloc(vertex_count*sizeof(int **));
	for(int i=0;i<vertex_count;i++)
        {
            adj_matrix[i] = (int **)malloc(vertex_count*sizeof(int *));
            for(int j=0;j<vertex_count;j++)
            {
                adj_matrix[i][j] = (int *)calloc(6,sizeof(int));
            }
        }


    katz_matrix = (float **)malloc(vertex_count*sizeof(float *));
    for(int i=0;i<vertex_count;i++)
        {
            katz_matrix[i] = (float *)calloc(vertex_count,sizeof(float));
        }

    // created adjacency matrix
    for(vptr=header;vptr!=NULL;vptr=vptr->Next)
    {
        for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
        {
            adj_matrix[(vptr->vertex_no)-1][(eptr->vertex_no)-1][0] = 1;
        }
    }


    // matrices indicating no of paths for l=2 to 6
    for(int k=1;k<6;k++)
    {
        for(int i=0;i<vertex_count;i++)
        {
            for(int j=0;j<vertex_count;j++)
            {
                for(int l=0;l<vertex_count;l++)
                {
                    adj_matrix[i][j][k] += adj_matrix[i][l][0]*adj_matrix[l][j][k-1];
                }
            }
        }
    }

    // katz matrix computation
    for(int i=0;i<vertex_count;i++)
    {
        for(int j=0;j<vertex_count;j++)
        {
            b=0.01;
            for(int l=1;l<6;l++){katz_matrix[i][j] += b*adj_matrix[i][j][l]; b=b*0.1;}
        }
    }

    return katz_matrix;
}


void jc_sort(struct Vertex *header, int vertex_count, int K)
{
    int flag[vertex_count][vertex_count]; // flag matrix(1/0) to indicate whether i-j edge is already taken as max score edge
    struct Vertex *vptr;
    struct Adj_list *eptr;
    float max=0;
    int s,d;
    FILE *fp;


    // initialize flag matrix as 0
    fp = fopen("Jaccard.txt","w");
    for(int i=0;i<vertex_count;i++)
    {
        for(int j=0;j<vertex_count;j++)
        {
            flag[i][j] = 0;
        }
    }


    for(int k=0;k<K;k++)
    {
        max = 0;
        for(vptr=header;vptr!=NULL;vptr=vptr->Next)
        {
            for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
            {
                if(flag[vptr->vertex_no-1][eptr->vertex_no-1]==0)
                {
                    if(eptr->jc > max) {max = eptr->jc; s = vptr->vertex_no; d = eptr->vertex_no;}
                    else if(eptr->jc == max)
                    {
                        if(vptr->vertex_no < s) {s = vptr->vertex_no; d = eptr->vertex_no;}
                        else if(vptr->vertex_no == s && eptr->vertex_no < d) {d = eptr->vertex_no;}
                    }
                }
            }
        }
        flag[s-1][d-1] = 1;
        flag[d-1][s-1] = 1;
        fprintf(fp,"%d %d %f\n", s,d,max);
    }
    fclose(fp);
}


void ks_sort(struct Vertex *header, int vertex_count, int K)
{
    int flag[vertex_count][vertex_count];
    struct Vertex *vptr;
    struct Adj_list *eptr;
    float max=0;
    int s,d;
    FILE *fp;

    fp = fopen("Katz.txt","w");
    for(int i=0;i<vertex_count;i++)
    {
        for(int j=0;j<vertex_count;j++)
        {
            flag[i][j] = 0;
        }
    }

    for(int k=0;k<K;k++)
    {
        max = 0;
        for(vptr=header;vptr!=NULL;vptr=vptr->Next)
        {
            for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
            {
                if(flag[vptr->vertex_no-1][eptr->vertex_no-1]==0)
                {
                    if(eptr->ks > max) {max = eptr->ks; s = vptr->vertex_no; d = eptr->vertex_no;}
                    else if(eptr->ks == max)
                    {
                        if(vptr->vertex_no < s) {s = vptr->vertex_no; d = eptr->vertex_no;}
                        else if(vptr->vertex_no == s && eptr->vertex_no < d) {d = eptr->vertex_no;}
                    }
                }
            }
        }
        flag[s-1][d-1] = 1;
        flag[d-1][s-1] = 1;
        fprintf(fp,"%d %d %f\n", s,d,max);
    }
    fclose(fp);
}


void ct_sort(struct Vertex *header, int vertex_count, int K)
{
    int flag[vertex_count][vertex_count];
    struct Vertex *vptr;
    struct Adj_list *eptr;
    float max=-9999;
    int s,d;
    FILE *fp;

    fp=fopen("HittingTime.txt","w");
    for(int i=0;i<vertex_count;i++)
    {
        for(int j=0;j<vertex_count;j++)
        {
            flag[i][j] = 0;
        }
    }

    for(int k=0;k<K;k++)
    {
        max = -9999;
        for(vptr=header;vptr!=NULL;vptr=vptr->Next)
        {
            for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
            {
                if(flag[vptr->vertex_no-1][eptr->vertex_no-1]==0)
                {
                    if(eptr->ct > max) {max = eptr->ct; s = vptr->vertex_no; d = eptr->vertex_no;}
                    else if(eptr->ct == max)
                    {
                        if(vptr->vertex_no < s) {s = vptr->vertex_no; d = eptr->vertex_no;}
                        else if(vptr->vertex_no == s && eptr->vertex_no < d) {d = eptr->vertex_no;}
                    }
                }
            }
        }
        flag[s-1][d-1] = 1;
        flag[d-1][s-1] = 1;
        fprintf(fp,"%d %d %f\n", s,d,max);
    }
    fclose(fp);
}

void sat_ct_sort(struct Vertex *header, int vertex_count, int K)
{
    int flag[vertex_count][vertex_count];
    struct Vertex *vptr;
    struct Adj_list *eptr;
    float max=-9999;
    int s,d,max_ct;
    FILE *fp;

    fp=fopen("HittingTimeAccurate.txt","w");
    for(int i=0;i<vertex_count;i++)
    {
        for(int j=0;j<vertex_count;j++)
        {
            flag[i][j] = 0;
        }
    }

    for(int k=0;k<K;k++)
    {
        max = -9999;
        for(vptr=header;vptr!=NULL;vptr=vptr->Next)
        {
            for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
            {
                if(flag[vptr->vertex_no-1][eptr->vertex_no-1]==0)
                {
                    if(eptr->sat_ct_score > max) {max = eptr->sat_ct_score; max_ct = eptr->sat_ct; s = vptr->vertex_no; d = eptr->vertex_no;}
                    else if(eptr->sat_ct_score == max)
                    {
                        if(vptr->vertex_no < s) {s = vptr->vertex_no; d = eptr->vertex_no; max_ct = eptr->sat_ct;}
                        else if(vptr->vertex_no == s && eptr->vertex_no < d) {d = eptr->vertex_no; max_ct = eptr->sat_ct;}
                    }
                }
            }
        }
        flag[s-1][d-1] = 1;
        flag[d-1][s-1] = 1;
        fprintf(fp,"%d %d %f %d\n", s,d,max,max_ct);
    }
    fclose(fp);
}

int main()
{
    struct Vertex *header,*vptr,*complement_header;
    struct Adj_list *eptr;
    int vertex_count,K;


    header = Create_graph("../contact-high-school-proj-graph.txt");


    for(vptr=header;vptr!=NULL;vptr=vptr->Next) {if (vptr->Next == NULL) vertex_count = vptr->vertex_no;}


    float** katz_matrix;
    katz_matrix = katz(header,vertex_count);


    float*** P;
    float*** HT;
    float*** CT;

    P = Hitting_time_prob(20,header,vertex_count);

    HT = (float ***)malloc(vertex_count*sizeof(float **));
    CT = (float ***)malloc(vertex_count*sizeof(float **));
    for(int i=0;i<vertex_count;i++)
    {
        HT[i] = (float **)malloc(vertex_count*sizeof(float *));
        CT[i] = (float **)malloc(vertex_count*sizeof(float *));
        for(int j=0;j<vertex_count;j++)
        {
            HT[i][j] = (float *)calloc(22,sizeof(float));
            CT[i][j] = (float *)calloc(22,sizeof(float));
        }
    }

    // hitting time computation
    for (int k=2;k<=21;k++)
    {
        for(int i=0;i<vertex_count;i++)
        {
            for(int j=0;j<vertex_count;j++)
            {
                for(int m=0;m<k;m++)
                {
                    HT[i][j][k] += (m+1)*P[i][j][m];
                }
            }
        }
    }


    for(int k=2;k<=21;k++)
    {
        for(int i=0;i<vertex_count;i++)
        {
            for(int j=0;j<vertex_count;j++)
            {
                CT[i][j][k] = HT[i][j][k] + HT[j][i][k];
            }
        }
    }


    complement_header = complement_graph(header);


    for(vptr=complement_header;vptr!=NULL;vptr=vptr->Next)
    {
        for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
        {
            eptr->jc = Jaccard(vptr->vertex_no, eptr->vertex_no, header);
            eptr->ks = katz_matrix[(vptr->vertex_no)-1][(eptr->vertex_no)-1];
            eptr->ct = CT[(vptr->vertex_no)-1][(eptr->vertex_no)-1][6];
            eptr->ct = -(eptr->ct);
            eptr->sat_ct = 20;
            eptr->flag = 0;
        }
    }

    for(int k=2;k<=20;k++)
    {
        for(vptr=complement_header;vptr!=NULL;vptr=vptr->Next)
        {
            for(eptr=vptr->Neighbours;eptr!=NULL;eptr=eptr->Next)
            {
                if (eptr->flag==0)
                {
                    if(CT[(vptr->vertex_no)-1][(eptr->vertex_no)-1][k] - CT[(vptr->vertex_no)-1][(eptr->vertex_no)-1][k+1] <= 0.01 && CT[(vptr->vertex_no)-1][(eptr->vertex_no)-1][k] - CT[(vptr->vertex_no)-1][(eptr->vertex_no)-1][k+1] >= -0.01)
                    {
                        eptr->flag = 1;
                        eptr->sat_ct = k;
                        eptr->sat_ct_score = -CT[(vptr->vertex_no)-1][(eptr->vertex_no)-1][k];
                    }
                }
            }
        }
    }


    printf("enter no. of non-existent edges to be predicted: ");
    scanf("%d", &K);

//    printf("\n\n\nJaccard coefficient:\n");
    jc_sort(complement_header, vertex_count, K);

//    printf("\n\n\nKatz's Score:\n");
    ks_sort(complement_header, vertex_count, K);

//    printf("\n\n\nCommute time:\n");
    ct_sort(complement_header, vertex_count, K);

    sat_ct_sort(complement_header, vertex_count, K);

	return 0;
}
