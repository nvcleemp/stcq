/* PLUGIN file to use with plantri.c 

   To use this, compile plantri.c using 
       cc -o plantri_fkt -O4 '-DPLUGIN="fkt.c"' plantri.c -lm

   This plug-in uses the FKT algorithm to calculate the number
   of perfect matchings in the dual of each generated planar graph.
*/

#define FILTER count_perfect_matchings
#define SUMMARY perfect_matchings_summary

#include <math.h>

static int DFS_seen[MAXF];
static int DFS_stack[NUMEDGES];
static EDGE *DFS_stack_edges[NUMEDGES];

static int unorientedEdgeCount[MAXN];

struct list_el {
    int key;
    int value;
    struct list_el * next;
    struct list_el * prev;
};

typedef struct list_el item;

item *perfect_matchings_counts = NULL;

static int make_dual(void);

item* increment(item* head, int key){
    //first check whether the list is empty
    if(head==NULL){
        item *new = (item *)malloc(sizeof(item));
        new->key = key;
        new->value = 1;
        new->next = NULL;
        new->prev = NULL;
        return new;
    }
    
    //find the position where the new value should be added
    item *currentItem = head;
    while(currentItem->key < key && currentItem->next!=NULL){
        currentItem = currentItem->next;
    }
    
    if(currentItem->key == key){
        currentItem->value++;
        return head;
    } else if(currentItem->key < key){
        item *new = (item *)malloc(sizeof(item));
        new->key = key;
        new->value = 1;
        new->next = NULL;
        new->prev = currentItem;
        currentItem->next = new;
        return head;
    } else if(currentItem == head){
        item *new = (item *)malloc(sizeof(item));
        new->key = key;
        new->value = 1;
        new->next = head;
        new->prev = NULL;
        head->prev = new;
        return new;
    } else {
        item *new = (item *)malloc(sizeof(item));
        new->key = key;
        new->value = 1;
        new->next = currentItem;
        new->prev = currentItem->prev;
        currentItem->prev->next = new;
        currentItem->prev = new;
        return head;
    }
}

void DFS_spanning_tree(int nf){
    int i;
    EDGE *e, *elast;
    
    for(i=0; i<nf; i++){
        DFS_seen[i]=FALSE;
    }
    
    //push first vertex on the stack
    int stackSize = 1;
    DFS_stack[0] = 0;
    DFS_stack_edges[0] = NULL;
    
    while(stackSize>0){
        stackSize--;
        int currentFace = DFS_stack[stackSize];
        if(!DFS_seen[currentFace]){
            //mark vertex as seen
            DFS_seen[currentFace] = TRUE;
            
            if(DFS_stack_edges[stackSize]!=NULL){
                //mark edges of spanning tree
                MARKHI(DFS_stack_edges[stackSize]);
                MARKLO(DFS_stack_edges[stackSize]->invers);
            }
            
            //push next vertices on stack
            e = elast = facestart[currentFace];
            do {
                if(!DFS_seen[e->invers->rightface]){
                    DFS_stack[stackSize]=e->invers->rightface;
                    DFS_stack_edges[stackSize]=e;
                    stackSize++;
                }
                e = e->invers->prev;
            } while (e != elast);
        }
    }
}

long determinant(int n, int mat[n][n]){
    int i,j,i_count,j_count, count=0;
    int array[n-1][n-1];
    long det=0;

    if(n==1) return mat[0][0];
    if(n==2) return (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]);
    if(n==3) return (mat[0][0]*((mat[1][1]*mat[2][2])-(mat[1][2]*mat[2][1])))
                  - (mat[0][1]*((mat[1][0]*mat[2][2])-(mat[1][2]*mat[2][0])))
                  + (mat[0][2]*((mat[1][0]*mat[2][1])-(mat[1][1]*mat[2][0])));
    if(n==4) {
        int m1 = mat[1][1]*(mat[2][2]*mat[3][3]-mat[2][3]*mat[3][2])
               - mat[1][2]*(mat[2][1]*mat[3][3]-mat[2][3]*mat[3][1])
               + mat[1][3]*(mat[2][1]*mat[3][2]-mat[2][2]*mat[3][1]);
        int m2 = mat[1][0]*(mat[2][2]*mat[3][3]-mat[2][3]*mat[3][2])
               - mat[1][2]*(mat[2][0]*mat[3][3]-mat[2][3]*mat[3][0])
               + mat[1][3]*(mat[2][0]*mat[3][2]-mat[2][2]*mat[3][0]);
        int m3 = mat[1][0]*(mat[2][1]*mat[3][3]-mat[2][3]*mat[3][1])
               - mat[1][1]*(mat[2][0]*mat[3][3]-mat[2][3]*mat[3][0])
               + mat[1][3]*(mat[2][0]*mat[3][1]-mat[2][1]*mat[3][0]);
        int m4 = mat[1][0]*(mat[2][1]*mat[3][2]-mat[2][2]*mat[3][1])
               - mat[1][1]*(mat[2][0]*mat[3][2]-mat[2][2]*mat[3][0])
               + mat[1][2]*(mat[2][0]*mat[3][1]-mat[2][1]*mat[3][0]);
        return mat[0][0] * m1
             - mat[0][1] * m2
             + mat[0][2] * m3
             - mat[0][3] * m4;
    }
 
    for(count=0; count<n; count++) {
        i_count=0;
        for(i=1; i<n; i++) {
            j_count=0;
            for(j=0; j<n; j++) {
                if(j == count) continue;
                array[i_count][j_count] = mat[i][j];
                j_count++;
            }
            i_count++;
        }
        if(mat[0][count]){
            if(count%2){
                det +=   mat[0][count] * determinant(n-1,array);	//Recursive call
            } else {
                det += - mat[0][count] * determinant(n-1,array);	//Recursive call
            }
        }
    }
    return det;
}

static int count_perfect_matchings(int nbtot, int nbop, int doflip) {
    int i, j;
    EDGE *e, *elast, *eunmarked;
    
    if(missing_vertex >= 0){
        //disk triangulations not supported
        return FALSE;
    }
    for(i=0; i<nv; i++) unorientedEdgeCount[i] = 0;
    
    //we need the dual graph for FKT
    //nf stores the number of faces
    //(i.e. the number of vertices in the dual graph) (Duh!)
    int nf = make_dual();
    
    RESETMARKS
    
    //construct spanning tree of graph using DFS
    DFS_spanning_tree(nf);
    
    //count the number of unmarked (directed) edges
    int unmarked = 0;
    for(i=0; i<nv; i++){
        e = elast = firstedge[i];
        do {
            if(!ISMARKED(e)){
                unmarked++;
                unorientedEdgeCount[i]++;
            }
            e = e->next;
        } while (e != elast);
    }
    
    while(unmarked>0){
        //find a leaf
        i = 0;
        while(i<nv && unorientedEdgeCount[i]!=1) i++;
        if(i==nv){
            //error
            fprintf(stderr, "Illegal state of program.\n");
            exit(1);
        }
        
        //find the unoriented edge and count the number of clock-wise oriented edges
        int clockwise = 0;
        e = elast = firstedge[i];
        do {
            if(ISMARKEDHI(e)){
                clockwise++;
            } else if(!ISMARKEDLO(e)){
                eunmarked = e;
            }
            e = e->next;
        } while(e != elast);
        
        if(clockwise%2){
            MARKLO(eunmarked);
            MARKHI(eunmarked->invers);
        } else {
            MARKLO(eunmarked->invers);
            MARKHI(eunmarked);
        }
        
        unmarked-=2; //2 directed edges were marked
        
        //adjust the degrees of the dual tree
        unorientedEdgeCount[eunmarked->start]--;
        unorientedEdgeCount[eunmarked->end]--;
    }
    
    //calculate the pfaffian
    int am[nf][nf];
    
    for(i=0; i<nf; i++){
        for(j=0; j<nf; j++){
            am[i][j] = 0;
        }
    }
    for(i=0; i<nf; i++){
        e = elast = facestart[i];
        do {
            if(ISMARKEDHI(e)){
                am[i][e->invers->rightface] = 1;
            } else {
                am[i][e->invers->rightface] = -1;
            }
            e = e->invers->prev;
        } while (e != elast);
    }
    
    /*
    for(i=0; i<nf; i++){
        for(j=0; j<nf; j++){
            fprintf(stderr, "%2d ", am[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
     */
    
    long det = determinant(nf, am);
    
    int pmCount = (int)(sqrt(det)+0.5);
    
    perfect_matchings_counts = increment(perfect_matchings_counts, pmCount);
    
    return TRUE; //don't output graphs
}

void perfect_matchings_summary() {  
    item *currentItem = perfect_matchings_counts;
    fprintf(stderr, "Size   Count\n");
    fprintf(stderr, "------------\n");
    while(currentItem!=NULL){
        fprintf(stderr, "%4d : %5d\n", currentItem->key, currentItem->value);
        currentItem = currentItem->next;
    }
}