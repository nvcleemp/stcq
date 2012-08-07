/* PLUGIN file to use with plantri.c 

   To use this, compile plantri.c using 
       cc -o plantri_stcq -O4 '-DPLUGIN="stcq.c"' plantri.c

   This plug-in generates spherical tilings by congruent quadrangles, where
   the quadrangles are all of type 2, i.e. the length of three sides are 
   equal and the length of the fourth side is different.   
*/

#define FILTER generate_perfect_matchings_in_dual
#define PLUGIN_INIT init_plugin()
#define SUMMARY perfect_matchings_summary

int matched[MAXF];
int match[MAXF];
EDGE *matchingEdges[MAXF];

int alphaCount[MAXN];
int betaCount[MAXN];
int gammaCount[MAXN];
int deltaCount[MAXN];

/*
 * The following variable stores the direction in which the edges of the face
 * should be iterated over when assigning the angles alpha, beta, gamma and
 * delta. The assignment always starts with the edge stored in matchingEdges.
 * 
 * Possible directions are 0 and 1.
 */
int angleAssigmentDirection[MAXF];

static int make_dual(void);

void init_plugin(){
    qswitch = TRUE;
}


//////////////////////////////////////////////////////////////////////////////

struct list_el {
    int key;
    int value;
    struct list_el * next;
    struct list_el * prev;
};

typedef struct list_el item;

item *perfect_matchings_counts = NULL;

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

//////////////////////////////////////////////////////////////////////////////

unsigned long long int assignmentCount = 0;

void createSystem(){
    //clear systems
    int i;
    for(i=0; i<nv; i++){
        alphaCount[i] = betaCount[i] = gammaCount[i] = deltaCount[i] = 0;
    }
    
    //iterate over all faces
    for(i=0; i<nv-2; i++){
        EDGE *e1 = matchingEdges[i];
        EDGE *e2 = e1->invers->prev;
        EDGE *e3 = e2->invers->prev;
        EDGE *e4 = e3->invers->prev;

        //assert: e1 = e4->invers->prev;
        if(angleAssigmentDirection[i]){
            alphaCount[e1->end] += 1;
            betaCount[e2->end] += 1;
            gammaCount[e3->end] += 1;
            deltaCount[e4->end] += 1;
        } else {
            alphaCount[e4->end] += 1;
            betaCount[e3->end] += 1;
            gammaCount[e2->end] += 1;
            deltaCount[e1->end] += 1;
        }
    }
}

void printSystem(){
    int i;
    for(i=0; i<nv; i++){
        fprintf(stderr, "(%d,%d,%d,%d)\n", alphaCount[i], betaCount[i], gammaCount[i], deltaCount[i]);
    }
    fprintf(stderr, "\n");
}

void handleAngleAssignment(){
    assignmentCount++;
    createSystem();
#ifdef _DEBUG
    printSystem();
#endif
}

void assignAnglesForCurrentPerfectMatchingRecursion(int currentFace){
    if(currentFace==nv-2){
        handleAngleAssignment();
    } else {
            angleAssigmentDirection[currentFace] = 0;
            assignAnglesForCurrentPerfectMatchingRecursion(currentFace + 1);
            angleAssigmentDirection[currentFace] = 1;
            assignAnglesForCurrentPerfectMatchingRecursion(currentFace + 1);
    }
}    

void assignAnglesForCurrentPerfectMatching(){
    //we fix the direction of one face to prevent mirror images to both be generated.
    angleAssigmentDirection[0] = 0;
    assignAnglesForCurrentPerfectMatchingRecursion(1);
}

int matchingCount = 0;

void handlePerfectMatching(){
    matchingCount++;
    assignAnglesForCurrentPerfectMatching();
}

void matchNextFace(int lastFace, int matchingSize){
    if(matchingSize==(nv-2)/2){
        //Found a perfect matching
        handlePerfectMatching();
        return;
    }
    
    int nextFace = lastFace;
    while(nextFace<nv-2 && matched[nextFace]) nextFace++;
    
    if(nextFace == nv-2){
        fprintf(stderr, "Something went terribly wrong.");
        exit(1);
    }
    matched[nextFace] = TRUE;
    
    EDGE *e, *elast;
    
    e = elast = facestart[nextFace];
    do {
        int neighbour = e->invers->rightface;
        if(!matched[neighbour]){
            match[nextFace] = neighbour;
            match[neighbour] = nextFace;
            matched[neighbour] = TRUE;
            
            //store edge corresponding to the match for each face.
            matchingEdges[nextFace] = e;
            matchingEdges[neighbour] = e->invers;

            matchNextFace(nextFace, matchingSize+1);

            matched[neighbour] = FALSE;
        }
        e = e->invers->prev;
    } while(e != elast);
    
    matched[nextFace] = FALSE;
}

static int generate_perfect_matchings_in_dual(int nbtot, int nbop, int doflip) {
    int i, j;
    
    matchingCount = 0;
    
    
    //we need the dual graph
    //nf stores the number of faces
    //(i.e. the number of vertices in the dual graph) (Duh!)
    int nf = make_dual();
    if(nf != nv-2){
        fprintf(stderr, "Something went horribly wrong. Maybe some wrong parameter?\n");
        exit(1);
    }
    
    for(i=0; i<nv-2; i++){
        matched[i] = FALSE;
    }
    matched[0] = TRUE;
    
    EDGE *e, *elast;
    
    e = elast = facestart[0];
    do {
        int neighbour = e->invers->rightface;
        match[0] = neighbour;
        match[neighbour] = 0;
        matched[neighbour] = TRUE;
            
        //store edge corresponding to the match for each face.
        matchingEdges[0] = e;
        matchingEdges[neighbour] = e->invers;
        
        matchNextFace(0, 1);
        
        matched[neighbour] = FALSE;
        
        e = e->invers->prev;
    } while(e != elast);
    
    
    perfect_matchings_counts = increment(perfect_matchings_counts, matchingCount);
    
    return 0;
}

void perfect_matchings_summary() {  
    item *currentItem = perfect_matchings_counts;
    fprintf(stderr, "Size   Count\n");
    fprintf(stderr, "------------\n");
    while(currentItem!=NULL){
        fprintf(stderr, "%4d : %5d\n", currentItem->key, currentItem->value);
        currentItem = currentItem->next;
    }
    fprintf(stderr, "\nAssignments: %llu\n", assignmentCount);
}