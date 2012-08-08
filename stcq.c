/* PLUGIN file to use with plantri.c 

   To use this, compile plantri.c using 
       cc -o plantri_stcq -O4 '-DPLUGIN="stcq.c"' plantri.c

   This plug-in generates spherical tilings by congruent quadrangles, where
   the quadrangles are all of type 2, i.e. the length of three sides are 
   equal and the length of the fourth side is different.   
*/

#include "lp_lib.h"

#define FILTER generate_perfect_matchings_in_dual
#define PLUGIN_INIT init_plugin()
#define SUMMARY perfect_matchings_summary
#define PLUGIN_SWITCHES else if(arg[j]=='n'){\
                            unusedSwitch = TRUE;\
                        } else if(arg[j]=='z'){\
                            int switchvalue = getswitchvalue(arg,&j);\
                            if(switchvalue==1){\
                                printUnsolvableSystems = TRUE;\
                            } else if(switchvalue==2){\
                                printStatistics = TRUE;\
                            }\
                        } else if(arg[j]=='w'){\
                            int switchvalue = getswitchvalue(arg,&j);\
                            if(switchvalue==1){\
                                writeLpsolveUnsolvedSystems = TRUE;\
                            } else if(switchvalue==2){\
                                writeHammingDistanceUnsolvedSystems = TRUE;\
                            }\
                        } else if(arg[j]=='D'){\
                            printDuplicateEquations = FALSE;\
                        }

unsigned long long int rejectedByHammingDistance = 0;

int printDuplicateEquations = TRUE;

int unusedSwitch = FALSE; // if set to TRUE: unused graphs will be written to stdout

int printUnsolvableSystems = FALSE; //1
int printStatistics = FALSE; //2

int writeLpsolveUnsolvedSystems = FALSE; //1
int writeHammingDistanceUnsolvedSystems = FALSE; //2

int matched[MAXF];
int match[MAXF];
EDGE *matchingEdges[MAXF];

int alphaCount[MAXN];
int betaCount[MAXN];
int gammaCount[MAXN];
int deltaCount[MAXN];

int isDuplicateEquation[MAXN];
int duplicateEquationCount = 0;

/*
 * The following variable stores the direction in which the edges of the face
 * should be iterated over when assigning the angles alpha, beta, gamma and
 * delta. The assignment always starts with the edge stored in matchingEdges.
 * 
 * Possible directions are 0 and 1.
 */
int angleAssigmentDirection[MAXF];

unsigned long long int unusedGraphCount = 0;

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

void handleSolution(){
    
}

unsigned long long int solvable = 0;

void printSystem(){
    int i;
    for(i=0; i<nv; i++){
        if(printDuplicateEquations || !isDuplicateEquation[i]){
            fprintf(stderr, "(%d,%d,%d,%d)\n", alphaCount[i], betaCount[i], gammaCount[i], deltaCount[i]);
        }
    }
    fprintf(stderr, "\n");
}

void simplifySystem(){
    //nothing to do at the moment
}

void solveSystem(){
    lprec *lp;
    int *colno = NULL, i, j;
    REAL *row = NULL;
    
    lp = make_lp(0,4);
    /* There are 4 angles.
     */
    if(lp == NULL){
        exit(1);
    }
    resize_lp(lp, nv+1 - duplicateEquationCount, get_Ncolumns(lp));
    /* There nv + 1 equations: one for each vertex plus the extra equation
     * 
     *     alpha + beta + gamma + delta = 2 +4/F.
     * 
     * Of course we only add each distinct equation once, so we subtract the
     * number of duplicate equations.
     */
    
    //name the columns
    set_col_name(lp, 1, "alpha");
    set_col_name(lp, 2, "beta");
    set_col_name(lp, 3, "gamma");
    set_col_name(lp, 4, "delta");
    
    REAL epsilon = 0.0000001;
    REAL lowerBoundAngle = 0 + epsilon;
    REAL upperBoundAngle = 2 - epsilon;
    set_bounds(lp, 1, lowerBoundAngle, upperBoundAngle);
    set_bounds(lp, 2, lowerBoundAngle, upperBoundAngle);
    set_bounds(lp, 3, lowerBoundAngle, upperBoundAngle);
    set_bounds(lp, 4, lowerBoundAngle, upperBoundAngle);
    
    colno = (int *) malloc(4 * sizeof(*colno));
    row = (REAL *) malloc(4 * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        exit(1);
    }
    
    colno[0] = 1;
    row[0] = 1;
    
    if(!set_obj_fnex(lp, 1, row, colno)){
        exit(1);
    }
    
    set_add_rowmode(lp, TRUE); //start adding rows
    
    for(i=0; i<nv; i++){
        if(!isDuplicateEquation[i]) {
            j = 0;

            if(alphaCount[i]!=0){
                colno[j] = 1;
                row[j++] = alphaCount[i];
            }

            if(betaCount[i]!=0){
                colno[j] = 2;
                row[j++] = betaCount[i];
            }

            if(gammaCount[i]!=0){
                colno[j] = 3;
                row[j++] = gammaCount[i];
            }

            if(deltaCount[i]!=0){
                colno[j] = 4;
                row[j++] = deltaCount[i];
            }

            if(!add_constraintex(lp, j, row, colno, EQ, 2)){
                exit(1);
            }
        }
    }
    
    colno[0] = 1;
    row[0] = 1;
    colno[1] = 2;
    row[1] = 1;
    colno[2] = 3;
    row[2] = 1;
    colno[3] = 4;
    row[3] = 1;
    
    if(!add_constraintex(lp, 4, row, colno, EQ, 2 + 4.0/(nv-2))){
        exit(1);
    }
    
    set_add_rowmode(lp, FALSE); //stop adding rows
    
    set_maxim(lp);
    
#ifdef _DEBUG
    write_LP(lp, stderr);
#endif
    
    set_verbose(lp, IMPORTANT);
    
    int result = solve(lp);
    
    if(result == OPTIMAL){
        solvable++;
        handleSolution();
    } else if(printUnsolvableSystems || writeLpsolveUnsolvedSystems){
        printSystem();
    }
    
#ifdef _DEBUG
    fprintf(stderr, "Objective value: %f\n", get_objective(lp));

    get_variables(lp, row);
    for(j = 0; j < 4; j++){
      fprintf(stderr, "%s: %f\n", get_col_name(lp, j + 1), row[j]);
    }
    
    fprintf(stderr, "\n");
#endif
    
    free(row);
    free(colno);
    
    delete_lp(lp);
}

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

int firstCheckOfSystem(){
    int i, j;
    /* If the system contains two equations that are at hamming distance 1,
     * we return FALSE. At the same time we remove duplicate equations (i.e.
     * equations that are at Hamming distance 0).
     */

    //reset array
    for(i=0; i<nv; i++){
        isDuplicateEquation[i] = FALSE;
    }
    duplicateEquationCount = 0;
    
    for(i=0; i<nv-1; i++){
        if(isDuplicateEquation[i]) continue;
        /*
         * this equation is already a duplicate itself, so any equation
         * that would be at Hamming distance 1 of this equation is also
         * at Hamming Distance of that earlier equation.
         */
        for(j=i+1; j<nv; j++){
            if(isDuplicateEquation[j]) continue;
            /*
            * this equation is already a duplicate, so any equation
            * that would be at Hamming distance 1 of this equation is also
            * at Hamming Distance of that earlier equation.
            */
            int hammingDistance = 0;
            if(alphaCount[i]!=alphaCount[j]) hammingDistance++;
            if(betaCount[i]!=betaCount[j]) hammingDistance++;
            if(gammaCount[i]!=gammaCount[j]) hammingDistance++;
            if(deltaCount[i]!=deltaCount[j]) hammingDistance++;
            if(hammingDistance==1){
                return FALSE;
            } else if(hammingDistance==0){
                // if we get here, then equation j is a duplicate of i
                isDuplicateEquation[j] = TRUE;
                duplicateEquationCount++;
            }
        }
    }
    return TRUE;
}

void handleAngleAssignment(){
    assignmentCount++;
    createSystem();
    if(firstCheckOfSystem()){
        simplifySystem();
        solveSystem();
    } else {
        rejectedByHammingDistance++;
        if(printUnsolvableSystems || writeHammingDistanceUnsolvedSystems){
            simplifySystem();
            printSystem();
        }
    }
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
    int i;
    
    matchingCount = 0;
    
    unsigned long long int oldSolutionCount = solvable;
    
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
    
    if(oldSolutionCount == solvable){
        unusedGraphCount++;
        if(unusedSwitch){
            return 1;
        }
    } else if(unusedSwitch){
        return 0;
    }
    
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
    fprintf(stderr, "\nSolvable: %llu\n", solvable);
    fprintf(stderr, "\nNon-solvable: %llu\n", assignmentCount - solvable);
    fprintf(stderr, "\n%llu quadrangulations don't correspond to a tiling.\n", unusedGraphCount);
    if(printStatistics){
        fprintf(stderr, "\nRejected by Hamming distance: %llu\n", rejectedByHammingDistance);
        fprintf(stderr, "Rejected by lpsolve: %llu\n\n", assignmentCount - solvable - rejectedByHammingDistance);
    }
}