/* This program reads quadrangulations from standard in and
 * generates spherical tilings by congruent quadrangles, where
 * the quadrangles are all of type 2, i.e. the length of three sides are 
 * equal and the length of the fourth side is different.   
 * 
 * 
 * Compile with:
 *     
 *     cc -o stcq -O4 stcq_sa.c liblpsolve55.a
 */

#include "lp_lib.h"

#ifndef MAXN
#define MAXN 64            /* the maximum number of vertices */
#endif
#define MAXE (4*MAXN-8)    /* the maximum number of oriented edges */
#define MAXF (MAXN-2)      /* the maximum number of faces */
#define MAXVAL (MAXN-2)/2  /* the maximum degree of a vertex */
#define MAXCODELENGTH (MAXN+MAXE+3)

#undef FALSE
#undef TRUE
#define FALSE 0
#define TRUE  1

typedef struct e /* The data type used for edges */ {
    int start; /* vertex where the edge starts */
    int end; /* vertex where the edge ends */
    int rightface; /* face on the right side of the edge
                          note: only valid if make_dual() called */
    struct e *prev; /* previous edge in clockwise direction */
    struct e *next; /* next edge in clockwise direction */
    struct e *inverse; /* the edge that is inverse to this one */
    int mark,index;    /* two ints for temporary use;
                          Only access mark via the MARK macros. */

    int left_facesize; /* size of the face in prev-direction of the edge.
        		  Only used for -p option. */
    char angle; /* angle between this edge and next edge;
                   0: alpha, 1: beta, 2: gamma, 3: delta */
} EDGE;

EDGE *firstedge[MAXN]; /* pointer to arbitrary edge out of vertex i. */
int degree[MAXN];

EDGE *facestart[MAXF]; /* pointer to arbitrary edge of face i. */
int faceSize[MAXF]; /* pointer to arbitrary edge of face i. */

EDGE edges[MAXE];

static int markvalue = 30000;
#define RESETMARKS {int mki; if ((markvalue += 2) > 30000) \
       { markvalue = 2; for (mki=0;mki<MAXE;++mki) edges[mki].mark=0;}}
#define MARK(e) (e)->mark = markvalue
#define MARKLO(e) (e)->mark = markvalue
#define MARKHI(e) (e)->mark = markvalue+1
#define UNMARK(e) (e)->mark = markvalue-1
#define ISMARKED(e) ((e)->mark >= markvalue)
#define ISMARKEDLO(e) ((e)->mark == markvalue)
#define ISMARKEDHI(e) ((e)->mark > markvalue)


unsigned long long int numberOfQuadrangulations = 0;
unsigned long long int rejectedByHammingDistance = 0;

int printDuplicateEquations = TRUE;

int unusedSwitch = FALSE; // if set to TRUE: unused graphs will be written to stdout
int usedSwitch = FALSE; // if set to TRUE: used graphs will be written to stdout

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

int nv; //the number of vertices of the current quadrangulation
int nf; //the number of faces of the current quadrangulation

/*
 * The following variable stores the direction in which the edges of the face
 * should be iterated over when assigning the angles alpha, beta, gamma and
 * delta. The assignment always starts with the edge stored in matchingEdges.
 * 
 * Possible directions are 0 and 1.
 */
int angleAssigmentDirection[MAXF];

char angleAroundVertex[MAXN][MAXN];

unsigned long long int unusedGraphCount = 0;

//////////////////////////////////////////////////////////////////////////////

void printPlanarGraph(){
    int i;
    for(i=0; i<nv; i++){
        fprintf(stderr, "%d: ", i);
        EDGE *e, *elast;
    
        e = elast = firstedge[i];
        do {
            fprintf(stderr, "%d ", e->end);
            e = e->next;
        } while (e!=elast);
        fprintf(stderr, "\n");
    }
}

void printAngleAssignment(){
    int i;
    for(i=0; i<nv; i++){
        fprintf(stderr, "%d: ", i);
        EDGE *e, *elast;
    
        e = elast = firstedge[i];
        do {
            fprintf(stderr, "%d ", e->end);
            fprintf(stderr, "(%c) ", 'a' + e->angle);
            e = e->next;
        } while (e!=elast);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printFaceMatching(){
    printPlanarGraph();
    
    int i;
    for (i = 0; i < nv - 2; i++) {
        EDGE *e = matchingEdges[i];
        if(e->start < e->end){
            fprintf(stderr, "%d - %d\n", e->start, e->end);
        }
    }
    fprintf(stderr, "\n");
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

item* increment(item* head, int key) {
    //first check whether the list is empty
    if (head == NULL) {
        item *new = (item *) malloc(sizeof (item));
        new->key = key;
        new->value = 1;
        new->next = NULL;
        new->prev = NULL;
        return new;
    }

    //find the position where the new value should be added
    item *currentItem = head;
    while (currentItem->key < key && currentItem->next != NULL) {
        currentItem = currentItem->next;
    }

    if (currentItem->key == key) {
        currentItem->value++;
        return head;
    } else if (currentItem->key < key) {
        item *new = (item *) malloc(sizeof (item));
        new->key = key;
        new->value = 1;
        new->next = NULL;
        new->prev = currentItem;
        currentItem->next = new;
        return head;
    } else if (currentItem == head) {
        item *new = (item *) malloc(sizeof (item));
        new->key = key;
        new->value = 1;
        new->next = head;
        new->prev = NULL;
        head->prev = new;
        return new;
    } else {
        item *new = (item *) malloc(sizeof (item));
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

void handleSolution() {

}

unsigned long long int solvable = 0;

void printSystem() {
    int i;
    for (i = 0; i < nv; i++) {
        if (printDuplicateEquations || !isDuplicateEquation[i]) {
            fprintf(stderr, "(%d,%d,%d,%d)\n", alphaCount[i], betaCount[i], gammaCount[i], deltaCount[i]);
        }
    }
    fprintf(stderr, "\n");
}

void simplifySystem() {
    //nothing to do at the moment
}

void solveSystem() {
    lprec *lp;
    int *colno = NULL, i, j;
    REAL *row = NULL;

    lp = make_lp(0, 4);
    /* There are 4 angles.
     */
    if (lp == NULL) {
        exit(1);
    }
    resize_lp(lp, nv + 1 - duplicateEquationCount, get_Ncolumns(lp));
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

    colno = (int *) malloc(4 * sizeof (*colno));
    row = (REAL *) malloc(4 * sizeof (*row));
    if ((colno == NULL) || (row == NULL)) {
        exit(1);
    }

    colno[0] = 1;
    row[0] = 1;

    if (!set_obj_fnex(lp, 1, row, colno)) {
        exit(1);
    }

    set_add_rowmode(lp, TRUE); //start adding rows

    for (i = 0; i < nv; i++) {
        if (!isDuplicateEquation[i]) {
            j = 0;

            if (alphaCount[i] != 0) {
                colno[j] = 1;
                row[j++] = alphaCount[i];
            }

            if (betaCount[i] != 0) {
                colno[j] = 2;
                row[j++] = betaCount[i];
            }

            if (gammaCount[i] != 0) {
                colno[j] = 3;
                row[j++] = gammaCount[i];
            }

            if (deltaCount[i] != 0) {
                colno[j] = 4;
                row[j++] = deltaCount[i];
            }

            if (!add_constraintex(lp, j, row, colno, EQ, 2)) {
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

    if (!add_constraintex(lp, 4, row, colno, EQ, 2 + 4.0 / (nv - 2))) {
        exit(1);
    }

    set_add_rowmode(lp, FALSE); //stop adding rows

    set_maxim(lp);

#ifdef _DEBUG
    write_LP(lp, stderr);
#endif

    set_verbose(lp, SEVERE);

    int result = solve(lp);

    if (result == OPTIMAL) {
        solvable++;
        handleSolution();
    } else if (printUnsolvableSystems || writeLpsolveUnsolvedSystems) {
        printSystem();
    }

#ifdef _DEBUG
    fprintf(stderr, "Objective value: %f\n", get_objective(lp));

    get_variables(lp, row);
    for (j = 0; j < 4; j++) {
        fprintf(stderr, "%s: %f\n", get_col_name(lp, j + 1), row[j]);
    }

    fprintf(stderr, "\n");
#endif

    free(row);
    free(colno);

    delete_lp(lp);
}

unsigned long long int assignmentCount = 0;

void createSystem() {
    //clear systems
    int i;
    for (i = 0; i < nv; i++) {
        alphaCount[i] = betaCount[i] = gammaCount[i] = deltaCount[i] = 0;
    }

    //iterate over all faces
    for (i = 0; i < nv - 2; i++) {
        EDGE *e1 = matchingEdges[i];
        EDGE *e2 = e1->inverse->prev;
        EDGE *e3 = e2->inverse->prev;
        EDGE *e4 = e3->inverse->prev;

        //assert: e1 = e4->inverse->prev;
        if (angleAssigmentDirection[i]) {
            alphaCount[e1->end] += 1;
            betaCount[e2->end] += 1;
            gammaCount[e3->end] += 1;
            deltaCount[e4->end] += 1;
            e2->angle = 0;
            e3->angle = 1;
            e4->angle = 2;
            e1->angle = 3;
        } else {
            alphaCount[e4->end] += 1;
            betaCount[e3->end] += 1;
            gammaCount[e2->end] += 1;
            deltaCount[e1->end] += 1;
            e1->angle = 0;
            e4->angle = 1;
            e3->angle = 2;
            e2->angle = 3;
        }
    }
}

int firstCheckOfSystem() {
    int i, j;
    /* If the system contains two equations that are at hamming distance 1,
     * we return FALSE. At the same time we remove duplicate equations (i.e.
     * equations that are at Hamming distance 0).
     */

    //reset array
    for (i = 0; i < nv; i++) {
        isDuplicateEquation[i] = FALSE;
    }
    duplicateEquationCount = 0;

    for (i = 0; i < nv - 1; i++) {
        if (isDuplicateEquation[i]) continue;
        /*
         * this equation is already a duplicate itself, so any equation
         * that would be at Hamming distance 1 of this equation is also
         * at Hamming Distance of that earlier equation.
         */
        for (j = i + 1; j < nv; j++) {
            if (isDuplicateEquation[j]) continue;
            /*
             * this equation is already a duplicate, so any equation
             * that would be at Hamming distance 1 of this equation is also
             * at Hamming Distance of that earlier equation.
             */
            int hammingDistance = 0;
            if (alphaCount[i] != alphaCount[j]) hammingDistance++;
            if (betaCount[i] != betaCount[j]) hammingDistance++;
            if (gammaCount[i] != gammaCount[j]) hammingDistance++;
            if (deltaCount[i] != deltaCount[j]) hammingDistance++;
            if (hammingDistance == 1) {
                return FALSE;
            } else if (hammingDistance == 0) {
                // if we get here, then equation j is a duplicate of i
                isDuplicateEquation[j] = TRUE;
                duplicateEquationCount++;
            }
        }
    }
    return TRUE;
}

void handleAngleAssignment() {
    assignmentCount++;
    createSystem();
    if (firstCheckOfSystem()) {
        simplifySystem();
        solveSystem();
    } else {
        rejectedByHammingDistance++;
        if (printUnsolvableSystems || writeHammingDistanceUnsolvedSystems) {
            simplifySystem();
            printSystem();
        }
    }
#ifdef _DEBUG
    printSystem();
#endif
}

void assignAnglesForCurrentPerfectMatchingRecursion(int currentFace) {
    if (currentFace == nv - 2) {
        handleAngleAssignment();
    } else {
        angleAssigmentDirection[currentFace] = 0;
        assignAnglesForCurrentPerfectMatchingRecursion(currentFace + 1);
        angleAssigmentDirection[currentFace] = 1;
        assignAnglesForCurrentPerfectMatchingRecursion(currentFace + 1);
    }
}

void assignAnglesForCurrentPerfectMatching() {
    //we fix the direction of one face to prevent mirror images to both be generated.
    angleAssigmentDirection[0] = 0;
    assignAnglesForCurrentPerfectMatchingRecursion(1);
}

int matchingCount = 0;

void handlePerfectMatching() {
    matchingCount++;
    assignAnglesForCurrentPerfectMatching();
}

void matchNextFace(int lastFace, int matchingSize) {
    if (matchingSize == (nv - 2) / 2) {
        //Found a perfect matching
        handlePerfectMatching();
        return;
    }

    int nextFace = lastFace;
    while (nextFace < nv - 2 && matched[nextFace]) nextFace++;

    if (nextFace == nv - 2) {
        fprintf(stderr, "Something went terribly wrong.");
        exit(1);
    }
    matched[nextFace] = TRUE;

    EDGE *e, *elast;

    e = elast = facestart[nextFace];
    do {
        int neighbour = e->inverse->rightface;
        if (!matched[neighbour]) {
            match[nextFace] = neighbour;
            match[neighbour] = nextFace;
            matched[neighbour] = TRUE;

            //store edge corresponding to the match for each face.
            matchingEdges[nextFace] = e;
            matchingEdges[neighbour] = e->inverse;

            matchNextFace(nextFace, matchingSize + 1);

            matched[neighbour] = FALSE;
        }
        e = e->inverse->prev;
    } while (e != elast);

    matched[nextFace] = FALSE;
}

int generate_perfect_matchings_in_dual() {
    int i;

    matchingCount = 0;

    unsigned long long int oldSolutionCount = solvable;

    if (nf != nv - 2) {
        fprintf(stderr, "Something went horribly wrong. Maybe some wrong parameter?\nnf: %d, nv: %d\n", nf, nv);
        exit(1);
    }

    for (i = 0; i < nv - 2; i++) {
        matched[i] = FALSE;
    }
    matched[0] = TRUE;

    EDGE *e, *elast;

    e = elast = facestart[0];
    do {
        int neighbour = e->inverse->rightface;
        match[0] = neighbour;
        match[neighbour] = 0;
        matched[neighbour] = TRUE;

        //store edge corresponding to the match for each face.
        matchingEdges[0] = e;
        matchingEdges[neighbour] = e->inverse;

        matchNextFace(0, 1);

        matched[neighbour] = FALSE;

        e = e->inverse->prev;
    } while (e != elast);


    perfect_matchings_counts = increment(perfect_matchings_counts, matchingCount);

    if (oldSolutionCount == solvable) {
        unusedGraphCount++;
        if (unusedSwitch) {
            return 1;
        }
    } else if (unusedSwitch) {
        return 0;
    } else if (usedSwitch) {
        return 1;
    }

    return 0;
}

void perfect_matchings_summary() {
    unsigned long long int totalPerfectMatchingsCount = 0;
    item *currentItem = perfect_matchings_counts;
    fprintf(stderr, "Size   Count\n");
    fprintf(stderr, "------------\n");
    while (currentItem != NULL) {
        fprintf(stderr, "%4d : %5d\n", currentItem->key, currentItem->value);
        totalPerfectMatchingsCount += (currentItem->value) * (currentItem->key);
        currentItem = currentItem->next;
    }
    fprintf(stderr, "\nQuadrangulations: %llu\n", numberOfQuadrangulations);
    fprintf(stderr, "\nMatchings: %llu\n", totalPerfectMatchingsCount);
    fprintf(stderr, "\nAssignments: %llu\n", assignmentCount);
    fprintf(stderr, "\nSolvable: %llu\n", solvable);
    fprintf(stderr, "\nNon-solvable: %llu\n", assignmentCount - solvable);
    fprintf(stderr, "\n%llu quadrangulations don't correspond to a tiling.\n", unusedGraphCount);
    if (printStatistics) {
        fprintf(stderr, "\nRejected by Hamming distance: %llu\n", rejectedByHammingDistance);
        fprintf(stderr, "Rejected by lpsolve: %llu\n\n", assignmentCount - solvable - rejectedByHammingDistance);
    }
}

EDGE *findEdge(int from, int to){
    EDGE *e, *elast;
    
    e = elast = firstedge[from];
    do {
        if(e->end==to){
            return e;
        }
        e = e->next;
    } while (e!=elast);
    fprintf(stderr, "error while looking for edge from %d to %d.\n", from, to);
    exit(0);
}

 
/* Store in the rightface field of each edge the number of the face on
   the right hand side of that edge.  Faces are numbered 0,1,....  Also
   store in facestart[i] an example of an edge in the clockwise orientation
   of the face boundary, and the size of the face in facesize[i], for each i.
   Returns the number of faces. */
void makeDual(){
    register int i,sz;
    register EDGE *e,*ex,*ef,*efx;
 
    RESETMARKS;
 
    nf = 0;
    for (i = 0; i < nv; ++i){

        e = ex = firstedge[i];
        do
        {
            if (!ISMARKEDLO(e))
            {
                facestart[nf] = ef = efx = e;
                sz = 0;
                do
                {
                    ef->rightface = nf;
                    MARKLO(ef);
                    ef = ef->inverse->prev;
                    ++sz;
                } while (ef != efx);
                faceSize[nf] = sz;
                ++nf;
            }
            e = e->next;
        } while (e != ex);
    }
}


void decodePlanarCode(unsigned short* code) {
    /* complexity of method to determine inverse isn't that good, but will have to satisfy for now
     */
    int i, j, codePosition;
    int edgeCounter = 0;
    EDGE *inverse;

    nv = code[0];
    codePosition = 1;

    for (i = 0; i < nv; i++) {
        degree[i] = 0;
        firstedge[i] = edges + edgeCounter;
        edges[edgeCounter].start = i;
        edges[edgeCounter].end = code[codePosition] - 1;
        edges[edgeCounter].next = edges + edgeCounter + 1;
        if(code[codePosition] - 1 < i){
            inverse = findEdge(code[codePosition]-1, i);
            edges[edgeCounter].inverse = inverse;
            inverse->inverse = edges + edgeCounter;
        } else {
            edges[edgeCounter].inverse = NULL;
        }
        edgeCounter++;
        codePosition++;
        for (j = 1; code[codePosition]; j++, codePosition++) {
            if (j == MAXVAL) {
                fprintf(stderr, "MAXVAL too small: %d\n", MAXVAL);
                exit(0);
            }
            edges[edgeCounter].start = i;
            edges[edgeCounter].end = code[codePosition] - 1;
            edges[edgeCounter].prev = edges + edgeCounter - 1;
            edges[edgeCounter].next = edges + edgeCounter + 1;
            if(code[codePosition] - 1 < i){
                inverse = findEdge(code[codePosition]-1, i);
                edges[edgeCounter].inverse = inverse;
                inverse->inverse = edges + edgeCounter;
            } else {
                edges[edgeCounter].inverse = NULL;
            }
            edgeCounter++;
        }
        firstedge[i]->prev = edges + edgeCounter - 1;
        edges[edgeCounter-1].next = firstedge[i];
        degree[i] = j;
        
        codePosition++; /* read the closing 0 */
    }
    
    makeDual();
}

/**
 * 
 * @param code
 * @param laenge
 * @param file
 * @return returns 1 if a code was read and 0 otherwise. Exits in case of error.
 */
int readPlanarCode(unsigned short code[], int *length, FILE *file) {
    static int first = 1;
    unsigned char c;
    char testheader[20];
    int bufferSize, zeroCounter;


    if (first) {
        first = 0;

        if (fread(&testheader, sizeof (unsigned char), 15, file) != 15) {
            fprintf(stderr, "can't read header ((1)file too small)-- exiting\n");
            exit(1);
        }
        testheader[15] = 0;
        if (strcmp(testheader, ">>planar_code<<") == 0) {
            
        } else {
            fprintf(stderr, "No planarcode header detected -- exiting!\n");
            exit(1);
        }
    }

    /* possibly removing interior headers -- only done for planarcode */
    if (fread(&c, sizeof (unsigned char), 1, file) == 0){
        //nothing left in file
        return (0);
    }
    
    if (c == '>'){
        // could be a header, or maybe just a 62 (which is also possible for unsigned char
        code[0] = c;
        bufferSize = 1;
        zeroCounter = 0;
        code[1] = (unsigned short) getc(file);
        if (code[1] == 0) zeroCounter++;
        code[2] = (unsigned short) getc(file);
        if (code[2] == 0) zeroCounter++;
        bufferSize = 3;
        // 3 characters were read and stored in buffer
        if ((code[1] == '>') && (code[2] == 'p')) /*we are sure that we're dealing with a header*/ {
            while ((c = getc(file)) != '<');
            /* read 2 more characters: */
            c = getc(file);
            if (c != '<') {
                fprintf(stderr, "Problems with header -- single '<'\n");
                exit(1);
            }
            if (!fread(&c, sizeof (unsigned char), 1, file)){
                //nothing left in file
                return (0);
            }
            bufferSize = 1;
            zeroCounter = 0;
        }
    } else {
        //no header present
        bufferSize = 1;
        zeroCounter = 0;
    }

    if (c != 0) /* unsigned chars would be sufficient */ {
        code[0] = c;
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        while (zeroCounter < code[0]) {
            code[bufferSize] = (unsigned short) getc(file);
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    } else {
        fread(code, sizeof (unsigned short), 1, file);
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        bufferSize = 1;
        zeroCounter = 0;
        while (zeroCounter < code[0]) {
            fread(code + bufferSize, sizeof (unsigned short), 1, file);
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    }

    *length = bufferSize;
    return (1);


}

int main(int argc, char *argv[]){
    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        numberOfQuadrangulations++;
        generate_perfect_matchings_in_dual();
    }
    perfect_matchings_summary();

}