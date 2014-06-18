/* This program reads quadrangulations from standard in and
 * generates spherical tilings by congruent quadrangles, where
 * the quadrangles are all of type 2, i.e. the length of three sides are 
 * equal and the length of the fourth side is different.   
 * 
 * 
 * Compile with:
 *     
 *     cc -o stcq -O4 stcq_sa.c liblpsolve55.a
 * 
 * On some systems it might be necessary to compile with:
 *     
 *     cc -o stcq -O4 stcq_sa.c liblpsolve55.a -lm -ldl -lssp
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "lp_lib.h"

#ifndef MAXN
#define MAXN 64            /* the maximum number of vertices */
#endif
#define MAXE (4*MAXN-8)    /* the maximum number of oriented edges */
#define MAXF (MAXN-2)      /* the maximum number of faces */
#define MAXVAL (MAXN-2)/2  /* the maximum degree of a vertex */
#define MAXCODELENGTH (MAXN+MAXE+3)

#define INFI (MAXN + 1)

#undef FALSE
#undef TRUE
#define FALSE 0
#define TRUE  1

typedef int boolean;

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
    int allowedInFaceMatching;
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
unsigned long long int rejectedByCoefficientDiff = 0;

boolean printDuplicateEquations = FALSE;

unsigned long long int filterOnly = 0;

boolean onlyConvex = TRUE;

char outputFormat = 'n'; //defaults to no output

char generatedType = 't'; //defaults to spherical tilings

boolean unusedQuadrangulations = FALSE; // if set to TRUE: unused quadrangulations will be written to stdout
boolean usedQuadrangulations = FALSE; // if set to TRUE: used quadrangulations will be written to stdout

boolean outputSolution = TRUE; //by default we output the solution

boolean isEarlyFilteringEnabled = TRUE;
boolean generateAllMatchings = FALSE;
boolean boundAngleAssignments = TRUE;

boolean printUnsolvableSystems = FALSE; //1
boolean printStatistics = FALSE; //2

boolean writeLpsolveUnsolvedSystems = FALSE; //1
boolean writeHammingDistanceUnsolvedSystems = FALSE; //2

boolean oneBased = FALSE;

boolean includeGroup = FALSE;

FILE *latexSummaryFile = NULL;
boolean latexPerSolution = FALSE;
char *latexBaseName = NULL;
char latexFileNameBuffer[100];

boolean matched[MAXF];
int match[MAXF];
EDGE *matchingEdges[MAXF];

int alphaCount[MAXN];
int betaCount[MAXN];
int gammaCount[MAXN];
int deltaCount[MAXN];

boolean isDuplicateEquation[MAXN];
int duplicateEquationCount = 0;

int degreeThreeVertexTypes[MAXN];

int nv; //the number of vertices of the current quadrangulation
int nf; //the number of faces of the current quadrangulation
int ne; //the number of edges of the current quadrangulation

int orderedFaces[MAXF];
int faceRank[MAXF]; //inverse of orderedFaces
boolean checkVerticesAfterFace[MAXF]; //if true for index i, then the vertex restriction
                                  //should be tested after i faces have been assigned
boolean vertexCompletedAfterFace[MAXN];

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

unsigned long long int solvable = 0;

unsigned long long int solvableAndCanonical = 0;

int quadrangulationAutomorphisms[4*MAXE][MAXN]; //there are at most 4e automorphisms
int quadrangulationAutomorphismsCount;

int aaAutomorphisms[4*MAXE][MAXN]; //there are at most 4e automorphisms
int aaAutomorphismsCount;
boolean aaAutomorphismGroupContainsOrientationReversingSymmetry;

boolean generateSTCQ4 = FALSE;

boolean relabelInputQuadrangulation = FALSE;

boolean mirrorImagesAreDistinct = FALSE;

int splitting_res = 0;
int splitting_mod = 1;
int splitting_level = 10;
int splitting_count;

//////////////////////////////////////////////////////////////////////////////

void calculateAutomorphismGroupAngleAssignments();

//////////////////////////////////////////////////////////////////////////////

boolean degreeThreeTypesCompatibility[10][10] =
{
    { TRUE, FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,  TRUE, FALSE,  TRUE},
    {FALSE,  TRUE,  TRUE,  TRUE,  TRUE, FALSE, FALSE, FALSE,  TRUE, FALSE},
    {FALSE,  TRUE,  TRUE,  TRUE, FALSE,  TRUE, FALSE, FALSE, FALSE,  TRUE},
    { TRUE,  TRUE,  TRUE,  TRUE, FALSE, FALSE,  TRUE, FALSE, FALSE, FALSE},
    { TRUE,  TRUE, FALSE, FALSE,  TRUE,  TRUE, FALSE,  TRUE, FALSE, FALSE},
    {FALSE, FALSE,  TRUE, FALSE,  TRUE,  TRUE, FALSE, FALSE,  TRUE,  TRUE},
    {FALSE, FALSE, FALSE,  TRUE, FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE},
    { TRUE, FALSE, FALSE, FALSE,  TRUE, FALSE,  TRUE,  TRUE,  TRUE, FALSE},
    {FALSE,  TRUE, FALSE, FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE, FALSE},
    { TRUE, FALSE,  TRUE, FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,  TRUE}
};

int degreeThreeTypesCombinationVertexLowerBound[10][10] =
{
    {   0, INFI, INFI,    0,    0, INFI, INFI,    8, INFI,   10},
    {INFI,    0,   10,    0,    0, INFI, INFI, INFI,    0, INFI},
    {INFI,   10,    0,   10, INFI,    8, INFI, INFI, INFI,    0},
    {   0,    0,   10,    0, INFI, INFI,    0, INFI, INFI, INFI},
    {   0,    0, INFI, INFI,    0,   10, INFI,    0, INFI, INFI},
    {INFI, INFI,   10, INFI,   10,    0, INFI, INFI,    0,    0},
    {INFI, INFI, INFI,    0, INFI, INFI,    0,   10,    0,    0},
    {   8, INFI, INFI, INFI,    0, INFI,   10,    0,   10, INFI},
    {INFI,    0, INFI, INFI, INFI,    0,    0,   10,    0, INFI},
    {  10, INFI,    0, INFI, INFI,    0,    0, INFI, INFI,    0}
};

int degreeThreeTypesCombinationVertexUpperBound[10][10] =
{
    {INFI,    0,    0, INFI, INFI,    0,    0,    8,    0, INFI},
    {   0, INFI, INFI, INFI, INFI,    0,    0,    0, INFI,    0},
    {   0, INFI, INFI, INFI,    0,    8,    0,    0,    0, INFI},
    {INFI, INFI, INFI, INFI,    0,    0, INFI,    0,    0,    0},
    {INFI, INFI,    0,    0, INFI, INFI,    0, INFI,    0,    0},
    {   0,    0,    8,    0, INFI, INFI,    0,    0, INFI, INFI},
    {   0,    0,    0, INFI,    0,    0, INFI, INFI, INFI, INFI},
    {   8,    0,    0,    0, INFI,    0, INFI, INFI, INFI,    0},
    {   0, INFI,    0,    0,    0, INFI, INFI, INFI, INFI,    0},
    {INFI,    0, INFI,    0,    0, INFI, INFI,    0,    0, INFI}
};

int getDegreeThreeVertexType(int v){
    int a, b, c, d;
    a = alphaCount[v];
    b = betaCount[v];
    c = gammaCount[v];
    d = deltaCount[v];
    if(a+b+c+d!=3){
        fprintf(stderr, "Something went wrong. :-(\n");
        exit(0);
    }
    if(b==3){
        return 0;
    } else if(b==2 && c==1){
        return 1;
    } else if(a==1 && d==1 && b==1){
        return 2;
    } else if(a==2 && c==1){
        return 3;
    } else if(a==2 && b==1){
        return 4;
    } else if(c==3){
        return 5;
    } else if(c==2 && b==1){
        return 6;
    } else if(a==1 && d==1 && c==1){
        return 7;
    } else if(d==2 && b==1){
        return 8;
    } else{
        return 9;
    }
}

//////////////////////////////////////////////////////////////////////////////

void printGroupElement(FILE *f, int *groupElement, int offset){
    int i, next;
    boolean printed[MAXN];
    
    for(i=0; i<MAXN; i++) printed[i] = FALSE;

    for(i=0; i < nv; i++){
        if(!printed[i]){
            fprintf(f, "(%d", i + offset);
            printed[i] = 1;
            next = groupElement[i];
            while(!printed[next]){
                fprintf(f, " %d", next + offset);
                printed[next] = TRUE;
                next = groupElement[next];
            }
            fprintf(f, ") ");
        }
    }
}

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

void printAngleAssignmentLatex(){
    if(latexPerSolution){
        //open file
        int result = sprintf(latexFileNameBuffer, latexBaseName, solvableAndCanonical);
        if(result>0){
           latexSummaryFile = fopen(latexFileNameBuffer, "w"); 
        } else {
            fprintf(stderr, "Error creating filename for LaTeX output -- exiting!\n");
            exit(EXIT_FAILURE);
        }
    }
    
    if(latexSummaryFile==NULL) return;
    int i;
    
    //start with the group
    if(includeGroup){
        calculateAutomorphismGroupAngleAssignments();
        fprintf(latexSummaryFile, "automorphism count: %d\\\\\n", aaAutomorphismsCount);
        for(i=0; i<aaAutomorphismsCount; i++){
            printGroupElement(latexSummaryFile, aaAutomorphisms[i], oneBased);
            fprintf(latexSummaryFile, "\\\\\n");
        }
    }
    
    for(i=0; i<nv; i++){
        fprintf(latexSummaryFile, "%d: ", i + oneBased);
        EDGE *e, *elast;
    
        e = elast = firstedge[i];
        do {
            fprintf(latexSummaryFile, "%d ", e->end + oneBased);
            if(e->angle==0)
                fprintf(latexSummaryFile, "($\\alpha$) ");
            else if(e->angle==1)
                fprintf(latexSummaryFile, "($\\beta$) ");
            else if(e->angle==2)
                fprintf(latexSummaryFile, "($\\gamma$) ");
            else// (e->angle==3)
                fprintf(latexSummaryFile, "($\\delta$) ");
            e = e->next;
        } while (e!=elast);
        fprintf(latexSummaryFile, "\\\\\n");
    }
    fprintf(latexSummaryFile, "\\\\\n");
    
    if(latexPerSolution){
        fclose(latexSummaryFile);
        latexSummaryFile = NULL;
    }
}

void printSphericalTilingByCongruentQuadrangles(lprec *lp){
    int rows = get_Nrows(lp);
    int cols = get_Ncolumns(lp);
    REAL pv[1 + rows + cols];
    get_primal_solution(lp, pv);
    int i;
    for(i=0; i<cols; i++){
        fprintf(stderr, "%c = %f\n", 'a' + i, pv[1 + rows + i]);
    }
    printAngleAssignment();
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

void printQuadrangulationAutomorphismGroup(){    
    int i;
    fprintf(stderr, "automorphism count: %d\n", quadrangulationAutomorphismsCount);
    for (i = 0; i < quadrangulationAutomorphismsCount; i++) {
        printGroupElement(stderr, quadrangulationAutomorphisms[i], 0);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}


//////////////////////////////////////////////////////////////////////////////

/*
fills the array code with the angle assignment code of the current structure.
length will contain the length of the code. The maximum number of vertices is limited
to 255.
*/
void computeAngleAssignmentCode(unsigned char code[], int *length) {
    int i;
    unsigned char *codeStart;
    EDGE *e, *elast;

    codeStart = code;
    *code = (unsigned char) (nv);
    code++;
    for (i = 0; i < nv; i++) {
        e = elast = firstedge[i];
        do {
            *code = (unsigned char) (e->end + 1);
            code++;
            e = e->next;
        } while (e!=elast);
        *code = 0;
        code++;
        do {
            *code = (unsigned char) (e->angle + 1);
            code++;
            e = e->next;
        } while (e!=elast);
        *code = 0;
        code++;
    }
    *length = code - codeStart;
    return;
}

/*
fills the array code with the angle assignment code of the current structure.
length will contain the length of the code. The maximum number of vertices is limited
to 65535.
*/
void computeAngleAssignmentCodeShort(unsigned short code[], int *length) {
    int i;
    unsigned short *codeStart;
    EDGE *e, *elast;

    codeStart = code;
    *code = (unsigned short) (nv);
    code++;
    for (i = 0; i < nv; i++) {
        e = elast = firstedge[i];
        do {
            *code = (unsigned short) (e->end + 1);
            code++;
            e = e->next;
        } while (e!=elast);
        *code = 0;
        code++;
        do {
            *code = (unsigned short) (e->angle + 1);
            code++;
            e = e->next;
        } while (e!=elast);
        *code = 0;
        code++;
    }
    *length = code - codeStart;
    return;
}

void writeAngleAssignment(){
    static boolean first = TRUE;
    
    if(first){
        fprintf(stdout, ">>angle_assignment<<");
        first = FALSE;
    }
    
    int length;
    unsigned char code[MAXE * 2 + MAXN*2 + 1];
    unsigned short codeShort[MAXE * 2 + MAXN*2 + 1];

    if (nv + 1 <= 255) {
        computeAngleAssignmentCode(code, &length);
        if (fwrite(code, sizeof (unsigned char), length, stdout) != length) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    } else if (nv + 1 <= 65535){
        computeAngleAssignmentCodeShort(codeShort, &length);
        putc(0, stdout);
        if (fwrite(codeShort, sizeof (unsigned short), length, stdout) != length) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    } else {
        fprintf(stderr, "Graph too large for angle assignment format -- exiting!\n");
        exit(-1);
    }
}

/*
fills the array code with the planar code of the current quadrangulation.
length will contain the length of the code. The maximum number of vertices is limited
to 255.
*/
void computePlanarCode(unsigned char code[], int *length) {
    int i;
    unsigned char *codeStart;
    EDGE *e, *elast;

    codeStart = code;
    *code = (unsigned char) (nv);
    code++;
    for (i = 0; i < nv; i++) {
        e = elast = firstedge[i];
        do {
            *code = (unsigned char) (e->end + 1);
            code++;
            e = e->next;
        } while (e!=elast);
        *code = 0;
        code++;
    }
    *length = code - codeStart;
    return;
}

/*
fills the array code with the angle assignment code of the current structure.
length will contain the length of the code. The maximum number of vertices is limited
to 65535.
*/
void computePlanarCodeShort(unsigned short code[], int *length) {
    int i;
    unsigned short *codeStart;
    EDGE *e, *elast;

    codeStart = code;
    *code = (unsigned short) (nv);
    code++;
    for (i = 0; i < nv; i++) {
        e = elast = firstedge[i];
        do {
            *code = (unsigned short) (e->end + 1);
            code++;
            e = e->next;
        } while (e!=elast);
        *code = 0;
        code++;
    }
    *length = code - codeStart;
    return;
}

void writePlanarCode(){
    static boolean first = TRUE;
    
    if(first){
        fprintf(stdout, ">>planar_code<<");
        first = FALSE;
    }
    
    int length;
    unsigned char code[MAXE + MAXN + 1];
    unsigned short codeShort[MAXE + MAXN + 1];
    //each directed edge once, plus one zero per vertex, plus size of graphs

    if (nv + 1 <= 255) {
        computePlanarCode(code, &length);
        if (fwrite(code, sizeof (unsigned char), length, stdout) != length) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    } else if (nv + 1 <= 65535){
        computePlanarCodeShort(codeShort, &length);
        putc(0, stdout);
        if (fwrite(codeShort, sizeof (unsigned short), length, stdout) != length) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    } else {
        fprintf(stderr, "Graph too large for angle assignment format -- exiting!\n");
        exit(-1);
    }
}

void outputQuadrangulation(){
    if(outputFormat == 'c'){
        writePlanarCode();
    } else if (outputFormat == 'h'){
        printPlanarGraph();
    }
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

int cagqCertificate[MAXE+MAXN];
int cagqAlternateCertificate[MAXE+MAXN];
int cagqAlternateLabelling[MAXN];
EDGE *alternateFirstedge[MAXN];
int cagqQueue[MAXN];

void constructAlternateCertificate(EDGE *eStart){
    int i;
    for(i=0; i<MAXN; i++){
        cagqAlternateLabelling[i] = MAXN;
    }
    EDGE *e, *elast;
    int head = 1;
    int tail = 0;
    int vertexCounter = 1;
    int cagqAlternateCertificatePosition = 0;
    cagqQueue[0] = eStart->start;
    alternateFirstedge[eStart->start] = eStart;
    cagqAlternateLabelling[eStart->start] = 0;
    while(head>tail){
        int currentVertex = cagqQueue[tail++];
        e = elast = alternateFirstedge[currentVertex];
        do {
            if(cagqAlternateLabelling[e->end]==MAXN){
                cagqQueue[head++] = e->end;
                cagqAlternateLabelling[e->end] = vertexCounter++;
                alternateFirstedge[e->end] = e->inverse;
            }
            cagqAlternateCertificate[cagqAlternateCertificatePosition++] = cagqAlternateLabelling[e->end];
            e = e->next;
        } while (e!=elast);
        cagqAlternateCertificate[cagqAlternateCertificatePosition++] = MAXN;
    }
}

void constructAlternateCertificateOrientationReversing(EDGE *eStart){
    int i;
    for(i=0; i<MAXN; i++){
        cagqAlternateLabelling[i] = MAXN;
    }
    EDGE *e, *elast;
    int head = 1;
    int tail = 0;
    int vertexCounter = 1;
    int cagqAlternateCertificatePosition = 0;
    cagqQueue[0] = eStart->start;
    alternateFirstedge[eStart->start] = eStart;
    cagqAlternateLabelling[eStart->start] = 0;
    while(head>tail){
        int currentVertex = cagqQueue[tail++];
        e = elast = alternateFirstedge[currentVertex];
        do {
            if(cagqAlternateLabelling[e->end]==MAXN){
                cagqQueue[head++] = e->end;
                cagqAlternateLabelling[e->end] = vertexCounter++;
                alternateFirstedge[e->end] = e->inverse;
            }
            cagqAlternateCertificate[cagqAlternateCertificatePosition++] = cagqAlternateLabelling[e->end];
            e = e->prev;
        } while (e!=elast);
        cagqAlternateCertificate[cagqAlternateCertificatePosition++] = MAXN;
    }
}

void calculateAutomorphismGroupQuadrangulation(){
    quadrangulationAutomorphismsCount = 0;
    
    //construct certificate
    int pos = 0;
    int i, j;
    
    for(i=0; i<nv; i++){
        EDGE *e, *elast;

        e = elast = firstedge[i];
        do {
            cagqCertificate[pos++] = e->end;
            e = e->next;
        } while (e!=elast);
        cagqCertificate[pos++] = MAXN;
    }
    
    //construct alternate certificates
    EDGE *ebase = firstedge[0];
    
    for(i=0; i<nv; i++){
        if(degree[i]==degree[0]){
            EDGE *e, *elast;

            e = elast = firstedge[i];
            do {
                if(e!=ebase){
                    constructAlternateCertificate(e);
                    if(memcmp(cagqCertificate, cagqAlternateCertificate, sizeof(int)*pos) == 0) {
                        //store automorphism
                        memcpy(quadrangulationAutomorphisms[quadrangulationAutomorphismsCount], cagqAlternateLabelling, sizeof(int)*MAXN);
                        quadrangulationAutomorphismsCount++;
                    }
                }
                if(!mirrorImagesAreDistinct){
                    constructAlternateCertificateOrientationReversing(e);
                    if(memcmp(cagqCertificate, cagqAlternateCertificate, sizeof(int)*pos) == 0) {
                        //store automorphism
                        memcpy(quadrangulationAutomorphisms[quadrangulationAutomorphismsCount], cagqAlternateLabelling, sizeof(int)*MAXN);
                        quadrangulationAutomorphismsCount++;
                    }
                }
                e = e->next;
            } while (e!=elast);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

int aaCertificate[MAXE+MAXN];
int aaAnglesCertificate[MAXE+MAXN];
int aaAlternateCertificate[MAXE+MAXN];
int aaAnglesAlternateCertificate[MAXE+MAXN];
int aaAlternateLabelling[MAXN];
EDGE *aaAlternateFirstedge[MAXN];
int aaQueue[MAXN];

void constructAlternateAngleAssignmentCertificate(EDGE *eStart){
    int i;
    for(i=0; i<MAXN; i++){
        aaAlternateLabelling[i] = MAXN;
    }
    EDGE *e, *elast;
    int head = 1;
    int tail = 0;
    int vertexCounter = 1;
    int aaAlternateCertificatePosition = 0;
    aaQueue[0] = eStart->start;
    aaAlternateFirstedge[eStart->start] = eStart;
    aaAlternateLabelling[eStart->start] = 0;
    while(head>tail){
        int currentVertex = aaQueue[tail++];
        e = elast = aaAlternateFirstedge[currentVertex];
        do {
            if(aaAlternateLabelling[e->end]==MAXN){
                aaQueue[head++] = e->end;
                aaAlternateLabelling[e->end] = vertexCounter++;
                aaAlternateFirstedge[e->end] = e->inverse;
            }
            aaAlternateCertificate[aaAlternateCertificatePosition] = aaAlternateLabelling[e->end];
            aaAnglesAlternateCertificate[aaAlternateCertificatePosition++] = e->angle;
            e = e->next;
        } while (e!=elast);
        aaAlternateCertificate[aaAlternateCertificatePosition] = MAXN;
        aaAnglesAlternateCertificate[aaAlternateCertificatePosition++] = MAXN;
    }
}

void constructAlternateAngleAssignmentCertificateOrientationReversing(EDGE *eStart){
    int i;
    for(i=0; i<MAXN; i++){
        aaAlternateLabelling[i] = MAXN;
    }
    EDGE *e, *elast;
    int head = 1;
    int tail = 0;
    int vertexCounter = 1;
    int aaAlternateCertificatePosition = 0;
    aaQueue[0] = eStart->start;
    aaAlternateFirstedge[eStart->start] = eStart;
    aaAlternateLabelling[eStart->start] = 0;
    while(head>tail){
        int currentVertex = aaQueue[tail++];
        e = elast = aaAlternateFirstedge[currentVertex];
        do {
            if(aaAlternateLabelling[e->end]==MAXN){
                aaQueue[head++] = e->end;
                aaAlternateLabelling[e->end] = vertexCounter++;
                aaAlternateFirstedge[e->end] = e->inverse;
            }
            aaAlternateCertificate[aaAlternateCertificatePosition] = aaAlternateLabelling[e->end];
            e = e->prev;
            aaAnglesAlternateCertificate[aaAlternateCertificatePosition++] = e->angle;
        } while (e!=elast);
        aaAlternateCertificate[aaAlternateCertificatePosition] = MAXN;
        aaAnglesAlternateCertificate[aaAlternateCertificatePosition++] = MAXN;
    }
}

/**
 * Checks whether the current angle assignment is canonical.
 */
boolean isCanonicalAngleAssignment(){
    //construct certificate
    int pos = 0;
    int i, j;
    
    for(i=0; i<nv; i++){
        EDGE *e, *elast;

        e = elast = firstedge[i];
        do {
            aaCertificate[pos] = e->end;
            aaAnglesCertificate[pos++] = e->angle;
            e = e->next;
        } while (e!=elast);
        aaCertificate[pos] = MAXN;
        aaAnglesCertificate[pos++] = MAXN;
    }
    
    //construct alternate certificates
    EDGE *ebase = firstedge[0];
    
    for(i=0; i<nv; i++){
        if(degree[i]==degree[0]){
            EDGE *e, *elast;

            e = elast = firstedge[i];
            do {
                if(e!=ebase){
                    constructAlternateAngleAssignmentCertificate(e);
                    //if the vertex certificates are equal, we compare the angle certificates
                    if(memcmp(aaCertificate, aaAlternateCertificate, sizeof(int)*pos) == 0) {
                        //compare angle certificates
                        for(j = 0; j < pos; j++){
                            if(aaAnglesCertificate[j] < aaAnglesAlternateCertificate[j]){
                                break;
                            } else if(aaAnglesCertificate[j] > aaAnglesAlternateCertificate[j]){
                                return FALSE;
                            }
                        }
                        //compare angle-reversed certificates
                        if(generateSTCQ4){
                            //when generating STCQ4 we can interchange alpha <-> gamma
                            //compare angle-reversed certificates
                            for(j = 0; j < pos; j++){
                                if(aaAnglesAlternateCertificate[j]==0){ //alpha
                                    if(aaAnglesCertificate[j] < 2){
                                        break;
                                    } else if(aaAnglesCertificate[j] > 2){
                                        return FALSE;
                                    }
                                } else if(aaAnglesAlternateCertificate[j]==2){ //gamma
                                    if(aaAnglesCertificate[j] > 0){
                                        return FALSE;
                                    }
                                } else { //beta or delta
                                    if(aaAnglesCertificate[j] < aaAnglesAlternateCertificate[j]){
                                        break;
                                    } else if(aaAnglesCertificate[j] > aaAnglesAlternateCertificate[j]){
                                        return FALSE;
                                    }
                                }
                            }
                        } else {
                            //when generating STCQ4 we can't interchange alpha <-> delta and beta <-> gamma
                            //compare angle-reversed certificates
                            for(j = 0; j < pos; j++){
                                if(aaAnglesCertificate[j] < 3 - aaAnglesAlternateCertificate[j]){
                                    break;
                                } else if(aaAnglesCertificate[j] > 3 - aaAnglesAlternateCertificate[j]){
                                    return FALSE;
                                }
                            }
                        }
                    }
                }
                if(!mirrorImagesAreDistinct){
                    constructAlternateAngleAssignmentCertificateOrientationReversing(e);
                    //if the vertex certificates are equal, we compare the angle certificates
                    if(memcmp(aaCertificate, aaAlternateCertificate, sizeof(int)*pos) == 0) {
                        //compare angle certificates
                        for(j = 0; j < pos; j++){
                            if(aaAnglesCertificate[j] < aaAnglesAlternateCertificate[j]){
                                break;
                            } else if(aaAnglesCertificate[j] > aaAnglesAlternateCertificate[j]){
                                return FALSE;
                            }
                        }
                        if(generateSTCQ4){
                            //when generating STCQ4 we can interchange alpha <-> gamma
                            //compare angle-reversed certificates
                            for(j = 0; j < pos; j++){
                                if(aaAnglesAlternateCertificate[j]==0){ //alpha
                                    if(aaAnglesCertificate[j] < 2){
                                        break;
                                    } else if(aaAnglesCertificate[j] > 2){
                                        return FALSE;
                                    }
                                } else if(aaAnglesAlternateCertificate[j]==2){ //gamma
                                    if(aaAnglesCertificate[j] > 0){
                                        return FALSE;
                                    }
                                } else { //beta or delta
                                    if(aaAnglesCertificate[j] < aaAnglesAlternateCertificate[j]){
                                        break;
                                    } else if(aaAnglesCertificate[j] > aaAnglesAlternateCertificate[j]){
                                        return FALSE;
                                    }
                                }
                            }
                        } else {
                            //when generating STCQ4 we can't interchange alpha <-> delta and beta <-> gamma
                            //compare angle-reversed certificates
                            for(j = 0; j < pos; j++){
                                if(aaAnglesCertificate[j] < 3 - aaAnglesAlternateCertificate[j]){
                                    break;
                                } else if(aaAnglesCertificate[j] > 3 - aaAnglesAlternateCertificate[j]){
                                    return FALSE;
                                }
                            }
                        }
                    }
                }
                e = e->next;
            } while (e!=elast);
        }
    }
    return TRUE;
}

void calculateAutomorphismGroupAngleAssignments(){
    aaAutomorphismsCount = 0;
    aaAutomorphismGroupContainsOrientationReversingSymmetry = FALSE;
    
    //construct certificate
    int pos = 0;
    int i;
    
    for(i=0; i<nv; i++){
        aaAutomorphisms[aaAutomorphismsCount][i] = i; //include identity
        
        EDGE *e, *elast;

        e = elast = firstedge[i];
        do {
            aaCertificate[pos] = e->end;
            aaAnglesCertificate[pos++] = e->angle;
            e = e->next;
        } while (e!=elast);
        aaCertificate[pos] = MAXN;
        aaAnglesCertificate[pos++] = MAXN;
    }
    
    aaAutomorphismsCount++;
    
    //construct alternate certificates
    EDGE *ebase = firstedge[0];
    
    for(i=0; i<nv; i++){
        if(degree[i]==degree[0]){
            EDGE *e, *elast;

            e = elast = firstedge[i];
            do {
                if(e!=ebase){
                    constructAlternateAngleAssignmentCertificate(e);
                    if(memcmp(aaCertificate, aaAlternateCertificate, sizeof(int)*pos) == 0 &&
                            memcmp(aaAnglesCertificate, aaAnglesAlternateCertificate, sizeof(int)*pos) == 0) {
                        //store automorphism
                        memcpy(aaAutomorphisms[aaAutomorphismsCount], aaAlternateLabelling, sizeof(int)*MAXN);
                        aaAutomorphismsCount++;
                    }
                }
                constructAlternateAngleAssignmentCertificateOrientationReversing(e);
                if(memcmp(aaCertificate, aaAlternateCertificate, sizeof(int)*pos) == 0 &&
                        memcmp(aaAnglesCertificate, aaAnglesAlternateCertificate, sizeof(int)*pos) == 0) {
                    //store automorphism
                    memcpy(aaAutomorphisms[aaAutomorphismsCount], aaAlternateLabelling, sizeof(int)*MAXN);
                    aaAutomorphismsCount++;
                    aaAutomorphismGroupContainsOrientationReversingSymmetry = TRUE;
                }
                e = e->next;
            } while (e!=elast);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void handleSolution(lprec *lp) {
    if(!isCanonicalAngleAssignment()) return;
    solvableAndCanonical++;
    if(outputSolution){
        if(outputFormat == 'h'){
            //human-readable output
            printSphericalTilingByCongruentQuadrangles(lp);
        } else if(outputFormat == 'c'){
            //code

        }
        if(latexPerSolution || latexSummaryFile!=NULL){
            //output to LaTeX
            printAngleAssignmentLatex();
        }
    }
}

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
    if(onlyConvex){
        if(generateSTCQ4){
            resize_lp(lp, nv + 3 - duplicateEquationCount, get_Ncolumns(lp));
        } else {
            resize_lp(lp, nv + 5 - duplicateEquationCount, get_Ncolumns(lp));
        }
        /* There nv + 1 equations: one for each vertex plus the extra equation
         * 
         *     alpha + beta + gamma + delta = 2 +4/F.
         * 
         * Of course we only add each distinct equation once, so we subtract the
         * number of duplicate equations.
         * For STCQ2 there are also 4 inequalities:
         *     alpha - beta + delta < 1
         *     alpha + beta - delta < 1
         *     alpha - gamma + delta < 1
         *    -alpha + gamma + delta < 1
         * For STCQ4 there are 2 inequalities:
         *     alpha - gamma + delta < 1
         *    -alpha + gamma + delta < 1
         */
    } else {
        resize_lp(lp, nv + 1 - duplicateEquationCount, get_Ncolumns(lp));
        /* There nv + 1 equations: one for each vertex plus the extra equation
         * 
         *     alpha + beta + gamma + delta = 2 +4/F.
         * 
         * Of course we only add each distinct equation once, so we subtract the
         * number of duplicate equations.
         */
    }

    //name the columns
    set_col_name(lp, 1, "alpha");
    set_col_name(lp, 2, "beta");
    set_col_name(lp, 3, "gamma");
    set_col_name(lp, 4, "delta");

    REAL epsilon = 0.000001;
    REAL lowerBoundAngle = 0 + epsilon;
    REAL upperBoundAngle = 2 - epsilon;
    if(onlyConvex){
        upperBoundAngle = 1 - epsilon;
    }
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
    
    if(onlyConvex){

        if(!generateSTCQ4){
            //the following constraints are only applicable for STCQ2
            colno[0] = 1;
            row[0] = 1;
            colno[1] = 2;
            row[1] = -1;
            colno[2] = 3;
            row[2] = 0;
            colno[3] = 4;
            row[3] = 1;

            //in case of onlyConvex: upperBoundAngle = 1 - epsilon
            if (!add_constraintex(lp, 4, row, colno, LE, upperBoundAngle)) {
                exit(1);
            }


            colno[0] = 1;
            row[0] = 1;
            colno[1] = 2;
            row[1] = 1;
            colno[2] = 3;
            row[2] = 0;
            colno[3] = 4;
            row[3] = -1;

            //in case of onlyConvex: upperBoundAngle = 1 - epsilon
            if (!add_constraintex(lp, 4, row, colno, LE, upperBoundAngle)) {
                exit(1);
            }
        }

        colno[0] = 1;
        row[0] = 1;
        colno[1] = 2;
        row[1] = 0;
        colno[2] = 3;
        row[2] = -1;
        colno[3] = 4;
        row[3] = 1;
        
        //in case of onlyConvex: upperBoundAngle = 1 - epsilon
        if (!add_constraintex(lp, 4, row, colno, LE, upperBoundAngle)) {
            exit(1);
        }

        colno[0] = 1;
        row[0] = -1;
        colno[1] = 2;
        row[1] = 0;
        colno[2] = 3;
        row[2] = 1;
        colno[3] = 4;
        row[3] = 1;
        
        //in case of onlyConvex: upperBoundAngle = 1 - epsilon
        if (!add_constraintex(lp, 4, row, colno, LE, upperBoundAngle)) {
            exit(1);
        }
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
        handleSolution(lp);
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
        EDGE *e1 = matchingEdges[orderedFaces[i]];
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
void createPartialSystem(int currentFaceCount) {
    //clear systems
    int i;
    for (i = 0; i < nv; i++) {
        alphaCount[i] = betaCount[i] = gammaCount[i] = deltaCount[i] = 0;
    }

    //iterate over all faces
    for (i = 0; i < currentFaceCount; i++) {
        EDGE *e1 = matchingEdges[orderedFaces[i]];
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

boolean firstCheckOfSystem() {
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
         * at Hamming distance 1 of that earlier equation.
         */
        for (j = i + 1; j < nv; j++) {
            if (isDuplicateEquation[j]) continue;
            /*
             * this equation is already a duplicate, so any equation
             * that would be at Hamming distance 1 of this equation is also
             * at Hamming distance 1 of that earlier equation.
             */
            int diffAlpha = alphaCount[i] - alphaCount[j];
            int diffBeta = betaCount[i] - betaCount[j];
            int diffGamma = gammaCount[i] - gammaCount[j];
            int diffDelta = deltaCount[i] - deltaCount[j];
            
            int hammingDistance = 0;
            if (diffAlpha != 0) hammingDistance++;
            if (diffBeta != 0) hammingDistance++;
            if (diffGamma != 0) hammingDistance++;
            if (diffDelta != 0) hammingDistance++;
            if (hammingDistance == 0) {
                // if we get here, then equation j is a duplicate of i
                isDuplicateEquation[j] = TRUE;
                duplicateEquationCount++;
            } else {
                if((diffAlpha<=0 && diffBeta<=0 && diffGamma<=0 && diffDelta<=0) || (diffAlpha>=0 && diffBeta>=0 && diffGamma>=0 && diffDelta>=0)){
                    return FALSE;
                }
                if(diffAlpha==-diffGamma && diffBeta==0 && diffDelta==0){
                    //alpha==gamma
                    return FALSE;
                }
                if(diffAlpha==-diffDelta && diffGamma==0 && diffBeta==0){
                    //alpha==delta
                    return FALSE;
                }
                if(diffBeta==-diffGamma && diffAlpha==0 && diffDelta==0){
                    //beta==gamma
                    return FALSE;
                }
                if(diffBeta==-diffDelta && diffGamma==0 && diffAlpha==0){
                    //beta==delta
                    return FALSE;
                }
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
        rejectedByCoefficientDiff++;
        if (printUnsolvableSystems || writeHammingDistanceUnsolvedSystems) {
            simplifySystem();
            printSystem();
        }
    }
#ifdef _DEBUG
    printSystem();
#endif
}

boolean checkPartialSystem(int currentFace) {
    int i, j;

    //reset array and count number of degree 3 types
    int degreeThreeVertexTypeCount = 0;
    int type1, type2;
    for (i = 0; i < nv; i++) {
        isDuplicateEquation[i] = FALSE;
        if(degree[i]==3 && vertexCompletedAfterFace[i]<currentFace){
            int currentType = degreeThreeVertexTypes[i] = getDegreeThreeVertexType(i);
            if(degreeThreeVertexTypeCount == 0){
                type1 = currentType;
                degreeThreeVertexTypeCount = 1;
            } else if(degreeThreeVertexTypeCount == 1){
                if(type1 != currentType){
                    type2 = currentType;
                    degreeThreeVertexTypeCount = 2;
                }
            } else { //(degreeThreeVertexTypeCount == 2)
                if(type1 != currentType && type2 != currentType){
                    //at most two different types of degree three vertices
                    return FALSE;
                }
            }
        } else {
            degreeThreeVertexTypes[i] = -1;
        }
    }
    duplicateEquationCount = 0;
    
    if(degreeThreeVertexTypeCount == 2){
        if(!degreeThreeTypesCompatibility[type1][type2]){
            return FALSE;
        }
        if(degreeThreeTypesCombinationVertexLowerBound[type1][type2]>nv){
            return FALSE;
        }
        if(degreeThreeTypesCombinationVertexUpperBound[type1][type2]<nv){
            return FALSE;
        }
    }

    for (i = 0; i < nv - 1; i++) {
        if (vertexCompletedAfterFace[i]>=currentFace) continue;
        if (isDuplicateEquation[i]) continue;
        /*
         * this equation is already a duplicate itself, so any equation
         * that would be at Hamming distance 1 of this equation is also
         * at Hamming Distance of that earlier equation.
         */
        for (j = i + 1; j < nv; j++) {
            if (vertexCompletedAfterFace[j]>=currentFace) continue;
            if (isDuplicateEquation[j]) continue;
            /*
             * this equation is already a duplicate, so any equation
             * that would be at Hamming distance 1 of this equation is also
             * at Hamming Distance of that earlier equation.
             */
            
            int diffAlpha = alphaCount[i] - alphaCount[j];
            int diffBeta = betaCount[i] - betaCount[j];
            int diffGamma = gammaCount[i] - gammaCount[j];
            int diffDelta = deltaCount[i] - deltaCount[j];
            
            int hammingDistance = 0;
            if (diffAlpha != 0) hammingDistance++;
            if (diffBeta != 0) hammingDistance++;
            if (diffGamma != 0) hammingDistance++;
            if (diffDelta != 0) hammingDistance++;
            if (hammingDistance == 0) {
                // if we get here, then equation j is a duplicate of i
                isDuplicateEquation[j] = TRUE;
                duplicateEquationCount++;
            } else {
                if((diffAlpha<=0 && diffBeta<=0 && diffGamma<=0 && diffDelta<=0) || (diffAlpha>=0 && diffBeta>=0 && diffGamma>=0 && diffDelta>=0)){
                    return FALSE;
                }
                if(diffAlpha==-diffGamma && diffBeta==0 && diffDelta==0){
                    //alpha==gamma
                    return FALSE;
                }
                if(diffAlpha==-diffDelta && diffGamma==0 && diffBeta==0){
                    //alpha==delta
                    return FALSE;
                }
                if(diffBeta==-diffGamma && diffAlpha==0 && diffDelta==0){
                    //beta==gamma
                    return FALSE;
                }
                if(diffBeta==-diffDelta && diffGamma==0 && diffAlpha==0){
                    //beta==delta
                    return FALSE;
                }
            }
        }
    }
    return TRUE;
}

boolean checkSTCQ4Assignment(int currentFace){
    //TODO: avoid checking the same face multiple times
    //TODO: check efficiency of this method!!
    
    //check that a c-edge is always next to a c-edge
    int i;
    
    for(i = 0; i < currentFace; i++){
        EDGE *em = matchingEdges[orderedFaces[i]]; //matching edge
        int fn;
        if(angleAssigmentDirection[i]){
            fn = em->next->rightface; // neighbouring face along c-edge
        } else {
            fn = em->inverse->prev->inverse->rightface; // neighbouring face along c-edge
        }
        if(faceRank[fn] < currentFace && faceRank[fn] < i){
            //if angles are fixed for the neighbouring face
            //second part is to avoid comparing the same pair twice

            EDGE *efn = matchingEdges[fn];
            //find face next to neighbouring face
            if(angleAssigmentDirection[faceRank[fn]]){
                if(orderedFaces[i] != efn->next->rightface){
                    return FALSE;
                }
            } else {
                if(orderedFaces[i] != efn->inverse->prev->inverse->rightface){
                    return FALSE;
                }
            }
        }
    }
    //no violation of STCQ4 edge assignment found
    return TRUE;
}

void assignAnglesForCurrentPerfectMatchingRecursion(int currentFace) {
    if(generateSTCQ4 && !checkSTCQ4Assignment(currentFace)){
        return;
    }
    if (currentFace == nv - 2) {
        handleAngleAssignment();
    } else {
        if(boundAngleAssignments && checkVerticesAfterFace[currentFace]){
            createPartialSystem(currentFace);
            if(!checkPartialSystem(currentFace)){
                return;
            }
        }
        angleAssigmentDirection[currentFace] = 0;
        assignAnglesForCurrentPerfectMatchingRecursion(currentFace + 1);
        angleAssigmentDirection[currentFace] = 1;
        assignAnglesForCurrentPerfectMatchingRecursion(currentFace + 1);
    }
}

void assignAnglesForCurrentPerfectMatching() {
    angleAssigmentDirection[0] = 0;
    assignAnglesForCurrentPerfectMatchingRecursion(1);
    angleAssigmentDirection[0] = 1;
    assignAnglesForCurrentPerfectMatchingRecursion(1);
}

int matchingCount = 0;

void handlePerfectMatching() {
    matchingCount++;
    assignAnglesForCurrentPerfectMatching();
}

void matchNextFace(int lastFace, int matchingSize) {
    //handle splitting
    if(matchingSize == splitting_level){
        
        if (splitting_count-- != 0){
            return;
        }
        splitting_count = splitting_mod - 1;
    }
        
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
        if (!matched[neighbour] && e->allowedInFaceMatching) {
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

void markEdgesAtCubicTristar(){
    int i;
    for(i=0; i<nv; i++){
        if(degree[i]==3){
            EDGE *e, *elast;

            int isTristar = TRUE;

            e = elast = firstedge[i];
            do {
                if(degree[e->end]!=3){
                    isTristar = FALSE;
                }
                e = e->next;
            } while (e!=elast);

            if(isTristar){
                e = elast = firstedge[i];
                do {
                    e->allowedInFaceMatching = FALSE;
                    e->inverse->prev->inverse->next->allowedInFaceMatching = FALSE;
                    e->inverse->next->inverse->prev->allowedInFaceMatching = FALSE;
                    e = e->next;
                } while (e!=elast);
            }
        }
    }
}

int generate_perfect_matchings_in_dual() {
    int i;

    matchingCount = 0;

    unsigned long long int oldSolutionCount = solvable;

    if (nf != nv - 2) {
        fprintf(stderr, "Something went horribly wrong. Maybe some wrong parameter?\nnf: %d, nv: %d\n", nf, nv);
        exit(1);
    }
    
    if(!generateAllMatchings){
        if(nv>8 && onlyConvex){
            markEdgesAtCubicTristar();
        }
    }

    for (i = 0; i < nv - 2; i++) {
        matched[i] = FALSE;
    }
    matched[0] = TRUE;

    EDGE *e, *elast;

    e = elast = facestart[0];
    do {
        if(e->allowedInFaceMatching){
            int neighbour = e->inverse->rightface;
            match[0] = neighbour;
            match[neighbour] = 0;
            matched[neighbour] = TRUE;

            //store edge corresponding to the match for each face.
            matchingEdges[0] = e;
            matchingEdges[neighbour] = e->inverse;

            matchNextFace(0, 1);

            matched[neighbour] = FALSE;
        }
        e = e->inverse->prev;
    } while (e != elast);


    perfect_matchings_counts = increment(perfect_matchings_counts, matchingCount);

    if (oldSolutionCount == solvable) {
        unusedGraphCount++;
        if (unusedQuadrangulations) {
            outputQuadrangulation();
        }
    } else if (usedQuadrangulations) {
        outputQuadrangulation();
    }

    return 0;
}

void printSummary() {
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
    if(filterOnly){
        fprintf(stderr, "Only quadrangulation %llu was used\n", filterOnly);
    }
    fprintf(stderr, "\nMatchings: %llu\n", totalPerfectMatchingsCount);
    fprintf(stderr, "\nAssignments: %llu\n", assignmentCount);
    fprintf(stderr, "\nSolvable: %llu\n", solvable);
    fprintf(stderr, "\nSolvable and canonical: %llu\n", solvableAndCanonical);
    if (printStatistics) {
        fprintf(stderr, "\nNon-solvable: %llu\n", assignmentCount - solvable);
        fprintf(stderr, "\n%llu quadrangulations do not correspond to a tiling.\n", unusedGraphCount);
        fprintf(stderr, "%llu quadrangulations can correspond to a tiling.\n", numberOfQuadrangulations - unusedGraphCount);
        fprintf(stderr, "\nRejected by coefficient diff: %llu\n", rejectedByCoefficientDiff);
        fprintf(stderr, "Rejected by lpsolve: %llu\n\n", assignmentCount - solvable - rejectedByCoefficientDiff);
    }
}

//////////////////////////////////////////////////////////////////////////////

void orderFaces(){
    int numberedFacesAt[MAXN];
    int isNumbered[MAXF];
    int faceCounter = 0;
    EDGE *edge, *edgelast;
    
    int i;
    
    for(i=0; i<MAXN; i++){
        numberedFacesAt[i] = 0;
    }
    for(i=0; i<MAXF; i++){
        isNumbered[i] = FALSE;
    }
    
    //number faces at cubic edges
    for(i=0; i<ne; i++){
        EDGE *e = edges+i;
        if(degree[e->end]==3 && degree[e->start]==3){
            int face1 = e->rightface;
            int face2 = e->inverse->rightface;
            //number faces at both side of the edge
            if(!isNumbered[face1]){
                //number face
                isNumbered[face1] = TRUE;
                orderedFaces[faceCounter] = face1;
                
                //mark changes for vertices around face
                edge = edgelast = facestart[face1];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);
                
                //increase faceCounter
                faceCounter++;
            }
            if(!isNumbered[face2]){
                //number face
                isNumbered[face2] = TRUE;
                orderedFaces[faceCounter] = face2;
                
                //mark changes for vertices around face
                edge = edgelast = facestart[face2];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);
                
                //increase faceCounter
                faceCounter++;
            }
            //number remaining faces incident to the end points
            int face3 = e->next->rightface;
            int face4 = e->inverse->next->rightface;
            if(!isNumbered[face3]){
                //number face
                isNumbered[face3] = TRUE;
                orderedFaces[faceCounter] = face3;
                
                //mark changes for vertices around face
                edge = edgelast = facestart[face3];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);
                
                //increase faceCounter
                faceCounter++;
            }
            if(!isNumbered[face4]){
                //number face
                isNumbered[face4] = TRUE;
                orderedFaces[faceCounter] = face4;
                
                //mark changes for vertices around face
                edge = edgelast = facestart[face4];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);
                
                //increase faceCounter
                faceCounter++;
            }
        }
    }
    
    //number faces at cubic vertices
    for(i=0; i<nv; i++){
        if(degree[i]==3){
            EDGE *e = firstedge[i];
            //number faces incident to vertex
            int face1 = e->rightface;
            int face2 = e->next->rightface;
            int face3 = e->next->next->rightface;
            //number faces at both side of the edge
            if(!isNumbered[face1]){
                //number face
                isNumbered[face1] = TRUE;
                orderedFaces[faceCounter] = face1;

                //mark changes for vertices around face
                edge = edgelast = facestart[face1];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);

                //increase faceCounter
                faceCounter++;
            }
            if(!isNumbered[face2]){
                //number face
                isNumbered[face2] = TRUE;
                orderedFaces[faceCounter] = face2;

                //mark changes for vertices around face
                edge = edgelast = facestart[face2];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);

                //increase faceCounter
                faceCounter++;
            }
            if(!isNumbered[face3]){
                //number face
                isNumbered[face3] = TRUE;
                orderedFaces[faceCounter] = face3;

                //mark changes for vertices around face
                edge = edgelast = facestart[face3];
                do {
                    numberedFacesAt[edge->start]++;
                    if(numberedFacesAt[edge->start]==degree[edge->start]){
                        //vertex completed
                        checkVerticesAfterFace[faceCounter] = TRUE;
                        vertexCompletedAfterFace[edge->start] = faceCounter;
                    }
                    edge = edge->inverse->prev;
                } while (edge!=edgelast);

                //increase faceCounter
                faceCounter++;
            }
        }
    }
    
    //number remaining faces
    for(i=0; i<nv; i++){
        if(degree[i]!=numberedFacesAt[i]){
            EDGE *e, *elast;
            e = elast = firstedge[i];
            do {
                //number faces incident to vertex
                int face = e->rightface;
                //number faces at both side of the edge
                if(!isNumbered[face]){
                    //number face
                    isNumbered[face] = TRUE;
                    orderedFaces[faceCounter] = face;

                    //mark changes for vertices around face
                    edge = edgelast = facestart[face];
                    do {
                        numberedFacesAt[edge->start]++;
                        if(numberedFacesAt[edge->start]==degree[edge->start]){
                            //vertex completed
                            checkVerticesAfterFace[faceCounter] = TRUE;
                            vertexCompletedAfterFace[edge->start] = faceCounter;
                        }
                        edge = edge->inverse->prev;
                    } while (edge!=edgelast);

                    //increase faceCounter
                    faceCounter++;
                }
                e = e->next;
            } while(e!=elast);
        }
    }
    
    //construct faceRank
    for(i = 0; i < nf; i++){
        faceRank[orderedFaces[i]] = i;
    }
}

//////////////////////////////////////////////////////////////////////////////

/*
 * This method returns TRUE if the current quadrangulation contains a cubic
 * quadrangle, i.e., a quadrangle with 4 cubic vertices, and FALSE otherwise.
 */
boolean cubicQuadSearch(){
    int i;
    for(i=0; i<nf; i++){
        EDGE *e, *elast;
        
        boolean isCubicQuad = TRUE;
        
        e = elast = facestart[i];
        do {
            if(degree[e->start]!=3){
                isCubicQuad = FALSE;
            }
            e = e->inverse->prev;
        } while (e!=elast);
        
        if(isCubicQuad) return TRUE;
    }
    return FALSE;
}

/*
 * Method that returns FALSE if can be decided that this quadrangulation does not
 * admit a STCQ2. At the moment this is only the case when the quadrangulation
 * contains a cubic quadrangle, i.e., a quadrangle with 4 cubic vertices.
 */
boolean earlyFilterQuadrangulations(){
    if(!onlyConvex){
        //we currently have no early filter for the concave case
        return TRUE;
    }
    if(nf>6 && cubicQuadSearch()){
        return FALSE;
    }
    return TRUE;
}

/* Sets the labelling of the quadrangulation to a BFS-labelling
 */
void relabelQuadrangulation(){
    int newLabelling[MAXN];
    int newDegree[MAXN];
    EDGE *newFirstEdge[MAXN];
    int queue[MAXN];
    int i;
    for(i=0; i<MAXN; i++){
        newLabelling[i] = MAXN;
    }
    EDGE *e, *elast;
    int head = 1;
    int tail = 0;
    int vertexCounter = 1;
    queue[0] = 0;
    newFirstEdge[0] = firstedge[0];
    newLabelling[0] = 0;
    newDegree[0] = degree[0];
    while(head>tail){
        int currentVertex = queue[tail++];
        e = elast = newFirstEdge[currentVertex];
        do {
            if(newLabelling[e->end]==MAXN){
                //unlabelled vertex gets label vertexCounter
                newLabelling[e->end] = vertexCounter;
                //newly labelled vertex is placed on queue and initial edge is stored
                queue[head++] = vertexCounter;
                newFirstEdge[vertexCounter] = e->inverse;
                //the degree is copied
                newDegree[vertexCounter] = degree[e->end];
                //increase counter
                vertexCounter++;
            }
            e = e->next;
        } while (e!=elast);
    }
    
    //apply new labelling
    for(i = 0; i < nv; i++){
        e = elast = firstedge[i] = newFirstEdge[i];
        degree[i] = newDegree[i];
        do {
            e->start = i;
            e->end = newLabelling[e->end];
            e = e->next;
        } while (e!=elast);
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
        edges[edgeCounter].allowedInFaceMatching = TRUE;
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
            edges[edgeCounter].allowedInFaceMatching = TRUE;
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
    
    ne = edgeCounter;
    
    makeDual();
    
    // nv - ne/2 + nf = 2
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
    
    int readCount;


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
        readCount = fread(code, sizeof (unsigned short), 1, file);
        if(!readCount){
            fprintf(stderr, "Unexpected EOF.\n");
            exit(1);
        }
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant N too small %d > %d \n", code[0], MAXN);
            exit(1);
        }
        bufferSize = 1;
        zeroCounter = 0;
        while (zeroCounter < code[0]) {
            readCount = fread(code + bufferSize, sizeof (unsigned short), 1, file);
            if(!readCount){
                fprintf(stderr, "Unexpected EOF.\n");
                exit(1);
            }
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    }

    *length = bufferSize;
    return (1);


}

//====================== USAGE =======================

void help(char *name){
    fprintf(stderr, "The program %s calculates spherical tilings by congruent qaudrangles\n", name);
    fprintf(stderr, "that have any of the input graphs as underlying graph.\n\n");
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "Without any options, this program will generate spherical tilings by\n");
    fprintf(stderr, "congruent convex without writing any output.\n\n");
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices. Recompile if you need larger\n", MAXN);
    fprintf(stderr, "graphs.\n\n");
    fprintf(stderr, "General options\n===============\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
    fprintf(stderr, "    -4\n");
    fprintf(stderr, "       Generate STCQ4 tilings instead of STCQ2. (experimental)\n");
    fprintf(stderr, "    -t, --type type\n");
    fprintf(stderr, "       Specifies the generated type where type is one of\n");
    fprintf(stderr, "           e, edge    edge assignments\n");
    fprintf(stderr, "           a, angle   angle assignments\n");
    fprintf(stderr, "           t, tiling  spherical tilings (default)\n");
    fprintf(stderr, "    -c, --concave\n");
    fprintf(stderr, "       Also allow concave quadrangles (currently not supported)\n");
    fprintf(stderr, "    -s, --statistics\n");
    fprintf(stderr, "       Print extra statistics\n");
    fprintf(stderr, "    -f, --filter number\n");
    fprintf(stderr, "       Only perform the calculations for the graph with the given number.\n");
    fprintf(stderr, "    -r, --relabel\n");
    fprintf(stderr, "       Relabel the quadrangulations that are used as input. The program requires\n");
    fprintf(stderr, "       the graphs to have a BFS-labelling compatible with the embedding. If the\n");
    fprintf(stderr, "       input comes from plantri, then relabelling is not necessary.\n");
    fprintf(stderr, "    --mirror\n");
    fprintf(stderr, "       Makes the program consider mirror images as distinct.\n");
    fprintf(stderr, "    -m, --modulo r:n\n");
    fprintf(stderr, "       Split the generation in multiple parts. The generation is split into n\n");
    fprintf(stderr, "       parts and only part r is generated. The number n needs to be an integer\n");
    fprintf(stderr, "       larger than 0 and r should be a positive integer smaller than n.\n");
    fprintf(stderr, "    --splitlevel l\n");
    fprintf(stderr, "       Sets the level at which point the generation will be split. By default\n");
    fprintf(stderr, "       this is set to 10. The value l should lie between 1 and half the number\n");
    fprintf(stderr, "       of faces in the quadrangulation. If l is less than or equal to 0, then\n");
    fprintf(stderr, "       the splitting will be disabled. If l is more than half the number of\n");
    fprintf(stderr, "       faces, then all work will be done in part 0.\n");
    fprintf(stderr, "\nOutput options\n==============\n");
    fprintf(stderr, "    -o, --output format\n");
    fprintf(stderr, "       Specifies the export format where format is one of\n");
    fprintf(stderr, "           c, code    code depends on the generated type\n");
    fprintf(stderr, "           h, human   human-readable output\n");
    fprintf(stderr, "           n, none    no output: only count (default)\n");
    fprintf(stderr, "    --usedquadrangulations\n");
    fprintf(stderr, "       Only output quadrangulations that might be used in a STCQ.\n");
    fprintf(stderr, "    --unusedquadrangulations\n");
    fprintf(stderr, "       Only output quadrangulations that cannot be used in a STCQ.\n");
    fprintf(stderr, "    --latex filename\n");
    fprintf(stderr, "       Writes the solutions to a file with the given name as a LaTeX fragment.\n");
    fprintf(stderr, "       Note that this option cancels any previous --latex-per-solution.\n");
    fprintf(stderr, "    --latex-per-solution basename\n");
    fprintf(stderr, "       Writes the solutions to files formed with the given base name. Each file\n");
    fprintf(stderr, "       contains a solution as a LaTeX fragment. basename is given as a format.\n");
    fprintf(stderr, "       A valid use is, e.g.:\n\n");
    fprintf(stderr, "            --latex-per-solution result_%%02d.tex\n\n");
    fprintf(stderr, "       If basename contains no format tag, then all solution will be written to\n");
    fprintf(stderr, "       the same file, and only the last solution will be present in the file.\n");
    fprintf(stderr, "       Note that this option cancels any previous --latex.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options] n\n", name);
    fprintf(stderr, "       %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]){
    
    int i;
    for(i = 0; i < argc; i++){
        fprintf(stderr, "%s ", argv[i]);
    }
    fprintf(stderr, "\n");

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
        {"usedquadrangulations", no_argument, &usedQuadrangulations, TRUE},
        {"unusedquadrangulations", no_argument, &unusedQuadrangulations, TRUE},
        {"latex", required_argument, NULL, 0},
        {"onebased", no_argument, &oneBased, TRUE},
        {"group", no_argument, &includeGroup, TRUE},
        {"latex-per-solution", required_argument, NULL, 0},
        {"mirror", no_argument, NULL, 0},
        {"splitlevel", required_argument, NULL, 0},
        {"help", no_argument, NULL, 'h'},
        {"concave", no_argument, NULL, 'c'},
        {"statistics", no_argument, NULL, 's'},
        {"type", required_argument, NULL, 't'},
        {"output", required_argument, NULL, 'o'},
        {"filter", required_argument, NULL, 'f'},
        {"relabel", no_argument, NULL, 'r'},
        {"modulo", required_argument, NULL, 'm'}
    };
    int option_index = 0;

    char *splitting_string;
    while ((c = getopt_long(argc, argv, "hcst:o:f:4rm:", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                switch (option_index) {
                    case 0:
                    case 1:
                        //the correct flags are already set to true
                        //we just need to disable the default output
                        outputSolution = FALSE;
                        break;
                    case 2:
                        latexSummaryFile = fopen(optarg, "w");
                        latexPerSolution = FALSE;
                        break;
                    case 3:
                    case 4:
                        break;
                    case 5:
                        latexPerSolution = TRUE;
                        latexBaseName = optarg;
                        if(latexSummaryFile != NULL){
                            fclose(latexSummaryFile);
                        }
                        break;
                    case 6:
                        mirrorImagesAreDistinct = TRUE;
                        break;
                    case 7:
                        splitting_level =  atoi(optarg);
                        break;
                    default:
                        fprintf(stderr, "Illegal option.\n");
                        usage(name);
                        return EXIT_FAILURE;
                }
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'c':
                onlyConvex = FALSE;
                boundAngleAssignments = FALSE; //the known bounds only apply to convex stcq
                break;
            case 's':
                printStatistics = TRUE;
                break;
            case '4':
                generateSTCQ4 = TRUE;
                break;
            case 'o':
                outputFormat = optarg[0];
                switch (outputFormat) {
                    case 'n': //no output (default)
                    case 'c': //computer code
                    case 'h': //human-readable
                        break;
                    default:
                        fprintf(stderr, "Illegal output format %c.\n", c);
                        usage(name);
                        return 1;
                }
                break;
            case 't':
                generatedType = optarg[0];
                switch (generatedType) {
                    case 'e': //edge assignments
                    case 'a': //angle assignments
                    case 't': //spherical tilings
                        break;
                    default:
                        fprintf(stderr, "Illegal generated type %c.\n", c);
                        usage(name);
                        return 1;
                }
                break;
            case 'f':
                filterOnly = atoi(optarg);
                break;
            case 'r':
                relabelInputQuadrangulation = TRUE;
                break;
            case 'm':
                //modulo
                splitting_string = optarg;
                splitting_res = atoi(splitting_string);
                splitting_string = strchr(splitting_string, ':');
                if(splitting_string==NULL){
                    fprintf(stderr, "Illegal format for modulo.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                splitting_mod = atoi(splitting_string+1);
                if (splitting_res >= splitting_mod) {
                    fprintf(stderr, "Illegal format for modulo: rest must be smaller than mod.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                if (splitting_res < 0) {
                    fprintf(stderr, "Illegal format for modulo: rest must be positive.\n");
                    usage(name);
                    return EXIT_FAILURE;
                }
                break;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }
    
    //check splitting
    if(splitting_mod == 1){
        splitting_level = 0;
    }
    if(splitting_level <= 0 && splitting_res > 0){
        fprintf(stderr, "Nothing to do for this part.\n");
        return EXIT_SUCCESS;
    }

    /*=========== read quadrangulations ===========*/
    
    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        //reset splitting
        splitting_count = splitting_res;
        
        //decode the graph
        decodePlanarCode(code);
        
        //check splitting
        if(splitting_level > (nv - 2) / 2 && splitting_res > 0){
            //nothing to do for this part
            numberOfQuadrangulations++;
            continue;
        }
        
        if(relabelInputQuadrangulation){
            relabelQuadrangulation();
        }
        numberOfQuadrangulations++;
        if(filterOnly==0 || numberOfQuadrangulations==filterOnly){
            if(!isEarlyFilteringEnabled || earlyFilterQuadrangulations()){
                orderFaces();
                generate_perfect_matchings_in_dual(); 
            } else if(unusedQuadrangulations){
                unusedGraphCount++;
                outputQuadrangulation();
            } else {
                unusedGraphCount++;
            }
        }
    }
    //close any possible open files
    if(latexSummaryFile != NULL){
        fclose(latexSummaryFile);
    }
    printSummary();

}
