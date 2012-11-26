/* This program reads quadrangulations from standard in and
* calculates perfect matchings in the dual of the quadrangulation.
*
*
* Compile with:
*
* cc -o facematching -O4 facematching.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#ifndef MAXN
#define MAXN 64 /* the maximum number of vertices */
#endif
#define MAXE (4*MAXN-8) /* the maximum number of oriented edges */
#define MAXF (MAXN-2) /* the maximum number of faces */
#define MAXVAL (MAXN-2)/2 /* the maximum degree of a vertex */
#define MAXCODELENGTH (MAXN+MAXE+3)

#undef FALSE
#undef TRUE
#define FALSE 0
#define TRUE 1

typedef struct e /* The data type used for edges */ {
    int start; /* vertex where the edge starts */
    int end; /* vertex where the edge ends */
    int rightface; /* face on the right side of the edge
note: only valid if make_dual() called */
    struct e *prev; /* previous edge in clockwise direction */
    struct e *next; /* next edge in clockwise direction */
    struct e *inverse; /* the edge that is inverse to this one */
    int mark,index; /* two ints for temporary use;
Only access mark via the MARK macros. */

    int left_facesize; /* size of the face in prev-direction of the edge.
Only used for -p option. */
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


unsigned long long int numberOfGraphs = 0;

int matched[MAXF];
int match[MAXF];
EDGE *matchingEdges[MAXF];

int nv; //the number of vertices of the current quadrangulation
int nf; //the number of faces of the current quadrangulation

int automorphisms[4*MAXE][MAXN]; //there are at most 4e automorphisms
int automorphismsCount;

int calculateAllPerfectMatchings = FALSE;
int outputPerfectMatchings = FALSE;

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

void printAutomorphismGroup(){
    int i, j, printed[MAXN], next;
    fprintf(stderr, "automorphism count: %d\n", automorphismsCount);
    for (i = 0; i < automorphismsCount; i++) {
        for(j=0; j < MAXN; j++){
            printed[j] = 0;
        }
        for(j=0; j < nv; j++){
            if(!printed[j]){
                fprintf(stderr, "(%d", j);
                printed[j] = 1;
                next = automorphisms[i][j];
                while(!printed[next]){
                    fprintf(stderr, " %d", next);
                    printed[next] = 1;
                    next = automorphisms[i][next];
                }
                fprintf(stderr, ") ");
            }
        }
        fprintf(stderr, "\n");
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

void calculateAutomorphismGroup(){
    automorphismsCount = 0;
    
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
                        memcpy(automorphisms[automorphismsCount], cagqAlternateLabelling, sizeof(int)*MAXN);
                        automorphismsCount++;
                    }
                }
                constructAlternateCertificateOrientationReversing(e);
                if(memcmp(cagqCertificate, cagqAlternateCertificate, sizeof(int)*pos) == 0) {
                    //store automorphism
                    memcpy(automorphisms[automorphismsCount], cagqAlternateLabelling, sizeof(int)*MAXN);
                    automorphismsCount++;
                }
                e = e->next;
            } while (e!=elast);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////


long long int matchingCount = 0;

int isCanonicalMatching() {
    if(automorphismsCount==0) return 1;
    
    int currentMatching[nf/2][2];
    int currentMatchingCount = 0;
    
    int i, j;
    for (i = 0; i < nv - 2; i++) {
        EDGE *e = matchingEdges[i];
        if(e->start < e->end){
            int pos = currentMatchingCount;
            while(pos>0 && (currentMatching[pos-1][0]>e->start || (currentMatching[pos-1][0]==e->start && currentMatching[pos-1][1]>e->end))){
                currentMatching[pos][0] = currentMatching[pos-1][0];
                currentMatching[pos][1] = currentMatching[pos-1][1];
                pos--;
            }
            currentMatching[pos][0] = e->start;
            currentMatching[pos][1] = e->end;
            currentMatchingCount++;
        }
    }
    
    int alternateMatching[nf/2][2];
    
    for(i=0; i<automorphismsCount; i++){
        for(j=0; j<nf/2; j++){
            int newStart = automorphisms[i][currentMatching[j][0]];
            int newEnd = automorphisms[i][currentMatching[j][1]];
            if(newEnd<newStart){
                int t = newStart;
                newStart = newEnd;
                newEnd = t;
            }
            int pos = j;
            while(pos>0 && (alternateMatching[pos-1][0]>newStart || (alternateMatching[pos-1][0]==newStart && alternateMatching[pos-1][1]>newEnd))){
                alternateMatching[pos][0] = alternateMatching[pos-1][0];
                alternateMatching[pos][1] = alternateMatching[pos-1][1];
                pos--;
            }
            alternateMatching[pos][0] = newStart;
            alternateMatching[pos][1] = newEnd;
        }
        //compare matchings
        j = 0;
        while(j < nf/2 && alternateMatching[j][0] == currentMatching[j][0] && alternateMatching[j][1] == currentMatching[j][1])
            j++;
        if(j<nf/2 && (alternateMatching[j][0] < currentMatching[j][0] || (alternateMatching[j][0] == currentMatching[j][0] && alternateMatching[j][1] < currentMatching[j][1]))){
            return 0;
        }
    }
    return 1;
}

void handlePerfectMatching() {
    if(calculateAllPerfectMatchings || isCanonicalMatching()){
        matchingCount++;
    }
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

void generate_perfect_matchings_in_dual() {
    int i;

    matchingCount = 0;

    if (nf != nv - 2) {
        fprintf(stderr, "Something went horribly wrong. Maybe some wrong parameter?\nnf: %d, nv: %d\n", nf, nv);
        exit(1);
    }

    if(nv%2==0){
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
    }
    perfect_matchings_counts = increment(perfect_matchings_counts, matchingCount);
}

void perfect_matchings_summary() {
    unsigned long long int totalPerfectMatchingsCount = 0;
    if(numberOfGraphs>1){
        item *currentItem = perfect_matchings_counts;
        fprintf(stderr, "Size   Count\n");
        fprintf(stderr, "------------\n");
        while (currentItem != NULL) {
            fprintf(stderr, "%4d : %5d\n", currentItem->key, currentItem->value);
            totalPerfectMatchingsCount += (currentItem->value) * (currentItem->key);
            currentItem = currentItem->next;
        }
    } else {
        item *currentItem = perfect_matchings_counts;
        while (currentItem != NULL) {
            totalPerfectMatchingsCount += (currentItem->value) * (currentItem->key);
            currentItem = currentItem->next;
        }
    }
    fprintf(stderr, "\nQuadrangulations: %llu\n", numberOfGraphs);
    fprintf(stderr, "\nMatchings: %llu\n", totalPerfectMatchingsCount);
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
the right hand side of that edge. Faces are numbered 0,1,.... Also
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


int relabelling[MAXN];
int reverseRelabelling[MAXN];
EDGE *relabellingFirstedge[MAXN];
int relabellingQueue[MAXN];
int relabellingDegree[MAXN];

void constructCanonicallyLabeledGraph(EDGE *eStart){
    //calculate new labbeling
    int i;
    for(i=0; i<MAXN; i++){
        relabelling[i] = MAXN;
    }
    EDGE *e, *elast;
    int head = 1;
    int tail = 0;
    int vertexCounter = 1;
    relabellingQueue[0] = eStart->start;
    relabellingFirstedge[eStart->start] = eStart;
    relabelling[eStart->start] = 0;
    reverseRelabelling[0] = eStart->start;
    relabellingDegree[0] = degree[eStart->start];
    while(head>tail){
        int currentVertex = relabellingQueue[tail++];
        e = elast = relabellingFirstedge[currentVertex];
        do {
            if(relabelling[e->end]==MAXN){
                relabellingQueue[head++] = e->end;
                reverseRelabelling[vertexCounter] = e->end;
                relabellingDegree[vertexCounter] = degree[e->end];
                relabelling[e->end] = vertexCounter++;
                relabellingFirstedge[e->end] = e->inverse;
            }
            e = e->next;
        } while (e!=elast);
    }

    //perform relabelling
    for(i=0; i<nv; i++){
        EDGE *e, *elast;
        
        degree[i] = relabellingDegree[i];
    	
	firstedge[i] = relabellingFirstedge[reverseRelabelling[i]];
        e = elast = firstedge[i];
        do {
            e->start = relabelling[e->start];
            e->end = relabelling[e->end];
            e = e->next;
        } while (e!=elast);
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

    constructCanonicallyLabeledGraph(firstedge[0]);
    
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
//====================== USAGE =======================

void help(char *name){
    fprintf(stderr, "The program %s calculates perfect face matchings in quadrangulations.\n", name);
    fprintf(stderr, "This are perfect matchings in teh dual of a quadrangulation.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "Without any options, this program will calculate the number of perfect matchings\n");
    fprintf(stderr, "modulo the automoprhism group of the quadrangulation and report the number of perfect matchings.\n\n");
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices. Recompile if you need larger\n", MAXN);
    fprintf(stderr, "graphs.\n\n");
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, " -h, --help\n");
    fprintf(stderr, " Print this help and return.\n");
    fprintf(stderr, " -a, --all\n");
    fprintf(stderr, " Count all perfect matching (isomorphic copies will be counted)\n");
    fprintf(stderr, " -o, --output\n");
    fprintf(stderr, " Print the perfect matchings to stdout. (Currently not yet supported)\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options] n\n", name);
    fprintf(stderr, " %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]){
    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"all", no_argument, NULL, 'a'},
        {"output", no_argument, NULL, 'o'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hao", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'a':
                calculateAllPerfectMatchings = TRUE;
                break;
            case 'o':
                outputPerfectMatchings = TRUE;
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



    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        numberOfGraphs++;
        if(!calculateAllPerfectMatchings){
            calculateAutomorphismGroup();
        }
        generate_perfect_matchings_in_dual();
    }
    perfect_matchings_summary();

}
