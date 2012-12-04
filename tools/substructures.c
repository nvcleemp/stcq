/* This program reads quadrangulations from standard in and
 * looks for certain substructures.
 *
 *
 * Compile with:
 *
 * cc -o substructures -O4 substructures.c
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
unsigned long long int numberOfGraphsWithSubstructure = 0;

int nv; //the number of vertices of the current quadrangulation
int nf; //the number of faces of the current quadrangulation

//////////////////////////////////////////////////////////////////////////////

int cubeSearch(){
    int i;
    for(i=0; i<nf; i++){
        EDGE *e, *elast;
        
        int isCube = TRUE;
        
        e = elast = facestart[i];
        do {
            if(degree[e->start]!=3){
                isCube = FALSE;
            }
            e = e->inverse->prev;
        } while (e!=elast);
        
        if(isCube) return TRUE;
    }
    return FALSE;
}


int tristarSearch(){
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

            if(isTristar) return TRUE;
        }
    }
    return FALSE;
}


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
//====================== USAGE =======================

void help(char *name){
    fprintf(stderr, "The program %s looks for substructures in quadrangulations.\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices. Recompile if you need larger\n", MAXN);
    fprintf(stderr, "graphs.\n\n");
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, " -h, --help\n");
    fprintf(stderr, " Print this help and return.\n");
    fprintf(stderr, " -c, --cube\n");
    fprintf(stderr, " Looks for a quadrangle for which all vertices have degree 3.\n");
    fprintf(stderr, " -t, --tristar\n");
    fprintf(stderr, " Looks for a vertex of degree 3 for which all neighbours have degree 3.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]){
    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"cube", no_argument, NULL, 'c'},
        {"tristar", no_argument, NULL, 't'}
    };
    int option_index = 0;
    
    int (*searchFunction)() = NULL;

    while ((c = getopt_long(argc, argv, "hct", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case 'c':
                searchFunction = &cubeSearch;
                break;
            case 't':
                searchFunction = &tristarSearch;
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
    
    if(searchFunction==NULL){
        fprintf(stderr, "No substructure was specified.\n");
        usage(name);
        return EXIT_FAILURE;
    }



    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        numberOfGraphs++;
        if(searchFunction()){
            numberOfGraphsWithSubstructure++;
        }
    }
    
    fprintf(stderr, "%llu quadrangulations read.\n", numberOfGraphs);
    fprintf(stderr, "%llu quadrangulations have the requested substructure. (%3.4f%%)\n", numberOfGraphsWithSubstructure, (100.0*numberOfGraphsWithSubstructure)/numberOfGraphs);
    fprintf(stderr, "%llu quadrangulations do not have the requested substructure.  (%3.4f%%)\n", numberOfGraphs - numberOfGraphsWithSubstructure,  (100.0*(numberOfGraphs - numberOfGraphsWithSubstructure))/numberOfGraphs);
    

}
