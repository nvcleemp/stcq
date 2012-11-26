/* This program constructs the planar graph for gTRPZ_n.
*
*
* Compile with:
*
* cc -o gTRPZ -O4 gTRPZ.c
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

typedef int boolean;

#define TRUE 1
#define FALSE 0

void printAdjacencyList(int faces){
    int n = (faces-2)/4;
    int v = 1;
    int u = 2*n+3;

    int i, j;

    int currentVertex = 1;
    
    //v_0
    fprintf(stderr, "%2d) ", currentVertex++);
    fprintf(stderr, "%2d ", v+1);
    for (i = 0; i <= 2*n; i+=2){
        fprintf(stderr, "%2d ", u+i);
    }
    fprintf(stderr, "\n");
    
    //v_j (j=1,...,2n)
    for (j = 1; j <= 2*n; j++){
        fprintf(stderr, "%2d) ", currentVertex++);
        if(j%2){
            fprintf(stderr, "%2d ", v+j-1);
            fprintf(stderr, "%2d ", u+2*n+1);
            fprintf(stderr, "%2d ", v+j+1);
            fprintf(stderr, "\n");
        } else {
            fprintf(stderr, "%2d ", v+j+1);
            fprintf(stderr, "%2d ", u+0);
            fprintf(stderr, "%2d ", v+j-1);
            fprintf(stderr, "\n");
        }
    }
    
    //v_(2n+1)
    fprintf(stderr, "%2d) ", currentVertex++);
    fprintf(stderr, "%2d ", v+2*n);
    for (i = 2*n+1; i >= 1; i-=2){
        fprintf(stderr, "%2d ", u+i);
    }
    fprintf(stderr, "\n");
    
    //u_0
    fprintf(stderr, "%2d) ", currentVertex++);
    fprintf(stderr, "%2d ", u+1);
    for (i = 0; i <= 2*n; i+=2){
        fprintf(stderr, "%2d ", v+i);
    }
    fprintf(stderr, "\n");
    
    //u_j (j=1,...,2n)
    for (j = 1; j <= 2*n; j++){
        fprintf(stderr, "%2d) ", currentVertex++);
        if(j%2){
            fprintf(stderr, "%2d ", u+j-1);
            fprintf(stderr, "%2d ", v+2*n+1);
            fprintf(stderr, "%2d ", u+j+1);
            fprintf(stderr, "\n");
        } else {
            fprintf(stderr, "%2d ", u+j+1);
            fprintf(stderr, "%2d ", v+0);
            fprintf(stderr, "%2d ", u+j-1);
            fprintf(stderr, "\n");
        }
    }
    
    //u_(2n+1)
    fprintf(stderr, "%2d) ", currentVertex++);
    fprintf(stderr, "%2d ", u+2*n);
    for (i = 2*n+1; i >= 1; i-=2){
        fprintf(stderr, "%2d ", v+i);
    }
    fprintf(stderr, "\n");
}

/*
fills the array code with the planar code of the structure to which start belongs.
length will contain the length of the code. The maximum number of vertices is limited
to 255.
*/
void computePlanarCode(unsigned char code[], int *length, int faces) {
    int n = (faces-2)/4;
    int vertices = faces + 2;
    int edges = 2*faces;
    int v = 1;
    int u = 2*n+3;

    int i, j;
    unsigned char *codeStart;

    codeStart = code;
    //write order to code
    *code = (unsigned char) (vertices);
    code++;
    
    //v_0
    *code = (unsigned char) (v+1);
    code++;
    for (i = 0; i <= 2*n; i+=2){
        *code = (unsigned char) (u+i);
        code++;
    }
    *code = (unsigned char) (0);
    code++;
    
    //v_j (j=1,...,2n)
    for (j = 1; j <= 2*n; j++){
        if(j%2){
            *code = (unsigned char) (v+j-1);
            code++;
            *code = (unsigned char) (u+2*n+1);
            code++;
            *code = (unsigned char) (v+j+1);
            code++;
            *code = (unsigned char) (0);
            code++;
        } else {
            *code = (unsigned char) (v+j+1);
            code++;
            *code = (unsigned char) (u+0);
            code++;
            *code = (unsigned char) (v+j-1);
            code++;
            *code = (unsigned char) (0);
            code++;
        }
    }
    
    //v_(2n+1)
    *code = (unsigned char) (v+2*n);
    code++;
    for (i = 2*n+1; i >= 1; i-=2){
        *code = (unsigned char) (u+i);
        code++;
    }
    *code = (unsigned char) (0);
    code++;
    
    //u_0
    *code = (unsigned char) (u+1);
    code++;
    for (i = 0; i <= 2*n; i+=2){
        *code = (unsigned char) (v+i);
        code++;
    }
    *code = (unsigned char) (0);
    code++;
    
    //u_j (j=1,...,2n)
    for (j = 1; j <= 2*n; j++){
        if(j%2){
            *code = (unsigned char) (u+j-1);
            code++;
            *code = (unsigned char) (v+2*n+1);
            code++;
            *code = (unsigned char) (u+j+1);
            code++;
            *code = (unsigned char) (0);
            code++;
        } else {
            *code = (unsigned char) (u+j+1);
            code++;
            *code = (unsigned char) (v+0);
            code++;
            *code = (unsigned char) (u+j-1);
            code++;
            *code = (unsigned char) (0);
            code++;
        }
    }
    
    //u_(2n+1)
    *code = (unsigned char) (u+2*n);
    code++;
    for (i = 2*n+1; i >= 1; i-=2){
        *code = (unsigned char) (v+i);
        code++;
    }
    *code = (unsigned char) (0);
    code++;


    *length = code - codeStart;
    return;
}


/*
fills the array code with the planar code of the structure to which start belongs.
length will contain the length of the code. The maximum number of vertices is limited
to 65535.
*/
void computePlanarCodeShort(unsigned short code[], int *length, int faces) {
    int n = (faces-2)/4;
    int vertices = faces + 2;
    int edges = 2*faces;
    int v = 1;
    int u = 2*n+3;

    int i, j;
    unsigned short *codeStart;

    codeStart = code;
    //write order to code
    *code = (unsigned short) (vertices);
    code++;
    
    //v_0
    *code = (unsigned short) (v+1);
    code++;
    for (i = 0; i <= 2*n; i+=2){
        *code = (unsigned short) (u+i);
        code++;
    }
    *code = (unsigned short) (0);
    code++;
    
    //v_j (j=1,...,2n)
    for (j = 1; j <= 2*n; j++){
        if(j%2){
            *code = (unsigned short) (v+j-1);
            code++;
            *code = (unsigned short) (u+2*n+1);
            code++;
            *code = (unsigned short) (v+j+1);
            code++;
            *code = (unsigned short) (0);
            code++;
        } else {
            *code = (unsigned short) (v+j+1);
            code++;
            *code = (unsigned short) (u+0);
            code++;
            *code = (unsigned short) (v+j-1);
            code++;
            *code = (unsigned short) (0);
            code++;
        }
    }
    
    //v_(2n+1)
    *code = (unsigned short) (v+2*n);
    code++;
    for (i = 2*n+1; i >= 1; i-=2){
        *code = (unsigned short) (u+i);
        code++;
    }
    *code = (unsigned short) (0);
    code++;
    
    //u_0
    *code = (unsigned short) (u+1);
    code++;
    for (i = 0; i <= 2*n; i+=2){
        *code = (unsigned short) (v+i);
        code++;
    }
    *code = (unsigned short) (0);
    code++;
    
    //u_j (j=1,...,2n)
    for (j = 1; j <= 2*n; j++){
        if(j%2){
            *code = (unsigned short) (u+j-1);
            code++;
            *code = (unsigned short) (v+2*n+1);
            code++;
            *code = (unsigned short) (u+j+1);
            code++;
            *code = (unsigned short) (0);
            code++;
        } else {
            *code = (unsigned short) (u+j+1);
            code++;
            *code = (unsigned short) (v+0);
            code++;
            *code = (unsigned short) (u+j-1);
            code++;
            *code = (unsigned short) (0);
            code++;
        }
    }
    
    //u_(2n+1)
    *code = (unsigned short) (u+2*n);
    code++;
    for (i = 2*n+1; i >= 1; i-=2){
        *code = (unsigned short) (v+i);
        code++;
    }
    *code = (unsigned short) (0);
    code++;


    *length = code - codeStart;
    return;
}

/* void exportPlanarGraphCode(EDGE *start, int maxVertex) */

/*
exports the given structure as planarcode on stdout.
*/
void exportPlanarGraphCode(int faces, boolean header) {
    int vertices = faces + 2;
    int edges = 2*faces;

    int length;
    unsigned char code[2*edges+vertices+1];
    unsigned short codeShort[2*edges+vertices+1];
    static int first = TRUE;

    if (first && header) {
        fprintf(stdout, ">>planar_code<<");
        first = FALSE;
    }

    if (vertices + 1 <= 255) {
        computePlanarCode(code, &length, faces);
        if (fwrite(code, sizeof (unsigned char), length, stdout) != length) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    } else if (vertices + 1 <= 65535){
        computePlanarCodeShort(codeShort, &length, faces);
        putc(0, stdout);
        if (fwrite(codeShort, sizeof (unsigned short), length, stdout) != length) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    } else {
        fprintf(stderr, "Graph too large for planarcode -- exiting!\n");
        exit(-1);
    }
}
//====================== USAGE =======================

void help(char *name){
    fprintf(stderr, "The program %s constructs the planar graph for gTRPZ_f.\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options] f\n\n", name);
    fprintf(stderr, "Without any options, this program will construct the planar graph for gTRPZ_f.\n");
    fprintf(stderr, "This graph has f faces, and f is of the form 4n+2.\n\n");
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, " -h, --help\n");
    fprintf(stderr, " Print this help and return.\n");
    fprintf(stderr, " --human\n");
    fprintf(stderr, " Print human-readable adjacency list instead of planarcode.\n");
}

void usage(char *name){
    fprintf(stderr, "Usage: %s [options] f\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]){
    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
        {"human", no_argument, NULL, 0},
        {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;
    boolean human = FALSE;

    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                //handle long option with no alternative
                switch(option_index) {
                    case 0:
                        human = TRUE;
                        break;
                    default:
                        fprintf(stderr, "Illegal option index %d.\n", option_index);
                        usage(name);
                        return EXIT_FAILURE;
                }
                break;
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    // check the non-option arguments
    if (argc - optind != 1) {
        usage(name);
        return 1;
    }

    //parse the number of hexagons
    int faces = strtol(argv[optind], NULL, 10);

    if(faces%4 != 2){
        fprintf(stderr, "The number of faces is not of the form 4n+2.\n", name);
        usage(name);
        return 1;
    }

    if(!human){
        exportPlanarGraphCode(faces, TRUE);
    } else {
        printAdjacencyList(faces);
    }
}
