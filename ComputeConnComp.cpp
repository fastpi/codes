#include "mex.h"
#include "matrix.h"
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <cstring>

using namespace std;

#include <stdio.h>
#include <signal.h>

short INTERRUPTED;
struct sigaction act, oact;

void sigproc_ctrlC(int sig)
{
	extern short INTERRUPTED;
	INTERRUPTED = 1;
}

//[plhs] = ComputeConnComp(A, k)
// Input
//    * A: an m x n rectangular matrix (i.e., m >> n)
//    * k: a hub selection ratio between 0 and 1
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   //
   extern short INTERRUPTED;
	INTERRUPTED = 0;
	act.sa_handler = sigproc_ctrlC;
	sigaction(SIGINT, &act, &oact);

   // input
   // * A_*: arrays for A (searching for neighbors of column nodes)
   // * AT_*: arrays for A' (searching for neighbors of row nodes)
   double*   A_pr = mxGetPr(prhs[0]);
   mwIndex*  A_ir = mxGetIr(prhs[0]);
   mwIndex*  A_jc = mxGetJc(prhs[0]);
   double*  AT_pr = mxGetPr(prhs[1]);
   mwIndex* AT_ir = mxGetIr(prhs[1]);
   mwIndex* AT_jc = mxGetJc(prhs[1]);
   //double       k = mxGetScalar(prhs[2]);
   double* row_topks = mxGetPr(prhs[2]);
   double* col_topks = mxGetPr(prhs[3]);
   mwSize m = mxGetM(prhs[0]); // # of rows
   mwSize n = mxGetN(prhs[0]); // # of columns
   mwSize nRowTopks = mxGetM(prhs[2]);
   mwSize nColTopks = mxGetM(prhs[3]);

   int N = m + n;
   // 0 ~ m - 1: for row nodes
   // m ~ N - 1: for column nodes

   //mexPrintf("[mex] nRows: %d\n", m);
   //mexPrintf("[mex] nCols: %d\n", n);
   //mexPrintf("[mex] nRowTopks: %d\n", nRowTopks);
   //mexPrintf("[mex] nColTopks: %d\n", nColTopks);

   // output
   mxArray* mx_cclabels = mxCreateDoubleMatrix(N, 1, mxREAL);
   double* cclabels = mxGetPr(mx_cclabels);
   memset(cclabels, 0.0, N * sizeof(double));

   //perform BFS
   int label = 1;

   for(int i = 0; i < nRowTopks; i++){
      mwIndex u = mwIndex(row_topks[i]) - 1;
      cclabels[u] = label++;
   }

   for(int i = 0; i < nColTopks; i++){
      mwIndex u = mwIndex(col_topks[i]) - 1 + m;
      cclabels[u] = label++;
   }

   for(int u = 0; u < N; u++){
      if(cclabels[u] == 0.0){
         // start BFS from u
         std::queue<int> Q;
         Q.push(u);
         cclabels[u] = label;

         while(!Q.empty()){
            int v = Q.front(); Q.pop();
            bool is_row_node = (v < m)   ? true  : false;
            mwIndex* ir = (is_row_node)  ? AT_ir :  A_ir;
            mwIndex* jc = (is_row_node)  ? AT_jc :  A_jc;
            int base    = (is_row_node)  ?     m :     0;
            v           = (is_row_node)  ?     v : v - m;

            for(int i = jc[v]; i < jc[v+1]; i++){
               int w = ir[i] + base;
               if(cclabels[w] == 0.0){
                  Q.push(w);
                  cclabels[w] = label;
               }//if
            }//for

            if(INTERRUPTED) break;
         }//while

         label++;
      }//if

      if(INTERRUPTED == 1) {
			mexPrintf("Ctrl-C interrupted...\n");
			break;
		}
   }//for

   plhs[0] = mx_cclabels;

}
