#include "solver.h"

/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/

// Helper functions will be explained in detail further on. I am the only team member on this,
// and I recieved whiteboard-level MPI help from a.) ak07@cornell.edu b.) Peter Ludwig at URome.
	
void seq_solver_recursive(std::vector<std::vector<unsigned int> >& solns,
					      std::vector<unsigned int> &psoln,
		                  std::vector<bool> &CQ,
						  std::vector<bool> &XQ1,
						  std::vector<bool> &XQ2,
		                  unsigned int r, 
		                  unsigned int& n,
		                  bool EoF
		                );

bool seq_solver_rp(MPI_Request& request, 
			 MPI_Status& status,
			 std::vector<std::vector<unsigned int> >& CB,
			 std::vector<std::vector<unsigned int> >& psolns, 
			 std::vector<unsigned int>& ssolns,
			 int c, 
			 int n, 
			 int k, 
			 int& SC, 
			 int& RC, 
			 int& F, 
			 int& SS 
			);

/*************************** solver.h functions ************************/

// Here is the sequential solution to the problem. Here I have stored the values for if a queen is present
// in 3 vectors corresponding to column, and the two diagonals. The idea is to search along them for any
// conflicting queens, and then if a solution is found, to store it. Then we recursively backtrack to a 
// previous point in our search, and check a different solution. EoF is implemented by stopping after one
// solution is detected in the answer vector.

void seq_solver(unsigned int n, unsigned int exit_on_first, std::vector<std::vector<unsigned int> >& solns) {

	std::vector<unsigned int> psoln;
	std::vector<bool> CQ (n , 0);
	std::vector<bool> XQ1 (2*n - 1 , 0);
	std::vector<bool> XQ2 (2*n - 1 , 0);
	
	if (exit_on_first == false) seq_solver_recursive(solns, psoln, CQ, XQ1, XQ2, 0, n, false);
	else seq_solver_recursive(solns, psoln, CQ, XQ1, XQ2, 0, n, true);
	
	return;	

}

void seq_solver_recursive(std::vector<std::vector<unsigned int> >& solns,
					      std::vector<unsigned int> &psoln,
		                  std::vector<bool> &CQ,
						  std::vector<bool> &XQ1,
						  std::vector<bool> &XQ2,
		                  unsigned int r, 
		                  unsigned int& n,
		                  bool EoF
		                ) {

	if(r == n) {
		std::vector<unsigned int> temp;
		temp = psoln;
		solns.push_back(temp);
		return;
	}
	
	for(unsigned int c = 0; c < n; c++) {
		if(CQ[c] == false && XQ1[r + c] == false && XQ2[n - 1 + c - r] == false) {
			CQ[c] = true;
			XQ1[r + c] = true;
			XQ2[n - 1 + c - r] = true;
			psoln.push_back(c);
			seq_solver_recursive(solns, psoln, CQ, XQ1, XQ2, r + 1, n, EoF);
			if (EoF == true && solns.size()!= 0) return;
			psoln.pop_back();
			CQ[c] = false;
			XQ1[r + c] = false;
			XQ2[n - 1 + c - r] = false;
		}
	}
}

// For the parallel part I chose a different approach, storing all queen positions on chess board CB.
// My general approach is verbatim to the approach shown in the TAs comments.

bool seq_solver_rp(MPI_Request& request, 
			 MPI_Status& status,
			 std::vector<std::vector<unsigned int> >& CB,
			 std::vector<std::vector<unsigned int> >& psolns, 
			 std::vector<unsigned int>& ssolns,
			 int c, 
			 int n, 
			 int k, 
			 int& SC, 
			 int& RC, 
			 int& F, 
			 int& SS
			) {

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    static int r = 1;   
    
// Master continuously computes partial solutions and computes them regardless of workers.
	
    if (rank == 0 && size > 1 && c == k) {
			
            static int holder = 0;
			std::vector<unsigned int> tfiller(n);
			psolns.push_back(tfiller);
			for (int i = 0; i < n; i++) {
				psolns[holder].resize(n);
				for (int j = 0; j < n; j++) if (CB[j][i] == 1) psolns[holder][i] = j;
			}
			holder++;

            if (r < size) {
                MPI_Send(&psolns[SC][0], n, MPI_UNSIGNED, r, 111, MPI_COMM_WORLD);
                SC++;       
            }
			
            else {
                if (F != 0) {
                    MPI_Irecv(&SS, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);
                    F = 0;
                }
                MPI_Test(&request, &F, &status);
                if (F != 0) {
                    RC++;
                    if (SS != 0) {
                        ssolns.resize(ssolns.size() + SS);
                        MPI_Recv(&ssolns[ssolns.size() - SS], SS, MPI_UNSIGNED, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
                    }
                    MPI_Send(&psolns[SC][0], n, MPI_UNSIGNED, status.MPI_SOURCE, 111, MPI_COMM_WORLD);
                    SC++;  
                    F = -1;
                }
            }
			
            r++;      
            return true;

    }
	
// This does exactly the same backtracking logic as the sequential part.
	
    else if (c == n) {
            psolns.resize(0);
            psolns.push_back(std::vector<unsigned int>(n));
			for (int i = 0; i < n; i++) {
				psolns[0].resize(n);
				for (int j = 0; j < n; j++) if (CB[j][i] == 1) psolns[0][i] = j;
			}
            for (int j = 0; j < n; j++) ssolns.push_back(psolns[0][j]);
            return true;
    }
    
    bool pcheck = false;
    for (int p = 0; p < n; p++)
    {
		bool check = true;
		int i, j;
		for (i = 0; i < c; i++) if (CB[p][i]) check = false;
		for (i = p, j = c; i >= 0 && j >= 0; i--, j--) if (CB[i][j]) check = false;
		for (i = p, j = c; j >= 0 && i < n; i++, j--) if (CB[i][j]) check = false;
        if (check)
        {
            CB[p][c] = 1;
            pcheck = seq_solver_rp(request, status, CB, psolns, ssolns, c + 1, n, k, SC, RC, F, SS) || pcheck;
            CB[p][c] = 0;
        }
    }
    return pcheck;
}

// Master starts by first computing the first k partial sums.

void nqueen_master( unsigned int n,
                    unsigned int k,
                    unsigned int exit_on_first,
                    std::vector<std::vector<unsigned int>>& solns) {

    int size, SC, RC, ctr, SS, psolnsize;                         
    int F = -1;                                         
    
	std::vector<std::vector<unsigned int>> psolns;     
    std::vector<std::vector<unsigned int>> CB;          
    std::vector<unsigned int> ssolns;
	
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_Request request;
	
	if (exit_on_first == true) {
		seq_solver(n, exit_on_first, solns);
		for (int i = 1; i < size; i++) MPI_Send(0, 0, MPI_INT, i, 222, MPI_COMM_WORLD);
		return;
	}

    CB.resize(n);
    for (unsigned int i = 0; i < n; i++) CB[i].resize(n, 0);
    
    SC = 0;         
    RC = 0;
    
    seq_solver_rp(request, status, CB, psolns, ssolns, 0, n, k, SC, RC, F, SS);
    psolnsize = psolns.size();
	
// This loop sends all the workers to their work.
	
    while (RC < (psolnsize - size + 1)) {
        if (F != 0) {
            MPI_Irecv(&SS, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);
            F = 0;
        }
        
        MPI_Test(&request, &F, &status);
		
        if (F != 0) {
            RC++;
			
            if (SS != 0) {
                ssolns.resize(ssolns.size() + SS);
                MPI_Recv(&ssolns[ssolns.size() - SS], SS, MPI_UNSIGNED, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
            }
            
            MPI_Send(&psolns[SC][0], n, MPI_UNSIGNED, status.MPI_SOURCE, 111, MPI_COMM_WORLD);
            SC++; 
            F = -1;
        }
        
    }

// This loop in essence recieves all partial sums and then kills workers.
  
    (SC < size)?(ctr = SC + 1):(ctr = size);
    
    for(int FRec = 1; FRec < ctr; ) {
        if (F != 0) {
            MPI_Irecv(&SS, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &request);
            F = 0;
        }
		
        MPI_Test(&request, &F, &status);
		
        if (F != 0) {
            FRec++;
            if (SS != 0) {
                ssolns.resize(ssolns.size() + SS);
                MPI_Recv(&ssolns[ssolns.size() - SS], SS, MPI_UNSIGNED, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
                F = -1;
            }
        }
    }
	
    for (int i = 1; i < size; i++) MPI_Send(0, 0, MPI_INT, i, 222, MPI_COMM_WORLD);
    
//  This compiles all the answers into the proper form.
	
    solns.resize(ssolns.size() / n);       
    int place = 0;
    for (unsigned int i = 0; i < ssolns.size() / n; i++) {
        solns[i].resize(n);
        for (unsigned int j = 0; j < n; j++) {
            solns[i][j] = ssolns[place];
            place++;
        }
    }

}

void nqueen_worker( unsigned int n,
                    unsigned int k,
                    unsigned int exit_on_first) {

    int rank, SS, SC, RC;        
    int F = -1;
	
	std::vector<std::vector<unsigned int>> psolns;     
    std::vector<unsigned int> ssolns;                    
    std::vector<unsigned int> sinsolns;                 
    std::vector<std::vector<unsigned int>> CB;
	
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    MPI_Request request;
    
	CB.resize(n);
    for (unsigned int i = 0; i < n; i++) CB[i].resize(n, 0);
    sinsolns.resize(n);
    
// While loops rotate the work around until  workers are killed by the master.
	
    while(true) {
        ssolns.clear();
        sinsolns.clear();
        sinsolns.resize(n);
        CB.clear();
        CB.resize(n);
		
        for (unsigned int i = 0; i < n; i++) CB[i].resize(n, 0);
        MPI_Recv(&sinsolns[0], n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        
		if (status.MPI_TAG == 222) return;     
        
		for (unsigned int i = 0; i < k; i++) CB[sinsolns[i]][i] = 1;
        
		bool SF = true;
        
		if (!seq_solver_rp(request, status, CB, psolns, ssolns, k, n, k, SC, RC, F, SS)) {
            SS = 0;
            MPI_Send(&SS, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            SF = false ;
        }
        
		if (SF) {
            SS = ssolns.size();
            MPI_Send(&SS, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&ssolns[0], SS, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
        }
    }
	
    return;
	
}