#include "adapt_parstate.h"

ParallelState::ParallelState(int N, MPI_Comm c)
{
    Nel  = N;
    comm = c;
     
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int nloc             = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(N/size) + MIN(rank, N%size);
     
    int* proc_nlocs                 = new int[size];
    int* proc_offset                = new int[size];
    nlocs                           = new int[size];
    offsets                         = new int[size];
    //std::cout << offset << std::endl;
    for(int i=0;i<size;i++)
    {
        nlocs[i]   = 0;
        offsets[i] = 0;
         
        if(i==rank)
        {
            proc_nlocs[i]  = nloc;
            proc_offset[i] = offset;
        }
        else
        {
            proc_nlocs[i]  = 0;
            proc_offset[i] = 0;
        }
    }
     
    MPI_Allreduce(proc_nlocs,  nlocs,   size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(proc_offset, offsets, size, MPI_INT, MPI_SUM, comm);
     
}// This is the constructor

ParallelState::~ParallelState()
{
    delete[] offsets;
    delete[] nlocs;
}

int* ParallelState::getOffsets( void )
{
    return offsets;
}

int* ParallelState::getNlocs( void )
{
    return nlocs;
}

int ParallelState::getOffset( int rank )
{
    return offsets[rank];
}

int ParallelState::getNloc( int rank )
{
    int nloc = nlocs[rank];
    return nloc;
}

int ParallelState::getNel( void )
{
  return Nel;
}
