#include "adapt_operations.h"


using namespace std;




int FindRank(int* arr, int size, int val)
{
    int start = 0;
    int last  = size-1;
    
    int mid   = (start+last)/2;
    
    while (start<=last)
    {
        if (arr[mid]<val)
        {
            start = mid + 1;
        }
        else
        {
            last  = mid - 1;
        }
        mid = (start+last)/2;
    }

    return mid;
}



int FindBoundaryID(int* arr, int size, int val)
{
    int start = 0;
    int last  = size-1;
    
    int mid   = (start+last)/2;
    
    while (start<=last)
    {
        if (arr[mid]<=val)
        {
            start = mid + 1;
        }
        else
        {
            last  = mid - 1;
        }
        mid = (start+last)/2;
    }

    return mid;
}



std::vector<int> FindDuplicates(std::vector<int> arr)
{
    int N = arr.size();
    sort(arr.begin(),arr.end());
    std::vector<int> res;
    std::set<int> check;
    for(int i=0;i<N;i++)
    {
        if(arr[i+1]==arr[i])
        {
            if(check.find(arr[i])==check.end())
            {
                check.insert(arr[i]);
                res.push_back(arr[i]);
            }
        }
    }
    
    return res;
}


int binarySearch(int* arr, int low, int high, int key)
{
    if (high < low)
        return -1;
    int mid = (low + high) / 2; /*low + (high - low)/2;*/
    if (key == arr[mid])
        return mid;
    if (key > arr[mid])
        return binarySearch(arr, (mid + 1), high, key);
    
    return binarySearch(arr, low, (mid - 1), key);
}



std::vector<int> FindDuplicatesInParallel_VecV2(std::vector<int> arr, int arr_size, int glob_size, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int levels = log2(size);
    
    int N = arr_size;

    //if (rank == 0)
    //{
    //glob_arr = new int[glob_size];
    //}
    std::vector<int> glob_arr(glob_size);
    std::vector<int> sorted = mergeSort_vec(levels, rank, arr, N, comm, glob_arr);
    
    MPI_Bcast(&sorted[0], glob_size, MPI_INT, 0, MPI_COMM_WORLD);

    ParallelState* pv = new ParallelState(glob_size, comm);
    
    std::vector<int> res;
    std::set<int> check;

    for(int i=0;i<pv->getNloc(rank);i++)
    {
        if(sorted[pv->getOffset(rank)+i+1]==sorted[pv->getOffset(rank)+i])
        {
            check.insert(pv->getOffset(rank)+sorted[i]);
            res.push_back(pv->getOffset(rank)+sorted[i]);
        }
    }
  
    
    int* dupl_locs     = new int[size];
    int* red_dupl_locs = new int[size];

    
    for(int i=0;i<size;i++)
    {
        red_dupl_locs[i]  = 0;
        
        if(i==rank)
        {
            dupl_locs[i]  = res.size();
        }
        else
        {
            dupl_locs[i]  = 0;
        }
    }
    
    MPI_Allreduce(dupl_locs, red_dupl_locs, size, MPI_INT, MPI_SUM, comm);
    
    int* red_dupl_offsets = new int[size];
    red_dupl_offsets[0] = 0;
    for(int i=0;i<size-1;i++)
    {
        red_dupl_offsets[i+1]=red_dupl_offsets[i]+red_dupl_locs[i];
    }
    
    int tot_dupl = red_dupl_offsets[size-1]+red_dupl_locs[size-1];
    
    std::vector<int> duplicates(tot_dupl);
    
    MPI_Allgatherv(&res[0],
                   res.size(),
                   MPI_INT,
                   &duplicates[0],
                   red_dupl_locs,
                   red_dupl_offsets,
                   MPI_INT, comm);
    
    delete pv;
    return duplicates;
}


std::vector<int> FindDuplicatesInParallel(int* arr, int arr_size, int glob_size, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int levels = log2(size);
    
    int N = arr_size;

    int* glob_arr = new int[glob_size];
    
    int* sorted = mergeSort(levels, rank, arr, N, comm, glob_arr);
    
    MPI_Bcast(sorted, glob_size, MPI_INT, 0, MPI_COMM_WORLD);

    ParallelState* pv = new ParallelState(glob_size, comm);
    
    std::vector<int> res;
    std::set<int> check;

    for(int i=0;i<pv->getNloc(rank);i++)
    {
        if(sorted[pv->getOffset(rank)+i+1]==sorted[pv->getOffset(rank)+i])
        {
            check.insert(sorted[pv->getOffset(rank)+i]);
            res.push_back(sorted[pv->getOffset(rank)+i]);
        }
    }
    
    
    int* dupl_locs     = new int[size];
    int* red_dupl_locs = new int[size];

    
    for(int i=0;i<size;i++)
    {
        red_dupl_locs[i]  = 0;
        
        if(i==rank)
        {
            dupl_locs[i]  = res.size();
        }
        else
        {
            dupl_locs[i]  = 0;
        }
    }
    
    MPI_Allreduce(dupl_locs,  red_dupl_locs,  size, MPI_INT, MPI_SUM, comm);
    
    int* red_dupl_offsets = new int[size];
    red_dupl_offsets[0] = 0;
    for(int i=0;i<size-1;i++)
    {
        red_dupl_offsets[i+1]=red_dupl_offsets[i]+red_dupl_locs[i];
    }
    
    int tot_dupl = red_dupl_offsets[size-1]+red_dupl_locs[size-1];
    
    std::vector<int> duplicates(tot_dupl);
    MPI_Allgatherv(&res[0],
                   res.size(),
                   MPI_INT,
                   &duplicates[0],
                   red_dupl_locs,
                   red_dupl_offsets,
                   MPI_INT, comm);
    
    
    delete[] glob_arr;
    delete[] dupl_locs;
    delete[] red_dupl_locs;
    delete[] red_dupl_offsets;
    delete pv;
    return duplicates;
}


int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int* merge(int* a, int* b, int size_a, int size_b, int* merged, int size)
{
    
    int i=0;
    int j=0;
    int k=0;
    
    while(i<=size_a-1 && j<=size_b-1)
    {
        if(a[i] <= b[j])
        {
            merged[k++] = a[i++];
        }
        else
        {
            merged[k++] = b[j++];
        }
    }
    while(i <= size_a-1)
    {
        merged[k++] = a[i++];
    }
    while(j <= size_b-1)
    {
        merged[k++] = b[j++];
    }
    
    return merged;
}


std::vector<int> merge_vec(std::vector<int> a, std::vector<int> b)
{
    int n = a.size();
    int m = b.size();
    int size = n+m;
    
    std::vector<int> merged(size);
    int i=0;
    int j=0;
    int k=0;
    
    while(i<=n-1 && j<=m-1)
    {
        if(a[i] <= b[j])
        {
            merged[k++] = a[i++];
        }
        else
        {
            merged[k++] = b[j++];
        }
    }
    while(i <= n-1)
    {
        merged[k++] = a[i++];
    }
    while(j <= m-1)
    {
        merged[k++] = b[j++];
    }
    
    return merged;
}


int* mergeSort(int height, int id, int* localArray, int size, MPI_Comm comm, int* globalArray)
{
    int parent, rightChild, myHeight;
    int *half1, *half2, *mergeResult;

    myHeight = 0;
    qsort(localArray, size, sizeof(int), compare); // sort local array
    half1 = localArray;  // assign half1 to localArray
    int size_half1, size_half2;
    
    while (myHeight < height) { // not yet at top
        parent = (id & (~(1 << myHeight)));

        if (parent == id)
        { // left child
              rightChild = (id | (1 << myHeight));

              // allocate memory and receive array of right child
              
              MPI_Recv(&size_half2, 1, MPI_INT, rightChild, 1234, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
              half2 = new int[size_half2];
            
              MPI_Recv(half2, size_half2, MPI_INT, rightChild, 5678, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
              // allocate memory for result of merge
              
              // merge half1 and half2 into mergeResult
              size_half1 = size;
              
              //mergeResult = (int*) malloc (size_half2+size_half1);
              mergeResult = new int[size_half2+size_half1];
              mergeResult = merge(half1, half2, size_half1, size_half2, mergeResult, size);
              // reassign half1 to merge result
            half1 = mergeResult;
            size = size_half1+size_half2;  // double size
            
            delete[] half2;
            mergeResult = NULL;

            myHeight++;

        }
        else
        {
            // right child
            // send local array to parent
            MPI_Send(&size,    1, MPI_INT, parent, 1234, MPI_COMM_WORLD);
            MPI_Send(half1, size, MPI_INT, parent, 5678, MPI_COMM_WORLD);
            if(myHeight != 0)
            {
                //delete[] half1;
            }
            myHeight = height;
        }
    }

    if(id == 0){
        globalArray = half1;   // reassign globalArray to half1
    }
    
    return globalArray;
}

//std::vector<int> MergeAndSort(std::vector<int> localArray, std::vector<int> globalArray, MPI_Comm comm)
//{
//
//}

std::vector<int> mergeSort_vec(int height, int rank, std::vector<int> localArray, int size, MPI_Comm comm, std::vector<int> globalArray){
    
    int parent, rightChild, myHeight;

    myHeight = 0;
    sort(localArray.begin(),localArray.end());
    //qsort(localArray, size, sizeof(int), compare); // sort local array
    std::vector<int> half1 = localArray;  // assign half1 to localArray
    int size_half1, size_half2;
    while (myHeight < height) { // not yet at top
        parent = (rank & (~(1 << myHeight)));

        if (parent == rank) { // left child
            rightChild = (rank | (1 << myHeight));

              // allocate memory and receive array of right child
              //half2 = (int*) malloc (size * sizeof(int));
            MPI_Recv(&size_half2, 1, MPI_INT, rightChild, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::vector<int> half2(size_half2);
            MPI_Recv(&half2[0], size_half2, MPI_INT, rightChild, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            size_half1 = size;
            // allocate memory for result of merge
            //mergeResult = (int*) malloc (size * 2 * sizeof(int));
              
            // merge half1 and half2 into mergeResult
            std::vector<int> mergeResult = merge_vec(half1, half2);
            // reassign half1 to merge result
            half1 = mergeResult;
            
            
            size = size_half1+size_half2;  // the size of mergeResult is the size of half1 and half2 added togehter.
            half2.erase(half2.begin(),half2.end());
            //free(half2);
            mergeResult.erase(mergeResult.begin(),mergeResult.end());

            myHeight++;

        } else { // right child
              // send local array to parent
            MPI_Send(&size,    1, MPI_INT, parent, 0, MPI_COMM_WORLD);

            MPI_Send(&half1[0], size, MPI_INT, parent, 0, MPI_COMM_WORLD);
            if(myHeight != 0)
            {
                //half1.erase(half1.begin(),half1.end());
            }
            myHeight = height;
        }
    }

    if(rank == 0){
        globalArray = half1;   // reassign globalArray to half1
    }
    return globalArray;
}


int largest(int arr[], int n)
{
    int i;
      
    // Initialize maximum element
    int max = arr[0];
  
    // Traverse array elements
    // from second and compare
    // every element with current max
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
  
    return max;
}




void TestFindRank(MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    int* arr = new int[10];
    
    arr[0] = 0;
    arr[1] = 10;
    arr[2] = 12;
    arr[3] = 24;
    arr[4] = 42;
    arr[5] = 55;
    arr[6] = 65;
    arr[7] = 66;
    arr[8] = 78;
    arr[9] = 81;
    
    int value = 82;
    
    int res = FindRank(arr,10,value);
    
    if(res == 3)
    {
        std::cout << "TestFindRank() has passed! " << res << std::endl;
    }
    else
    {
        
        std::cout << rank << "  TestFindRank() has failed! " << res << std::endl;
        std::cout << std::endl;
    }
    
}













//void ParallelSortTest_v2()
//{
//    MPI_Comm comm = MPI_COMM_WORLD;
//    MPI_Info info = MPI_INFO_NULL;
//    int size;
//    MPI_Comm_size(comm, &size);
//    // Get the rank of the process
//    int rank;
//    MPI_Comm_rank(comm, &rank);
//    
//    int levels = log2(size);
//    
//    const char* fn_conn="grids/piston/conn.h5";
//    //const char* fn_grid="grids/piston/grid.h5";
//    
//    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
//    
//    int nrow = ief->getNrow();
//    int ncol = ief->getNcol();
//    
//    int *data = new int[nrow*(ncol-1)];
//    std::vector<int> ief_copy(nrow*(ncol-1));
//    int fid;
//    int cnt=0;
//    for(int i=0;i<nrow;i++)
//    {
//        for(int j=0;j<ncol-1;j++)
//        {
//            fid = fabs(ief->getVal(i,j+1))-1;
//            ief_copy[cnt] = fid;
//            cnt++;
//        }
//    }
//    int lsize = ief_copy.size();
//    int* ief_arr = new int[lsize];
//    for(int i=0;i<lsize;i++)
//    {
//        ief_arr[i] = ief_copy[i];
//    }
//    //int *glob_arr;
////    if (rank == 0)
////    {
////        glob_arr = new int[ief->nglob*6];
////
////    }
//    int nglob = ief->getNglob();
//    std::vector<int> glob_arr(nglob*6);
//    
//    std::vector<int> sorted = mergeSort_vec(levels, rank, ief_copy, lsize, comm, glob_arr);
//    
//    if(rank == 0)
//    {
//        std::cout << "hoi " << rank << " " << nglob*6 << std::endl;
//       for(int i=0;i<nglob*6;i++)
//        {
//           std::cout << "vec " << rank << " " << i << " " << sorted[i] << std::endl;
//        }
//    }
//    
//    delete[] ief_arr;
//    delete[] data;
//    delete ief;
//    
//    ief_copy.erase(ief_copy.begin(),ief_copy.end());
//    
//}




//void MergeTest()
//{
//    std::vector<int> vec;
//    vec.push_back(1);
//    vec.push_back(3);
//    vec.push_back(2);
//    vec.push_back(14);
//    vec.push_back(8);
//    vec.push_back(2);
//    vec.push_back(3);
//    vec.push_back(14);
//    vec.push_back(9);
//    vec.push_back(15);
//    vec.push_back(1);
//    vec.push_back(14);
//    sort(vec.begin(),vec.end());
//    std::cout << "==vec1==" << std::endl;
//    for(int i = 0;i<vec.size();i++)
//    {
//        std::cout << vec[i] << std::endl;
//    }
//    std::cout << "====" << std::endl;
//    std::vector<int> vec2;
//    vec2.push_back(2);
//    vec2.push_back(6);
//    vec2.push_back(5);
//    vec2.push_back(7);
//    vec2.push_back(9);
//    vec2.push_back(4);
//    vec2.push_back(31);
//    vec2.push_back(4);
//    sort(vec2.begin(),vec2.end());
//    std::cout << "==vec2==" << std::endl;
//    for(int i = 0;i<vec2.size();i++)
//    {
//        std::cout << vec2[i] << std::endl;
//    }
//    std::cout << "====" << std::endl;
//    
//    std::vector<int> res = merge_vec(vec,vec2);
//    
//    std::cout << "====" << std::endl;
//    for(int i = 0;i<res.size();i++)
//    {
//        std::cout << res[i] << std::endl;
//    }
//    std::cout << "====" << std::endl;
//}
//
//
//
//void FindDuplicateTest()
//{
//    std::vector<int> vec;
//    vec.push_back(1);
//    vec.push_back(3);
//    vec.push_back(2);
//    vec.push_back(14);
//    vec.push_back(8);
//    vec.push_back(2);
//    vec.push_back(3);
//    vec.push_back(14);
//    vec.push_back(9);
//    vec.push_back(15);
//    vec.push_back(1);
//    vec.push_back(14);
//    std::vector<int> res = FindDuplicates(vec);
//    
//    std::cout << "====" << std::endl;
//    for(int i = 0;i<res.size();i++)
//    {
//        std::cout << res[i] << std::endl;
//    }
//    std::cout << "====" << std::endl;
//}
//
//
//
//
//
//void ParallelSortTest()
//{
//    int i=0;
//    MPI_Comm comm = MPI_COMM_WORLD;
//    MPI_Info info = MPI_INFO_NULL;
//    int size;
//    MPI_Comm_size(comm, &size);
//    // Get the rank of the process
//    int rank;
//    MPI_Comm_rank(comm, &rank);
//    
//    int levels = log2(size);
//    
//    const char* fn_conn="grids/piston/conn.h5";
//    
//    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
//    
//    int nrow = ief->getNrow();
//    int ncol = ief->getNcol();
//    
//    std::vector<int> ief_copy(nrow*(ncol-1));
//    int fid;
//    int cnt=0;
//    for(i=0;i<nrow;i++)
//    {
//        for(int j=0;j<ncol-1;j++)
//        {
//            fid = fabs(ief->getVal(i,j+1))-1;
//            //std::cout << cnt << " " << fid << std::endl;
//            ief_copy[cnt] = fid;
//            cnt++;
//        }
//    }
//    
//    int lsize = ief_copy.size();
//    int* ief_arr = new int[lsize];
//    
//    for(i=0;i<lsize;i++)
//    {
//        ief_arr[i] = ief_copy[i];
//    }
//    
//    int *glob_arr = NULL;
//    int nglob = ief->getNglob();
//
//    if (rank == 0)
//    {
//        glob_arr = new int[nglob*6];
//    }
//    int* sorted = mergeSort(levels, rank, ief_arr, lsize, comm, glob_arr);
//    
//    if(rank == 0)
//    {
//       for(int i=0;i<nglob*6;i++)
//        {
//           std::cout << rank << " " << i << " " << sorted[i] << std::endl;
//        }
//    }
//    
//    
//    delete ief;
//    ief_copy.erase(ief_copy.begin(),ief_copy.end());
//    
//}



