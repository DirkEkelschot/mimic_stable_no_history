#include "adapt_geometry.h"
#include "adapt_partition.h"
#include "adapt_geometry.h"
#include "adapt_compute.h"

#ifndef ADAPT_TOPOLOGY_H
#define ADAPT_TOPOLOGY_H

using namespace std;

struct Mesh_Topology_BL{
    int Nprisms;
    std::map<int,std::vector<int> > BLlayers;
    std::map<int,std::vector<std::vector<int> > > BLlayersPrisms;
    std::map<int,std::vector<std::vector<int> > > bcQuad;
    std::map<int,std::vector<std::vector<int> > > bcTria;
    std::map<int,std::vector<Element* > > BLlayersElements;
    std::vector<std::vector<int> > BndFaces;

};



class Mesh_Topology {
    public:
        Mesh_Topology(){};
        Mesh_Topology(Partition* Pa, MPI_Comm comm);
        ~Mesh_Topology();
        void DetermineBoundaryLayerElements(Partition* Pa, Array<int>* ife_in, int nLayer, int bID, MPI_Comm comm);
        std::map<int,std::vector<int> > getScheme_E2V();
        std::map<int,vector<Vec3D*> > getNormals();
        std::map<int,vector<Vec3D*> > getRvectors();
        std::map<int,vector<Vec3D*> > getdXfXc();
        std::map<int,vector<double> > getdr();
        std::map<int,vector<double> > getdS();
        Array<int>* getIFN();
        std::map<int,double> getVol();
        std::map<int,std::map<int,double> > GetElement2VertexScheme();
        //std::map<int,double> ReduceFieldToVertices(Domain* pDom, std::map<int,double> Uelem);
        //std::map<int,Array<double>* > ReduceMetricToVertices(Domain* pDom, std::map<int,Array<double>* > Telem);
        std::map<int,int> getFace2Ref();
        std::map<int,std::vector<int> > getRef2Face();
        std::map<int,int> getVert2Ref();
        std::map<int,std::vector<int> > getRef2Vert();
        std::map<int,std::vector<Vert*> > getVfacevector();
        Mesh_Topology_BL* getBLMeshTopology();
        std::map<int,Vert*> getGhostVerts();
    
    private:
        std::map<int,std::vector<Vert*> > vfacevector;
        std::map<int,std::vector<Vec3D*> > normals;
        std::map<int,std::vector<Vec3D*> > rvector;
        std::map<int,std::vector<Vec3D*> > dxfxc;
        std::map<int,std::vector<double> > dr;
        std::map<int,std::vector<double> > dS;
        std::map<int,double> Vol;
        Array<int>* ifn;
        std::map<int,int> face2ref;
        std::map<int,std::vector<int> > ref2face;
        std::map<int,int> vert2ref;
        std::map<int,std::vector<int> > ref2vert;
        std::map<int,int> Bface2Element;
        std::map<int,int> Bface2LocID;
        std::map<int,Vec3D*> Bface2Normal;
        std::map<int,std::vector<int> > BLlayers; 
        Mesh_Topology_BL* mesh_topo_bl;
        std::map<int,std::vector<int> > E2V_scheme;
        std::map<int,Vert*> ghostVerts;
    
        MPI_Comm c;
};


#endif
