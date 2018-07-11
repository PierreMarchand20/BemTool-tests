#ifndef BEMTOOLTESTS_MISC_RENUM_HPP
#define BEMTOOLTESTS_MISC_RENUM_HPP

#include <string>
#include <fstream>
#include <bemtool/tools.hpp>

namespace bemtool{

template<int dim>
std::vector<std::pair<int,int>> boundary_edges(Dof< BasisFct<P1,dim>> dof){

    int nb_dofs = NbDof(dof);
    int nb_elts = NbElt(dof);

    std::vector<int> first(nb_dofs,-1);
    std::vector<int> next;
    std::vector<std::pair<int,int>> edges;
    std::vector<int> is_boundary;
    int NbNotBoundary=0;

    for (int i=0;i<nb_elts;i++){
        const array<dim+1,int>& jdof = dof[i];
        for (int j=0;j<dim+1;j++){
            int v0 = jdof[j%3];
            int v1 = jdof[(j+1)%3];
            if (v1<v0){
                int c =v0;
                v0=v1;
                v1=c;
            }

            std::pair<int,int> new_edge(v0,v1);
            int I= first[v0];
            bool exist=0;
            while (I!=-1){
                if (edges[I]==new_edge){
                    is_boundary[I]=0;
                    NbNotBoundary+=1;
                    exist=1;
                    break;
                }
                else{
                    I=next[I];
                }
            }
            if(exist==0){
                next.push_back(first[v0]);
                first[v0]= next.size()-1;
                edges.push_back(new_edge);
                is_boundary.push_back(1);
            }
        }
    }

    std::vector<std::pair<int,int>> boundary(edges.size()-NbNotBoundary);
    int count=0;
    for (int i=0;i<edges.size();i++){
        if(is_boundary[i]){
            boundary[count]=edges[i];
            count++;
        }
    }

    return boundary;
}

template<int dim>
std::vector<int> renum_wo_boundary(Dof< BasisFct<P1,dim>> dof){

    int nb_dofs = NbDof(dof);
    int nb_elts = NbElt(dof);

    std::vector<int> first(nb_dofs,-1);
    std::vector<int> next;
    std::vector<std::pair<int,int>> edges;
    std::vector<int> is_boundary;
    int NbNotBoundary=0;

    for (int i=0;i<nb_elts;i++){
        const array<dim+1,int>& jdof = dof[i];
        for (int j=0;j<dim+1;j++){
            int v0 = jdof[j%3];
            int v1 = jdof[(j+1)%3];
            if (v1<v0){
                int c =v0;
                v0=v1;
                v1=c;
            }

            std::pair<int,int> new_edge(v0,v1);
            int I= first[v0];
            bool exist=0;
            while (I!=-1){
                if (edges[I]==new_edge){
                    is_boundary[I]=0;
                    NbNotBoundary+=1;
                    exist=1;
                    break;
                }
                else{
                    I=next[I];
                }
            }
            if(exist==0){
                next.push_back(first[v0]);
                first[v0]= next.size()-1;
                edges.push_back(new_edge);
                is_boundary.push_back(1);
            }
        }
    }

    std::vector<std::pair<int,int>> boundary_edges(edges.size()-NbNotBoundary);
    int count=0;
    for (int i=0;i<edges.size();i++){
        if(is_boundary[i]){
            boundary_edges[count]=edges[i];
            count++;
        }
    }

    std::unordered_set<int> s;
    for (int i=0;i<boundary_edges.size();i++){
        s.insert(boundary_edges[i].first);
        s.insert(boundary_edges[i].second);
    }
    std::vector<int> boundary(s.size());
    boundary.assign( s.begin(), s.end() );
    sort( boundary.begin(), boundary.end() );
    std::vector<int> wo_boundary(nb_dofs-boundary.size());
    count=0;
    int count2=0;
    for (int i=0;i<wo_boundary.size();i++){
        while (count==boundary[count2]) {
            count++;
            count2++;
        }
        wo_boundary[i]=count;
        count++;
    }


    return wo_boundary;
}

}
#endif
