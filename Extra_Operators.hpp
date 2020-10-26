#include "VectorND.hpp"

using namespace std;

bool operator == (VectorND& u, VectorND& v){
    bool Entry_Equality;
    if( u.dim()==v.dim() ){
        Entry_Equality = true;
        int dim = u.dim();
        for(int index=0; index<dim; ++index)
            Entry_Equality = Entry_Equality && (u[index]==v[index]);
    }
    else
        cout << "Problemas de dimensiones en operador ==" << endl;
    return Entry_Equality;
}

bool operator != (VectorND& u, VectorND& v){
    bool Entry_Equality;
    if( u.dim()==v.dim() ){
        Entry_Equality = true;
        int dim = u.dim();
        for(int index=0; index<dim; ++index)
            Entry_Equality = Entry_Equality && (u[index]!=v[index]);
    }
    else
        cout << "Problemas de dimensiones en operador ==" << endl;
    return ~Entry_Equality;
}

ostream& operator << (ostream& os, VectorND& vec){
    os << "<";
    for(int index=0; index<vec.dim(); ++index)
        ( index!=(vec.dim()-1) ) ? os << vec[index] << "," : os << vec[index] << ">";
    return os;
}
