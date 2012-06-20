#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>

#include "constants.h"

using namespace std;

int main(int argc,char * argv[]){
    int totalpos;
    int totalpersons;
    ostringstream oss;
    srand(time(NULL));
    oss<<TEMPFILE_PREFIX<<rand();
    string tempfile=oss.str();

    double * positions;
    bool ** genotypes;

    string line;
    while(getline(cin,line)){
        if (line.compare(HAPLOBEGIN)==0){
            ofstream haploDataFile(tempfile.data());
            getline(cin,line);
            while(line.compare(HAPLOEND)!=0){
                haploDataFile<<line<<endl;
                getline(cin,line);
            }
            haploDataFile.close();
        }
        for(int i=0;i<2;++i){
          getline(cin,line);
          istringstream info(line);
          string totalstr;
          int total;
          info>>totalstr>>total;
          if (totalstr.compare(TOTALSAMPLES)==0){
             totalpersons = total;
          }else if(totalstr.compare(TOTALSITES)==0){
             totalpos=total;
          }
        }
        getline(cin,line);
        if (line.compare(SNPBEGIN)==0){
            positions = new double[totalpos];
            genotypes = new bool * [totalpersons];
            for (int i=0;i<totalpersons;++i){
                genotypes[i] = new bool[totalpos];
            }
            ifstream haploDataFile(tempfile.data());
            // read in the positions
            getline(cin,line);
            istringstream input(line);
            int cur_snp_index = -1;
            double cur_snp_pos = 0.;
            string cur_haplo;
            for (int i=0;i<totalpos;++i){
                int selectedsnp;
                input>>selectedsnp;
                string haploline;
                do{
                    getline(haploDataFile,haploline);
                    istringstream haploInput(haploline);
                    haploInput>>cur_snp_index>>cur_snp_pos>>cur_haplo;
                }while(selectedsnp!=cur_snp_index);
                positions[i] = cur_snp_pos;
                for (int j=0;j<totalpersons;++j){
                    char cur_char = cur_haplo[j];
                    switch(cur_char){
                        case '0':
                            genotypes[j][i] = 0;
                            break;
                        case '1':
                            genotypes[j][i] = 1;
                            break;
                    }
                }
            }
            haploDataFile.close();
            remove(tempfile.data());
        }
        getline(cin,line);
        if (line.compare(SNPEND)==0){
            // cleanup arrays
            delete [] positions;
            bool printed=false;
            for (int j=0;j<totalpersons;++j){
                for (int i=0;i<totalpos;++i){
                    if (printed){
                       cout<<"\t";
                    }
                    cout<<genotypes[j][i];
                    printed=true;
                }
                cout<<endl;
            }
            for (int i=0;i<totalpersons;++i){
              delete genotypes[i];
            }
            delete [] genotypes;
        }
    }
}
