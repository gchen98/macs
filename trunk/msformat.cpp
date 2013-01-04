#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<sstream>

#include "constants.h"

using namespace std;


int main(int argc,char * argv[]){
    int totalpos;
    int totalpersons;
    ostringstream oss,oss_tree;
    srand(time(NULL));
    double r = rand();
    oss<<TEMPFILE_PREFIX<<r;
    oss_tree<<TREEFILE_PREFIX<<r;
    string tempfile=oss.str();
    ofstream haploDataFile(tempfile.data());
    string treefile=oss_tree.str();
    ofstream treeDataFile(treefile.data());

    double * positions = NULL;
    bool ** genotypes = NULL;
    bool slashprint = true;
    bool bPrintNewick = false;

    string line,constant;
    while(getline(cin,line)){
        istringstream tokens(line);
        tokens>>constant;
        //cout<<line<<endl;
        if (constant.compare(COMMAND)==0){
           string commandline;
           bool printed=false;
           while(tokens){
             tokens>>commandline;
             if (printed){ cout<<" "; }
             cout<<commandline;
             printed=true;
           }
           cout<<endl;
           // print the MS output header here
        }else if(constant.compare(SEED)==0){
           string seed;
           tokens>>seed;
           //getline(cin,line);
           // print the MS output header here
           cout<<seed<<endl;
        }else if(constant.compare(NEWICKTREE)==0){
           bPrintNewick = true;
           string tree;
           tokens>>tree;
           treeDataFile<<tree<<endl;
        }else if(constant.compare(MUTATIONSITE)==0){
            string index,position,mutationTime,mutations;
            tokens>>index>>position>>mutationTime>>mutations;
            haploDataFile<<index<<FIELD_DELIMITER<<position<<FIELD_DELIMITER
            <<mutations<<endl;
            //getline(cin,line);
            //while(line.compare(HAPLOEND)!=0){
            //    haploDataFile<<line<<endl;
            //    getline(cin,line);
            //}
        }else if (constant.compare(TOTALSAMPLES)==0){
           treeDataFile.close();
           haploDataFile.close();
           tokens>>totalpersons;
        }else if (constant.compare(TOTALSITES)==0){
           tokens>>totalpos;
        }else if (line.compare(SNPBEGIN)==0){
            if (slashprint){
             cout<<endl<<"//"<<endl;
             slashprint = false;
            }
            if (bPrintNewick){
              ifstream treeDataIn(treefile.data());
              string line;
              while(getline(treeDataIn,line)){
                cout<<line<<endl;
              }
              treeDataIn.close();
            }
            cout<<"segsites: "<<totalpos<<endl;
            cout<<"positions:";
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
        //getline(cin,line);
        else if (line.compare(SNPEND)==0){
            // print output in MS format now
            for (int i=0;i<totalpos;++i){
                cout<<" "<< setprecision(14)<<positions[i];
            }
            cout<<endl;
            // cleanup arrays
            delete [] positions;
            for (int j=0;j<totalpersons;++j){
                for (int i=0;i<totalpos;++i){
                    cout<<genotypes[j][i];
                }
                cout<<endl;
            }
            for (int i=0;i<totalpersons;++i){
              genotypes[i] = new bool[totalpos];
            }
            delete [] genotypes;
            // we can print the slashes for the next iteration
            slashprint = true;
            haploDataFile.open(tempfile.data());
            treeDataFile.open(treefile.data());
        }
    }
    haploDataFile.close();
    remove(tempfile.data());
    treeDataFile.close();
    remove(treefile.data());
}
