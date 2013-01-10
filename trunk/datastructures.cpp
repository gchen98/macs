#include<iostream>
#include<math.h>
#include<time.h>
#include<sstream>
#include "simulator.h"

using namespace std;

ChrPosition::ChrPosition(unsigned long int iGraphIteration,
double position){
    this->iGraphIteration = iGraphIteration;
    this->position = position;
}

PtrRefCountable::PtrRefCountable(){
    references = 0;
}

PtrRefCountable::~PtrRefCountable(){
//    cout<<"Parent destructor\n";
}

Mutation::Mutation(double dLocation,double dFreq){
    this->dLocation = dLocation;
    this->dFreq = dFreq;
    this->bPrintOutput = false;
}

Mutation::~Mutation(){
//    cout<<"Mutation destructor\n";
}

GeneConversion::GeneConversion(double dEndPos):
PtrRefCountable(){
    this ->dEndPos = dEndPos;
}

GeneConversion::~GeneConversion(){
    #ifdef DIAG
//    cout<<"Gene conversion destructor\n";
    #endif
}

HotSpotBin::HotSpotBin(double dStart,
    double dEnd,double dRate){
        this->dStart = dStart;
        this->dEnd = dEnd;
        this->dRate = dRate;
}

AlleleFreqBin::AlleleFreqBin(double dStart,
double dEnd, double dFreq){
    this->dStart = dStart;
    this->dEnd = dEnd;
    this->dFreq = dFreq;
    this->iObservedCounts = 0;
}

AlleleFreqBin::~AlleleFreqBin(){
    #ifdef DIAG
//    cout<<"AlleleFreqBin destructor\n";
    #endif
}


Population::Population(){
    setChrSampled(0);
    setPopSize(1);
    dGrowthAlpha=0;
    dLastTime=0;
}


void Population::setLastTime(double time){
    this->dLastTime = time;
    #ifdef DIAG
    if (time<0) throw "setLastTime, negative time";
    #endif
}

double Population::getLastTime(){
    return this->dLastTime;
}

void Population::setGrowthAlpha(double alpha){
    this->dGrowthAlpha = alpha;
    #ifdef DIAG
    if (this->dGrowthAlpha<0) throw "setGrowthAlpha, negative alpha";
    #endif
}

double Population::getGrowthAlpha(){
    return this->dGrowthAlpha;
}

void Population::setChrSampled(int iChrSampled){
    this->iChrSampled = iChrSampled;
    #ifdef DIAG
    if (this->iChrSampled<0) throw "setChrSampled, negative chrs";
    #endif
}

void Population::setPopSize(double dPopSize){
    this->dPopSize = dPopSize;
    #ifdef DIAG
    if (this->dPopSize<0) throw "setPopSize, negative pop";
    #endif
}

Edge::Edge(NodePtr & topNode,NodePtr & bottomNode):
PtrRefCountable(){
    #ifdef DIAG

    if (topNode->getHeight()<bottomNode->getHeight()){
        cerr<<"Error in creating edge.  Top node must be higher or equal\n";
        throw "Error in creating edge";
    }
    if (topNode->getType()==Node::MIGRATION &&
    bottomNode->getPopulation()==topNode->getPopulation()){
        cerr<<"Error in creating edge.  Top node must be different from "<<
        "migration node's population\n";
        throw "Error in creating edge";
    }
    cerr<<"edge of constructed.\n";
    #endif
    this->dLength = topNode->getHeight() - bottomNode->getHeight();
    this->topNode = topNode;
    this->getBottomNodeRef() = bottomNode;
    this->bDeleted = false;
    this->bInQueue = false;
    this->bInCurrentTree = false;
    this->iGraphIteration = 0;
}

Edge::~Edge(){
    #ifdef DIAG
//    cerr<<"edge of length: "<<dLength<<" destructed.\n";
    #endif

}

void Edge::setBottomNode(NodePtr & bottomNode){
    #ifdef DIAG
    if (topNode->getHeight()<bottomNode->getHeight())
    throw "Error in modifying edge.  Top node must be higher or equal";
    #endif
    this->dLength = topNode->getHeight() - bottomNode->getHeight();
    this->bottomNode = bottomNode;
}

const char* Node::getTypeStr(){
    switch(iType){
    case COAL:
        return "coal";
    case XOVER:
        return "xover";
    case MIGRATION:
        return "migr";
    case SAMPLE:
        return "sampl";
    case QUERY:
        return "query";
    }
    return "undef";
}

void Node::addNewEdge(EdgeLocation iLocation,
EdgePtr & newEdge){
    switch (iLocation){
        case Node::TOP_EDGE:
            if (topEdgeSize) this->topEdge2 = newEdge;
            else this->topEdge1 = newEdge;
            ++topEdgeSize;
            break;
        case Node::BOTTOM_EDGE:
            if (bottomEdgeSize) this->bottomEdge2 = newEdge;
            else this->bottomEdge1 = newEdge;
            ++bottomEdgeSize;
            break;
        break;
    }
    #ifdef DIAG
    switch(iType){
        case Node::COAL:
            if (bottomEdgeSize>2||
            topEdgeSize>1)
                throw "Coal node has too many edges\n";
            break;
        case Node::XOVER:
            if (topEdgeSize>2||
            bottomEdgeSize>1)
                throw "Xover node has too many edges\n";
            break;
        case Node::MIGRATION:
            if (bottomEdgeSize>1||
            topEdgeSize>1)
                throw "Migration node has too many edges\n";
            break;
        default:
            break;
    }
    #endif
}

void Node::replaceOldWithNewEdge(EdgeLocation iLocation,
EdgePtr & oldEdge,EdgePtr & newEdge){
    bool found = false;
    unsigned int i=0;
    if (iLocation==Node::TOP_EDGE){
        while(!found && i<this->topEdgeSize){
            WeakEdgePtr & topEdge = i?topEdge2:topEdge1;
            if (topEdge.lock()==oldEdge){
                topEdge=WeakEdgePtr(newEdge);
                found = true;
            }else{
                ++i;
            }
        }
        #ifdef DIAG
        if (!found) throw "Can't find top edge in replace edge";
        #endif
    }else if (iLocation==Node::BOTTOM_EDGE){
        while(!found && i<this->bottomEdgeSize){
            WeakEdgePtr & bottomEdge = i?bottomEdge2:bottomEdge1;
            if (bottomEdge.lock()==oldEdge){
                bottomEdge = WeakEdgePtr(newEdge);
                found = true;
            }else{
                ++i;
            }
        }
        if (!found) throw "Can't find bottom edge in replace edge";
    }
}


Node::Node(NodeType iType,short int iPopulation,double dHeight):
PtrRefCountable(){
    this->iType = iType;
    this->iPopulation = iPopulation;
    this->dHeight = dHeight;
    this->bEventDefined = false;
    this->topEdgeSize=0;
    this->bottomEdgeSize=0;
    this->bDeleted = false;
    #ifdef DIAG
    cerr<<"Graph node at "<<dHeight<<" constructed."<<endl;
    #endif
}

Node::~Node(){
    #ifdef DIAG
//    cerr<<"Graph node at "<<dHeight<<" destructed."<<endl;
    #endif
}

SampleNode::SampleNode(short int iPopulation,int iId):
Node(Node::SAMPLE,iPopulation,0.0){
    this->bAffected = false;
    this->iId = iId;
}


Event::~Event(){
    #ifdef DIAG
//    cerr<<"At time "<<this->dTime<<" Graph event destructor\n";
    #endif
}

Event::Event(EventType iType,double dTime):
PtrRefCountable(){
    #ifdef DIAG
    if (dTime<0) throw "Error in creating event. Time must be positive";
    cerr<<"At time "<<this->dTime<<" Graph event constructor\n";
    #endif
    this->iType = iType;
    this->dTime = dTime;
    this->bMarkedForDelete = false;
    if (dTime>Node::MAX_HEIGHT) throw "Stop here\n";
}



GenericEvent::GenericEvent(EventType iType,double dTime,
double dParameterValue):
Event(iType,dTime){
    this->dParameterValue = dParameterValue;
    //cerr<<"destroying generic event at height "<<dTime<<endl;
}

double GenericEvent::getParamValue(){
    return this->dParameterValue;
}

MigrationEvent::MigrationEvent(EventType iType,double dTime,
short int iPopMigratedTo,short int iPopMigratedFrom):
Event(iType,dTime){
    this->iPopMigratedTo = iPopMigratedTo;
    this->iPopMigratedFrom = iPopMigratedFrom;
}

MigrationRateEvent::MigrationRateEvent(EventType iType,
double dTime,short int iSourcePop,short int iDestPop,
    double dRate):
Event(iType,dTime){
    this->iSourcePop=iSourcePop;
    this->iDestPop=iDestPop;
    this->dRate=dRate;
}

short int MigrationRateEvent::getSourcePop(){
    return this->iSourcePop;
}

short int MigrationRateEvent::getDestPop(){
    return this->iDestPop;
}

double MigrationRateEvent::getRate(){
    return this->dRate;
}

MigrationRateMatrixEvent::MigrationRateMatrixEvent(
EventType iType,double dTime,MatrixDouble dMigrationMatrix):
Event(iType,dTime){
    this->dMigrationMatrix = dMigrationMatrix;
}

MatrixDouble MigrationRateMatrixEvent::getMigrationMatrix(){
    return this->dMigrationMatrix;
}

CoalEvent::CoalEvent(EventType iType,double dTime,
short int iPopulation):
Event(iType,dTime){
    this->iPopulation = iPopulation;
}


XoverEvent::XoverEvent(EventType iType,double dTime,
short int iPopulation):Event(iType,dTime){
    this->iPopulation = iPopulation;
}

PopSizeChangeEvent::PopSizeChangeEvent(EventType iType,double dTime,
short int iPopulationIndex,double dPopChangeParam):
Event(iType,dTime){

    this->iPopulationIndex=iPopulationIndex;
    this->dPopChangeParam=dPopChangeParam;
}

short int PopSizeChangeEvent::getPopulationIndex(){
    return iPopulationIndex;
}

double PopSizeChangeEvent::getPopChangeParam(){
    return dPopChangeParam;
}


PopJoinEvent::PopJoinEvent(EventType iType,double dTime,
short int iSourcePop,short int iDestPop):
Event(iType,dTime){
    this->iSourcePop = iSourcePop;
    this->iDestPop = iDestPop;
}

short int PopJoinEvent::getSourcePop(){
    return iSourcePop;
}

short int PopJoinEvent::getDestPop(){
    return iDestPop;
}

GraphBuilder::~GraphBuilder(){
    cerr<<"Graphbuilder destructor\n";
    this->pConfig = NULL;
    this->pRandNumGenerator = NULL;

    delete this->pEdgeListInARG;
    delete this->pEdgeVectorByPop;

    delete this->pVectorIndicesToRecycle;
    delete [] pTreeEdgesToCoalesceArray;
    delete this->pEdgeVectorInTree;
    delete [] this->pSampleNodeArray;
    MutationPtrVector::iterator it;
    for(it=pMutationPtrVector->begin();it!=pMutationPtrVector->end();++it){
        delete(*it);
    }
    delete this->pMutationPtrVector;
    delete this->pEventList;
    delete this->pGeneConversionPtrSet;
    delete [] sites;
    delete pChrPositionQueue;
}

GraphBuilder::GraphBuilder(Configuration *pConfig,RandNumGenerator * pRG){
    this->iGraphIteration = 0;
    this->bIncrementHistory = false;
    this->iTotalTreeEdges = 0;
    this->dArgLength = 0.;
    this->pConfig = pConfig;
    this->pChrPositionQueue = new ChrPositionQueue;
    this->dTrailingGap = pConfig->dBasesToTrack/pConfig->dSeqLength;
    this->bEndGeneConversion = false;
    this->bBeginGeneConversion = false;
    this->dScaledRecombRate = pConfig->dRecombRateRAcrossSites*
    (pConfig->dGeneConvRatio+1.);
    this->dMigrationMatrix = pConfig->dMigrationMatrix;
    this->pRandNumGenerator = pRG;
    this->pEdgeVectorByPop = new EdgePtrVectorByPop;
    this->pVectorIndicesToRecycle = new EdgeIndexQueueByPop;
    this->pEdgeListInARG = new EdgePtrList;
    this->pTreeEdgesToCoalesceArray = new EdgePtr[pConfig->iSampleSize];
    for (int i=0;i<pConfig->iTotalPops;++i){
        this->pEdgeVectorByPop->push_back(EdgePtrVector());
        this->pVectorIndicesToRecycle->push_back(EdgeIndexQueue());
    }
    this->pSampleNodeArray = new NodePtr[pConfig->iSampleSize];
    sites = new bool[pConfig->iSampleSize];
    this->pEdgeVectorInTree = new EdgePtrVector;
    this->pMutationPtrVector = new MutationPtrVector;
    this->pEventList = new EventPtrList;
    EventPtrList::iterator it = pConfig->pEventList->begin();
    for (it=pConfig->pEventList->begin();it!=pConfig->pEventList->end();++it){
        this->pEventList->push_back(*it);
    }
    this->pGeneConversionPtrSet = new GeneConversionPtrSet;
    if (pConfig->pAlleleFreqBinPtrSet!=NULL){
      AlleleFreqBinPtrSet::iterator it2;
      for (it2=pConfig->pAlleleFreqBinPtrSet->begin();it2!=pConfig->pAlleleFreqBinPtrSet->end();++it2){
        AlleleFreqBinPtr ptr = *it2;
        ptr->iObservedCounts = 0;
      }
    }
}

void GraphBuilder::checkPopCountIntegrity(PopVector & pPopList,
double dTime){
    int iTotalPops = pPopList.size();
    vector<int> iCounts;
    vector <vector <int> > edgeIndices;
    for (int k=0;k<iTotalPops;++k){
        vector<int> temp;
        edgeIndices.push_back(temp);
        iCounts.push_back(0);    }
    for (int k=0;k<iTotalPops;++k){
        EdgePtrVector & pEdgeVector = this->pEdgeVectorByPop->at(k);
        for (unsigned int i=0;i<pEdgeVector.size();++i){
            if (!pEdgeVector[i]->bDeleted && pEdgeVector[i]
            ->getBottomNodeRef()->getHeight()<=dTime &&
            pEdgeVector[i]->getTopNodeRef()->getHeight()>dTime){
                int iPopulation = pEdgeVector[i]->getBottomNodeRef()->getPopulation();
                ++iCounts[iPopulation];
                edgeIndices[iPopulation].push_back(i);
            }
        }
        if (coalescingEdge->getBottomNodeRef()->getPopulation()==k &&
            coalescingEdge->getBottomNodeRef()->getHeight()<=dTime &&
            coalescingEdge->getTopNodeRef()->getHeight()>dTime){
                ++iCounts[k];
        }
        if (originExtension->getBottomNodeRef()->getPopulation()==k &&
            originExtension->getBottomNodeRef()->getHeight()<=dTime &&
            originExtension->getTopNodeRef()->getHeight()>dTime){
                ++iCounts[k];
        }
    }

    for (int j=0;j<iTotalPops;++j){
        if (pPopList[j].getChrSampled()!=iCounts[j]){
            cerr<<"At time: "<<dTime<<endl;
            cerr<<"pop:"<<j<<",size:"<<pPopList[j].getChrSampled()<<endl;
            cerr<<"pop:"<<j<<",found:"<<iCounts[j]<<endl;
            printDataStructures();
            throw "Mismatch in pop counts in CheckPopCountIntegrity" ;
        }
    }
}

void GraphBuilder::insertNodeInRunningEdge(NodePtr & newNode,EdgePtr & tempEdge){

    NodePtr & bottomNode = tempEdge->getBottomNodeRef();
    NodePtr & topNode = tempEdge->getTopNodeRef();
    #ifdef DIAG
    if (newNode->getHeight()<=bottomNode->getHeight() ||
    newNode->getHeight() >= topNode->getHeight()){
           cerr<<"Edge has heights "<<bottomNode->getHeight()<<
        " and "<<topNode->getHeight()<<endl;
        cerr<<"New node has height "<<newNode->getHeight()<<endl;
        throw "Node to insert does not fit within coal edge ";
    }
    #endif

    EdgePtr tempEdgeCopy = tempEdge;

    tempEdge = EdgePtr(new Edge(topNode,newNode));
    tempEdge->iGraphIteration = iGraphIteration;
    newNode->addNewEdge(Node::TOP_EDGE,tempEdge);
    tempEdge->getTopNodeRef()->replaceOldWithNewEdge(
    Node::BOTTOM_EDGE,tempEdgeCopy,tempEdge);

    EdgePtr newBottomEdge = EdgePtr(new Edge(newNode,bottomNode));
    newBottomEdge->iGraphIteration = iGraphIteration;
    newNode->addNewEdge(Node::BOTTOM_EDGE,newBottomEdge);
    addEdge(newBottomEdge);
    bottomNode->replaceOldWithNewEdge(Node::TOP_EDGE,tempEdgeCopy,
    newBottomEdge);
}


void GraphBuilder::mergeEdges(EdgePtr & topEdge,EdgePtr & bottomEdge){
    #ifdef DIAG
    if (topEdge->getBottomNodeRef()->getHeight()!=bottomEdge->getTopNodeRef()->getHeight())
    throw "The top edge and bottom edge are mismatched for the join";
    #endif
    // mark expired node as deleted
    bottomEdge->getTopNodeRef()->bDeleted = true;
    NodePtr & bottomNode = bottomEdge->getBottomNodeRef();
    bottomNode->replaceOldWithNewEdge(Node::TOP_EDGE,bottomEdge,topEdge);
    topEdge->setBottomNode(bottomNode);
    deleteEdge(bottomEdge);
}

void GraphBuilder::insertNodeInEdge(NodePtr & newNode,
EdgePtr & selectedEdge){

    NodePtr bottomNodeCopy = selectedEdge->getBottomNodeRef();
    #ifdef DIAG
    NodePtr & topNode = selectedEdge->getTopNodeRef();
    if (newNode->getHeight()<bottomNodeCopy->getHeight() ||
    newNode->getHeight() > topNode->getHeight()){
        cerr<<"Edge has heights "<<bottomNodeCopy->getHeight()<<
        " and "<<topNode->getHeight()<<endl;
        cerr<<"New node has height "<<newNode->getHeight()<<endl;
        throw "Node to insert does not fit within edge";
    }
    #endif
    int iGraphIteration = selectedEdge->iGraphIteration;
    selectedEdge->setBottomNode(newNode);
    newNode->addNewEdge(Node::TOP_EDGE,selectedEdge);
    EdgePtr newBottomEdge = EdgePtr(new Edge(newNode,bottomNodeCopy));
    newBottomEdge->iGraphIteration = iGraphIteration;
    addEdge(newBottomEdge);
    bottomNodeCopy->replaceOldWithNewEdge(Node::TOP_EDGE,selectedEdge,
    newBottomEdge);
    newNode->addNewEdge(Node::BOTTOM_EDGE,newBottomEdge);
}



void GraphBuilder::deleteEdge(EdgePtr & edge){
    if (!edge->bDeleted){
        edge->bDeleted = true;
        if (pConfig->bDebug){
            cerr<<"Deleting edge with hts "<<edge->getBottomNodeRef()->getHeight()<<" and "<<
            edge->getTopNodeRef()->getHeight()<<endl;
        }
    }
}

// Insert into EdgeVector, pop refers to bottom node
void GraphBuilder::addEdge(EdgePtr & edge){
    unsigned int iPopulation = edge->getBottomNodeRef()->getPopulation();

    this->pEdgeListInARG->push_back(edge);

    if (iPopulation>=pEdgeVectorByPop->size()){
        if (iPopulation==pEdgeVectorByPop->size()){
            this->pEdgeVectorByPop->push_back(EdgePtrVector());
            this->pVectorIndicesToRecycle->push_back(EdgeIndexQueue());
        }else{
            cerr<<"Attempting to add edge of population "<<iPopulation<<" when the number of available population edge pools is only "<<pEdgeVectorByPop->size()<<". It is recommended that you increase the migration rates and/or number of sampled chromosomes."<<endl;
            throw "Data structure integrity error.";
        }
    }

    // Insert into the vector allowing for random access
    EdgePtrVector & pEdgeVector = this->pEdgeVectorByPop->
    at(iPopulation);
    EdgeIndexQueue & pVectorIndicesToRecycle = this->pVectorIndicesToRecycle->
    at(iPopulation);
    if (pVectorIndicesToRecycle.empty()){
        pEdgeVector.push_back(edge);
    }else{
        int iIndex = pVectorIndicesToRecycle.front();
        pVectorIndicesToRecycle.pop();
        pEdgeVector[iIndex] = edge;
    }
}

void GraphBuilder::addEdgeToCurrentTree(EdgePtr & edge){
    edge->bInCurrentTree = true;
    edge->iGraphIteration = this->iGraphIteration;
    if (iTotalTreeEdges<pEdgeVectorInTree->size()){
        pEdgeVectorInTree->at(iTotalTreeEdges) = edge;
    }else{
        pEdgeVectorInTree->push_back(edge);
    }
    ++iTotalTreeEdges;
}


void GraphBuilder::initializeCurrentTree(){
    // this method does two things, traverses the edges
    // to compute the total length, and clears out the status for
    // in current tree
    dLastTreeLength = 0.0;
    dArgLength = 0.;
    if (iGraphIteration==0){
        EdgePtrList::iterator it;
        for (it=pEdgeListInARG->begin();it!=pEdgeListInARG->end();++it){
            EdgePtr curEdge=*it;
            dLastTreeLength+=curEdge->getLength();
            this->addEdgeToCurrentTree(curEdge);
            curEdge->bInCurrentTree = false;
        }
        dArgLength = dLastTreeLength;
    }else{
        EdgePtrList::iterator it1=pEdgeListInARG->begin();
        while(it1!=pEdgeListInARG->end()){
            EdgePtr curEdge = *it1;
            if (!curEdge->bDeleted){
              dArgLength+=curEdge->getLength();
              if (curEdge->bInCurrentTree){
                dLastTreeLength+=curEdge->getLength();
                curEdge->bInCurrentTree = false;
              }
              ++it1;
            }else{
              it1 = pEdgeListInARG->erase(it1);
            }
        }
    }
}


void GraphBuilder::printHaplotypes(){
    unsigned int iTotalSites = pMutationPtrVector->size();
    bool bZeroCellCount=false;
    if (pConfig->bDebug) cerr<<"Total sites: "<<iTotalSites<<endl;
    if (iTotalSites){
        int iReducedSites=iTotalSites;
        if (pConfig->bSNPAscertainment){
            // first see if any expected count exceed actual counts
            bool bSufficientObs=false;
            do{
                bSufficientObs=true;
                AlleleFreqBinPtrSet::iterator it=pConfig->pAlleleFreqBinPtrSet->begin();
                while(bSufficientObs && !bZeroCellCount && it!=pConfig->pAlleleFreqBinPtrSet->end()){
                    AlleleFreqBinPtr bin = *it;
                    int iExpectedCount = static_cast<int>(bin->dFreq * iReducedSites);
                    if (!iExpectedCount && bin->dFreq>0.){
                        bZeroCellCount = true;
                        if (pConfig->bDebug){
                            cerr<<"Setting zero cell count true because expecting a proportion of "<<
                            bin->dFreq<<" but found "<<iExpectedCount<<endl;
                        }
                    }
                    else if (bin->iObservedCounts<iExpectedCount){
                        bSufficientObs = false;
                        --iReducedSites;
                    }else{
                        ++it;
                    }
                }
            }while(!bSufficientObs && !bZeroCellCount);
            if (bZeroCellCount){
                cerr<<"Warning: Some observed SNP counts were zero when they should have been positive.\n"<<
                "No ascertainment correction was applied.\n"<<
                "Try expanding frequency bin sizes and/or increasing mutation rate.\n";
                iReducedSites = 0;
            }else{
                int tally=0;
                for (AlleleFreqBinPtrSet::iterator it=pConfig->pAlleleFreqBinPtrSet->begin();
                it!=pConfig->pAlleleFreqBinPtrSet->end();++it){
                    AlleleFreqBinPtr bin = *it;
                    double dStart = bin->dStart;
                    double dEnd = bin->dEnd;
                    int iExpectedCount = static_cast<int>(bin->dFreq * iReducedSites);

                    tally+=iExpectedCount;
                    if (pConfig->bDebug) {
                        cerr<<"Looking for SNPs in range "<<dStart<<" to "<<dEnd<<endl;
                        cerr<<"iObserved count: "<<bin->iObservedCounts<<endl;
                        cerr<<"Expecting "<<iExpectedCount<<" SNPS."<<endl;
                    }
                    #ifdef DIAG
                    if (bin->iObservedCounts<iExpectedCount){
                        throw "Too many expected counts";
                    }
                    #endif
                    while(iExpectedCount>0){
                        int iRandIndex = static_cast<int>(pRandNumGenerator->unifRV()*iTotalSites);
                        MutationPtr mutation = pMutationPtrVector->at(iRandIndex);
                        if (!mutation->bPrintOutput && mutation->dFreq>=dStart && mutation->dFreq<=dEnd){
                            mutation->bPrintOutput = true;
                            --iExpectedCount;
                        }
                    }
                }
                iReducedSites = tally;
                cerr<<"Total sites reduced from "<<iTotalSites<<" to "<<iReducedSites<<endl;
            }
        }
        if (iReducedSites){
            MutationPtrVector::iterator it;
            // copy to a temporary vector if ascertained
            cout<<TOTALSAMPLES<<FIELD_DELIMITER<<pConfig->iSampleSize<<endl;
            cout<<TOTALSITES<<FIELD_DELIMITER<<iReducedSites<<endl;
            cout<<SNPBEGIN<<endl;
            if (pConfig->bSNPAscertainment && !bZeroCellCount){
                int origIndex=0;
                bool indexPrinted=false;
                for (it = pMutationPtrVector->begin();
                it!=pMutationPtrVector->end();++it){
                    MutationPtr mutation = *it;
                    if (mutation->bPrintOutput){
                        if (indexPrinted) cout<<FIELD_DELIMITER;
                        cout<<origIndex;
                        indexPrinted=true;
                    }
                    ++origIndex;
                }
            }else{
                for(int i=0;i<iReducedSites;++i){
                  if (i) cout<<FIELD_DELIMITER;
                  cout<<i;
                }
            }
            cout<<endl<<SNPEND<<endl;
        }
    }
}


void GraphBuilder::printDataStructures(){
    cerr<<endl<<"*** Begin printing structures ***"<<endl;

    double trueLen = 0.0;
    cerr<<"Full ARG (list of edges)\n";
    trueLen = 0.0;
    for (EdgePtrList::iterator it=pEdgeListInARG->begin();it!=pEdgeListInARG->end();it++){
        EdgePtr curEdge = *it;
        cerr<<"low:ht:"<<curEdge->getBottomNodeRef()->getHeight()<<
        ",type:"<<curEdge->getBottomNodeRef()->getTypeStr()<<
        ",pop:"<<curEdge->getBottomNodeRef()->getPopulation()<<
        ";high:ht:"<<curEdge->getTopNodeRef()->getHeight()<<
        ",type:"<<curEdge->getTopNodeRef()->getTypeStr()<<
        ",pop:"<<curEdge->getTopNodeRef()->getPopulation()<<
        ",del:"<<curEdge->bDeleted<<
//        ",length:"<<curEdge->getLength()<<
        ";hist:"<<curEdge->iGraphIteration<<endl;
        trueLen+=curEdge->getLength();
    }


    cerr<<"Last tree (list of edges)\n";
//    it;
    trueLen = 0.0;
    EdgePtrVector::iterator it=pEdgeVectorInTree->begin();
    unsigned int count=0;
    while(count<iTotalTreeEdges){
//    for (EdgePtrVector::iterator it=pEdgeVectorInTree->begin();it!=pEdgeVectorInTree->end();++it){
        EdgePtr curEdge = *it;
        cerr<<"low_ht:"<<curEdge->getBottomNodeRef()->getHeight()<<
        ",type:"<<curEdge->getBottomNodeRef()->getTypeStr()<<
        ",pop:"<<curEdge->getBottomNodeRef()->getPopulation()<<
        ";high_ht:"<<curEdge->getTopNodeRef()->getHeight()<<
        ",type:"<<curEdge->getTopNodeRef()->getTypeStr()<<
        ",pop:"<<curEdge->getTopNodeRef()->getPopulation()<<endl;
//        ",intree:"<<curEdge->bInCurrentTree<<endl;
//        ",length:"<<curEdge->getLength()<<
//        ";hist:"<<curEdge->iGraphIteration<<endl;
        trueLen+=curEdge->getLength();
        ++count;
        ++it;
    }
    cerr<<"MRCA: "<<localMRCA->getHeight()<<endl;
    cerr<<"Graph grandMRCA: "<<grandMRCA->getHeight()<<endl;
//    cerr<<"true length:"<<trueLen<<endl;

    cerr<<"*** Done printing structures ***"<<endl;
}

string GraphBuilder::getNewickTree(double lastCoalHeight,NodePtr & curNode){
    ostringstream oss;
    if (curNode->getType()==Node::SAMPLE){
        // print the Newick format with the Sample ID followed by length
        SampleNode * pSampleNode = static_cast<SampleNode *>(curNode.get());
        oss<<pSampleNode->iId<<":"<<lastCoalHeight-curNode->getHeight();
    }else if (curNode->getType()==Node::COAL){
        // if the two lower branches are in the current tree, then print the
        // Newick format, and recursively repeat
        if (iGraphIteration==curNode->getBottomEdgeByIndex(0)->iGraphIteration &&
        iGraphIteration==curNode->getBottomEdgeByIndex(1)->iGraphIteration){
            oss<<"("<<getNewickTree(curNode->getHeight(),
            curNode->getBottomEdgeByIndex(0)->getBottomNodeRef())
            <<","<<getNewickTree(curNode->getHeight(),
            curNode->getBottomEdgeByIndex(1)->getBottomNodeRef())
            <<")";
            // if we are at the MRCA we don't need to output the branch length since
            // it is obviously 0 length
            if (curNode!=localMRCA) oss<<":"<<lastCoalHeight-curNode->getHeight();
//            cout<<"got here\n";
        }else{
            for(int i=0;i<2;++i){
                if (curNode->getBottomEdgeByIndex(i)->iGraphIteration==iGraphIteration){
                    oss<<getNewickTree(lastCoalHeight,
                    curNode->getBottomEdgeByIndex(i)->getBottomNodeRef());
                }
            }
        }
    }else {
        // this is for everything else(e.g. migration, xover nodes),
        // don't print anything, but recursively
        // descend down the graph
        oss<<getNewickTree(lastCoalHeight,
        curNode->getBottomEdgeByIndex(0)->getBottomNodeRef());
    }
    return oss.str();
}

